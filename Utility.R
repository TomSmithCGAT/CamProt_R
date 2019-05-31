# Utility functions useful across all analyses
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(grid))
suppressMessages(library(VennDiagram))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(MSnbase))
suppressMessages(library(robustbase))

my_theme <- theme_bw() + theme(
  text=element_text(size=20,  family="serif"),
  panel.grid=element_blank(), aspect.ratio=1,
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# ----------------------------------------------------------------------------------------------------------------------
# Function	: features
# Aim		: To aggregate features into master-proteins and carry out QC on dataset. 
# Input :
#    infile = file containing output from mascot server
#    silac = input is from a SILAC experiment
# Output  	: A dataframe with features associated to a unique master.protein, which is not cRAP or cRAP-associated. 
# ----------------------------------------------------------------------------------------------------------------------
print_n_feature <- function(features_df, message){
  cat(sprintf("%s\t%s\n", length(rownames(features_df)), message))
}

print_n_prot <- function(features_df){
  cat(sprintf("These features are associated with %s master proteins\n",
              length(unique(features_df$master_protein))))}

print_summaries <- function(features_df, message){
  print_n_feature(features_df, message)
  print_n_prot(features_df)
}

parse_features <- function(infile, unique_master=TRUE, silac=FALSE, TMT=FALSE,
                           level="peptide", filter_crap=TRUE){
  
  if(!level %in% c("PSM", "peptide")){
    stop("level must be PSM or peptide")
  }
  
  features_df <- read.delim(infile,  header=T, stringsAsFactors=FALSE)
  cat("Tally of features at each stage:\n")
  
  print_summaries(features_df, "All features with a PSM")
  
  features_df <- features_df %>% filter(master_protein!="")
  print_summaries(features_df, "Excluding features without a master protein")
  
  if(unique_master){
    features_df <- features_df %>% filter(unique==1)
    print_summaries(features_df, "Excluding features without a unique master protein")
  }
  
  if(filter_crap){
    features_df <- features_df %>% filter(crap_protein==0)
    print_summaries(features_df, "Excluding features matching a cRAP protein")
    
    features_df <- features_df %>% filter(associated_crap_protein==0)
    print_summaries(features_df, "Excluding features associated with a cRAP protein")
  }

  if(silac|TMT & level=="peptide"){
    features_df <- features_df %>% filter(Quan.Info!="")
    print_summaries(features_df, "Excluding features without quantification")
  }
  return(features_df)
}

removeCrap <- function(obj, protein_col="Protein.Accessions"){
  cat(sprintf("Input data: %s rows\n", nrow(obj)))
  obj <- obj %>% filter(!grepl("cRAP", !!as.symbol(protein_col), ignore.case=FALSE))
  cat(sprintf("Output data: %s rows\n", nrow(obj)))
  return(obj)
}

summariseMissing <- function(res){
  cat("\ntallies for missing data (# samples with missing)")
  print(table(rowSums(is.na(as.data.frame(exprs(res))))))
}

makeMSNSet <- function(obj, samples_inf, ab_col_ix=3, level="peptide", quant_name="Abundance"){
  # make dataframes for MSnset object
  rownames(obj) <- seq(1, length(obj[,1]))
  
  meta_columns <- colnames(obj)
  meta_columns <- meta_columns[grep("Found.*", meta_columns, invert=TRUE)]
  meta_columns <- meta_columns[grep(sprintf('%s.*', quant_name), meta_columns, invert=TRUE)]


  if(level=="PSM"){
    abundance_columns <- colnames(obj)[grep(sprintf('%s.*', quant_name), colnames(obj))]
    renamed_abundance_columns <- sapply(strsplit(abundance_columns, "\\."), "[[", ab_col_ix)
    }
  else if(level=="peptide"){
    abundance_columns <- colnames(obj)[grep(sprintf('%s.*', quant_name), colnames(obj))]
    abundance_columns <- abundance_columns[grep(sprintf('%s.Count', quant_name), abundance_columns, invert=TRUE)]
    renamed_abundance_columns <- sapply(strsplit(abundance_columns, "\\."), "[[", ab_col_ix)
  }
  else{
    stop("level must be PSM or peptide")
  }
  
  exprsCsv <- obj[,abundance_columns]
  colnames(exprsCsv) <- renamed_abundance_columns
  
  exprsCsv[exprsCsv==""] <- NA
  
  fdataCsv <- obj[,meta_columns]
  
  pdataCsv <- read.table(samples_inf, sep="\t", header=T, row.names = 1, colClasses="character")
  
  exprsCsv <- exprsCsv[,rownames(pdataCsv)]
  
  res <- MSnSet(as.matrix(sapply(exprsCsv, as.numeric)), fdataCsv, pdataCsv)
  
  summariseMissing(res)
  
  cat(sprintf("\n%s peptides do not have any quantification values\n\n", sum(rowSums(is.na(exprs(res)))==ncol(exprs(res)))))
  res <- res[rowSums(is.na(exprs(res)))!=ncol(exprs(res)),] # exclude peptides without any quantification
  return(res)
}

plotMissing <- function(obj, verbose=TRUE, ...){
  tmp_obj <- obj
  exprs(tmp_obj)[!is.na(exprs(tmp_obj))] <- 1
  exprs(tmp_obj)[is.na(exprs(tmp_obj))] <- 0
  
  missing <- exprs(tmp_obj)
  missing <- missing[rowSums(missing==0)>0,] # identify features without missing values
  
  if(verbose){
    cat(sprintf("Out of %s total features, %s (%s%%) have missing values\n",
                length(rownames(exprs(tmp_obj))), length(rownames(missing)),
                round(100*length(rownames(missing))/length(rownames(exprs(tmp_obj))),3)))
    
    print(table(rowSums(missing==0)))
  }
  
  if(length(rownames(missing))>0){
    missing_twice <- missing[rowSums(missing==0)>1,]
    
    if(verbose){
      cat(sprintf("And %s features have more than one missing value\n", length(rownames(missing_twice))))
    }
    
    colnames(missing) <- pData(tmp_obj)$Sample_name
    
    heatmap.2(missing, col = c("lightgray", "black"),
              scale = "none", dendrogram = "none",
              trace = "none", keysize = 0.5, key = FALSE,Colv = F, labRow=F,
              #RowSideColors = ifelse(fData(x)$randna, "orange", "brown"),
              ...)
  }
}


plotLabelQuant <- function(obj, log=F, print=TRUE){
  tmp_df <- data.frame(exprs(obj))
  colnames(tmp_df) <- pData(obj)$Sample_name
  tmp_df[tmp_df==""] <- NA
  tmp_df <- melt(tmp_df, id.vars=NULL)
  tmp_df$value <- as.numeric(as.character(tmp_df$value))
  
  if(log){
    tmp_df$value = log(tmp_df$value,2)  
  }
  
  p <- ggplot(tmp_df) + my_theme
  
  p1 <- p + geom_boxplot(aes(variable, value)) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
    ylab("Peptide intensity (log2) ") + xlab("") +
    scale_y_continuous(breaks=seq(-100,100,2))
  
  p2 <- p + geom_density(aes(value, col=variable)) +
    xlab("Peptide intensity (log2) ") + ylab("Density")
  
  if(print){
    print(p1)
    print(p2)
  }
  
  return(list("p1"=p1, "p2"=p2))
}

myAggFunction <- function(x, FUN=sum){
  if(sum(is.finite(x))==0){
    return(NA)
  }
  else{
    return(FUN(x[is.finite(x)]))
  }
}

agg_to_peptide_mod <- function(obj, gb=NULL, mod_col="Modifications", seq_col="Sequence", fun=sum){
  
  if(missing(gb)){
    gb <- paste(fData(obj)[[seq_col]], fData(obj)[[mod_col]])
  }
  # we'll use our own aggregation function here to obtain the sum over all the unique peptide sequence + modification
  
  pep_agg <- combineFeatures(obj, groupBy=gb, fun=function(x) myAggFunction(x, FUN=fun))
  pData(pep_agg) <- pData(obj)
  return(pep_agg)
}


agg_to_peptides <- function(obj, gb=NULL, seq_col="Sequence", fun=sum){
  
  # remove any "CV.Abundance" columns. Otherwise, combineFeatures will throw an error
  fData(obj) <- fData(obj)[,grep("CV.*", colnames(fData(obj)), invert=TRUE)]
  
  if(missing(gb)){
    gb <- fData(obj)[[seq_col]]
  }
  
  # we'll use our own aggregation function here to obtain the sum over all the modified peptides
  pep_agg <- combineFeatures(obj, groupBy=gb, fun=function(x) myAggFunction(x, FUN=fun))
  pData(pep_agg) <- pData(obj)
  return(pep_agg)
}


agg_to_protein <- function(obj, protein_col="Master.Protein.Accessions", gb=NULL, fun=median){
  
  if(missing(gb)){
    gb <- fData(obj)[[protein_col]]
  }
  
  # remove any "CV.Abundance" columns. Otherwise, combineFeatures will throw an error
  fData(obj) <- fData(obj)[,grep("CV.*", colnames(fData(obj)), invert=TRUE)]
  
  # we'll use our own aggregation function here to obtain the median over all the peptides
  pep_agg <- combineFeatures(obj, groupBy=gb, fun=function(x) myAggFunction(x, FUN=fun))
  pData(pep_agg) <- pData(obj)
  return(pep_agg)
}


getMedians <- function(obj){
  medians <- colMedians(exprs(obj), na.rm=TRUE)
  return(medians)
}


centerNormalise <- function(obj, medians=NULL){

  if(missing(medians)){
    medians <- getMedians(obj)
    medians <- medians/mean(medians)
  }
  else{
    observed_medians <- colMedians(exprs(obj), na.rm=TRUE)
    medians = medians / (mean(medians)/mean(observed_medians))
  }
  
  exprs(obj) <- t(t(exprs(obj))/medians)
  
  return(obj)
}


logCenterNormalise <- function(obj, medians=NULL){
  exprs(obj) <- log2(exprs(obj))
  
  if(missing(medians)){
    medians <- getMedians(obj)
  }
  else{
    observed_medians <- colMedians(exprs(obj), na.rm=TRUE)
    medians = medians / (mean(medians)/mean(observed_medians))
  }
  
  exprs(obj) <- t(t(exprs(obj))-medians)
  
  return(obj)
}

makeVennPlot <- function(set1, set2, set3,
                         set1_name="Identified proteins",
                         set2_name="RBPs",
                         set3_name="glycoproteins",
                         subset_to=NULL){
  if(!missing(subset_to)){
    set1 <- intersect(set1, subset_to)
    set2 <- intersect(set2, subset_to)
    set3 <- intersect(set3, subset_to)
  }
  
  grid.newpage()
  
  venn.plot <- draw.triple.venn(
    area1 = length(set1),
    area2 = length(set2),
    area3 = length(set3),
    n12 = length(intersect(set1, set2)),
    n13 = length(intersect(set1, set3)),
    n23 = length(intersect(set2, set3)),
    n123 = length(intersect(intersect(set1, set2), set3)),
    category = c(set1_name, set2_name, set3_name),
    fill = c("steelblue2", "orangered2", "green2"),
    alpha =.2)
  
  grid.draw(venn.plot)
}


runFishersEnrichment <- function(df, col, identified_col="identified"){
  
  print(table(df[[col]], df[[identified_col]]))
  print(fisher.test(table(df[[col]], df[[identified_col]])))
  
}


harm_mean <- function(values){
  return(1/mean(1/values))
}

normrank <- function(array){
  ra <- rank(array)
  return(ra/max(ra))
}


plotMultipleDistributions <- function(df, gene_sets, ylb="Abundance (log2)"){
  
  df <- melt(as.matrix(df))
  df <- separate(df, Var2, into=c("dosage", "replicate"), sep="_", remove=FALSE)
  
  tmp_dfs <- NULL
  for(cat in names(gene_sets)){
    genes = gene_sets[[cat]]
    tmp_df <- df[df$Var1 %in% genes,]
    intersect_genes <- (sum(df$Var1 %in% genes))/length(unique(df$Var2))
    lab <- sprintf("%s \n(%s)", cat, intersect_genes)
    tmp_df$lab <- lab
    tmp_dfs[[lab]] <- tmp_df
    
  }
  
  final_df <- do.call("rbind", tmp_dfs)
  final_df$lab <- factor(final_df$lab, levels=names(tmp_dfs))
  
  p <- ggplot(aggregate(value~Var1+dosage+lab, final_df, FUN=mean),
              aes(dosage, value, colour=lab)) +
    my_theme +
    scale_colour_discrete(guide=F) +
    xlab("")+ ylab(ylb) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
    facet_wrap(~lab)
  
  p2 <- p + geom_violin() + geom_boxplot(notch=TRUE, outlier.size=NA, width=0.1)
  suppressMessages(print(p2))
  suppressMessages(print(p2 %+% final_df + aes(x=Var2)))
  
  p3 <- p +
    aes(fill=lab) +
    scale_fill_discrete(guide=F) +
    stat_summary(geom="pointrange", colour="black", size=0.5)#width=0.3, lwd=1, size=1)
  print(p3)
  
  p4 <- p + stat_summary(aes(group=Var1), fun.y = mean, geom="line", alpha=0.1)
  print(p4)
  
}

writeListToFile <- function(obj, file){
  sink(file)
  writeLines(unlist(lapply(obj, paste, collapse=" ")))
  sink()
}

readListFromFile <- function(file){
  return(scan(file, what="", sep="\n"))
}


# Extract quantification info, with optional feature columns too
extractQuantData <- function(obj, feature_columns=NULL){
  
  exprs_df <- data.frame(exprs(obj))
  
  if(!missing(feature_columns)){
    colnames(exprs_df) <- pData(obj)$Sample_name
    exprs_df <- cbind(fData(obj)[feature_columns], exprs_df)
  }
  
  return(exprs_df)
}


four.venn <- function(a, b, c, d, ...){
  
  area1 = length(a)
  area2 = length(b)
  area3 = length(c)
  area4 = length(d)
  n12 = (intersect(a, b))
  n13 = (intersect(a,c))
  n14 = (intersect(a,d))
  n23 = (intersect(b,c))
  n24 = (intersect(b,d))
  n34 = (intersect(c,d))
  n123 = (intersect(n12, c))
  n124 = intersect(n12, d)
  n134 = intersect(n13, d)
  n234 = intersect(n23, d)
  n1234 = intersect(n123, d)
  
  draw.quad.venn(area1, area2, area3, area4, length(n12), length(n13), length(n14), 
                      length(n23), length(n24), length(n34), length(n123), 
                      length(n124), length(n134), length(n234), length(n1234), 
                      fill = cbPalette[2:5], 
                      alpha = rep(0.4, 4), 
                      lty = "blank", 
                      fontface = "bold",
                      fontfamily = "sans", 
                      cat.fontface = "bold",
                      cat.default.pos = "outer",
                      cat.fontfamily = "sans", ...)
}

singleLM <- function(df, dep_var="value", indep_var="stage"){
  model <- paste0(dep_var, " ~ ", indep_var)
  fit <- lm(formula=formula(model), data=df)
  
  p.value <- summary(fit)$coefficients[,4]
  p.value <- p.value[2:length(p.value)]
  #names(p.value) <- paste0("p.value_", names(p.value))
  
  coeffs <- coef(fit)
  coeffs <- coeffs[2:length(coeffs)]
  
  tmp_lm_df <- data.frame(cbind(coeffs, p.value))
  tmp_lm_df$name <- names(coeffs)
  
  return(tmp_lm_df)
}

plotP <- function(lm_df){
  
  p <- ggplot(lm_df, aes(p.value)) + geom_histogram(bins=100) + my_theme
  print(p)
}

runLM <- function(df, dep_var, indep_var, abundance_df, protein_col="Var1", abundance_merge_col="uniprot_ac"){
  
  lm_df <- ddply(df, protein_col, function(x) singleLM(x, dep_var, indep_var))
  
  plotP(lm_df)
  print(head(lm_df))
  lm_df$BH <- p.adjust(lm_df$p.value, method="BH")
  lm_df$sig <- lm_df$BH<0.01
  lm_df$sig_up <- (lm_df$sig & lm_df$coeffs > 0)
  lm_df$sig_dw <- (lm_df$sig & lm_df$coeffs < 0)
  
  print("Tallies for significant changes:")
  print(table(lm_df$sig))
  
  print("Tallies for significant changes(unique proteins):")
  print(table(unique(lm_df[[protein_col]]) %in% unique(lm_df[lm_df$sig==TRUE, protein_col])))
  
  means <- aggregate(df$value, by=list(df[[protein_col]]), FUN=mean)
  colnames(means) <- c(protein_col, "mean_abundance")
  
  lm_df <- merge(lm_df, means, by=protein_col)
  
  lm_df <- merge(lm_df, abundance_df, by.x=protein_col, by.y=abundance_merge_col)
  
  return(lm_df)
}

# ------------------------
# Function	: shortenTerm
# ------------------------
shortenTerm <- function(term, cut_at){
  
  if(is.na(term)){
    return(NA)
  }
  if(nchar(term)>cut_at){
    return(paste0(substr(term, 1, cut_at), " [..]"))
  }
  else{
    return(term)
  }
}

# ------------------------
# Function	: shortenTerm
# ------------------------
add_newlines <- function(input_string, max_length=20,sep=" "){
  if(is.na(input_string)){
    return(input_string)
  }
  input_string <- gsub(",", "'", input_string)
  
  if (nchar(input_string)<=max_length){
    return(input_string)
  }
  else{
    components <- strsplit(input_string, sep)
    if (any(sapply(components, FUN=function(x) nchar(x)>max_length))){
      stop(sprintf("impossible to break down below max_length, require at least %i characters", max(sapply(components ,FUN=nchar))))
    }
    else{
      output_string = ""
      char_since_newline = 0
      for (component in components[[1]]){
        if(char_since_newline + nchar(component)<=max_length){
          output_string <- paste0(output_string, sep, component)
          char_since_newline <- char_since_newline + nchar(component)
        }
        else{
          output_string <- paste0(output_string, "\n", component)
          char_since_newline <- nchar(component)
        }
      }
      return(substring(output_string, 2))
    }
  }
}
