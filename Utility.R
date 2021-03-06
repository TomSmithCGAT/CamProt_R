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
suppressMessages(library(biobroom))

my_theme_box <- theme_bw() + theme(
  text=element_text(size=20,  family="serif"),
  panel.grid=element_blank(), aspect.ratio=1,
  panel.background = element_blank())

my_theme_no_box <- my_theme_box + theme(
  panel.border = element_blank(),
  axis.line = element_line(colour = "black"))

my_theme <- my_theme_no_box

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



print_n_feature <- function(features_df, message){
  cat(sprintf("%s\t%s\n", length(rownames(features_df)), message))
}

print_n_prot <- function(features_df, master_protein_col){
  cat(sprintf("These features are associated with %s master proteins\n",
              length(unique(features_df[[master_protein_col]]))))}

print_summaries <- function(features_df, master_protein_col, message){
  print_n_feature(features_df, message)
  print_n_prot(features_df, master_protein_col)
}

removeCrap <- function(obj, protein_col="Protein.Accessions"){
  cat(sprintf("Input data: %s rows\n", nrow(obj)))
  obj <- obj %>% filter(!grepl("cRAP", !!as.symbol(protein_col), ignore.case=FALSE))
  cat(sprintf("Output data: %s rows\n", nrow(obj)))
  return(obj)
}

# Steps are:
# 1. Filter modifications
# 2. Center-median normalise PSM
# 3. Aggregate PSM level data to unique peptide sequence + modification
# 4. Sum normalise
# 5. Impute missing values
# 6. (Optional) Aggregate to unique peptide sequence, sum normalise and impute missing values (starts from step 4)
# 7. (Optional) Aggregate to unique protein ID, sum normalise and impute missing values (starts from step 6 prior to sum normalisation)
#
# Input :
#   raw_psm: [REQUIRED] The PSM level dataframe. This can be obtained by parsing PD output with `parse_features`. If PSMs contain PTMs,
#             you may want to further parse with `parsePTMScores`, `addPTMPositions` and `addSiteSequence` first.
#   sample_infile: [REQUIRED] Filename for table mapping TMT tags to sample names. e.g:
#                   Tag     Sample_name
#                   126     100
#                   127N    400
#                   [...]   [...]

# ----------------------------------------------------------------------------------------------------------------------
# Function	: parse_features
# Aim		: Parse output from PD and filter
#
# Steps are:
# 1. Read in PD data
# 2. Exclude features without a master protein
# 3. (Optional) Exclude features without a unique master protein (Number.of.Protein.Groups==1)
# 4. (Optional) Exclude features matching a cRAP protein
# 5. (Optional) Exclude features matching a proteins associated with a cRAP protein (see below)
# 6. Filter out features without quantification values (only if (TMT or SILAC) & peptide level input)
#
# "Associated cRAP"
# I've observed that the cRAP database does not contain all the possible cRAP proteins.
# E.g some features can be assigned to a keratin which is not in the cRAP database.
# In most cases, these "associated cRAP" proteins will have at least one feature shared with a cRAP
# protein. Thus, use this to remove them: see simplified example below
#
# feature Protein.Accessions          Master.Protein.Accessions
# f1      protein1, protein2, cRAP1   protein 1
# f2      protein1, protein3          protein 3
# f3      protein2                    protein 2
#
# Here, f1 indicates that protein1 and protein2 are associated with a cRAP protein.
# f2 and f3 are therefore filtered out as well, regardless of the Master.Protein.Accession column value
#
#
# Input:
#    infile: [REQUIRED] file containing output from Proteome Discoverer
#    master_protein_col: [REQUIRED; DEFUALT="Master.Protein.Accessions"]
#                        The column containing the master protein.
#    protein_col: [REQUIRED; DEFUALT="Protein.Accessions"]
#                 The column containing all the protein matches.
#    unique_master: [DEFAULT=TRUE] Filter out features without a unique master protein
#    silac: [DEFAULT=FALSE] SILAC experiment
#    TMT: [DEFAULT=FALSE] TMT experiment
#    level: [REQUIRED; DEFAULT="peptide"] Input level. Must be "peptide" or "PSM
#    filter_crap: [DEFAULT=TRUE] Filter out the features which match a cRAP protein
#    crap_fasta: [DEFAULT=NULL] Fasta file containing the cRAP proteins. Expects fasta header format thusly
#    >sp|cRAP002|P02768|ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2, e.g Uniprot in 3rd position ('|' delimited)
#    filter_associated_crap: [DEFAULT=TRUE] Filter out the features which match a cRAP associated protein
#
# Output: Filtered PD results
# ----------------------------------------------------------------------------------------------------------------------
parse_features <- function(infile,
                           sep="\t",
                           master_protein_col="Master.Protein.Accessions",
                           protein_col="Protein.Accessions",
                           unique_master=TRUE,
                           silac=FALSE,
                           TMT=FALSE,
                           level="peptide",
                           filter_crap=TRUE,
                           crap_fasta=NULL,
                           filter_associated_crap=TRUE){
  
  if(!level %in% c("PSM", "peptide")){
    stop("level must be PSM or peptide")
  }
  
  features_df <- read.delim(infile,  sep=sep, header=T, stringsAsFactors=FALSE)

  cat("Tally of features at each stage:\n")
  
  print_summaries(features_df, master_protein_col, "All features")
  
  features_df <- features_df %>% filter(UQ(as.name(master_protein_col))!="")
  print_summaries(features_df, master_protein_col, "Excluding features without a master protein")
  
  if(filter_crap){
    
    if(is.null(crap_fasta)){
      stop('must supply the crap fasta argument to filter cRAP proteins')
    }
    
    con <- file(crap_fasta, open="r")
    lines <- readLines(con)
    
    #print(lines %>% strsplit(split="\\|"))
    crap_proteins <- lines %>% strsplit("\\|") %>% lapply(function(x){
      if(substr(x[[1]],1,1)!=">"){
        return()
      }
      else{
        return(x[[3]])
      }
    }) %>% unlist()
    
    close(con)
    
    if(filter_associated_crap){
      associated_crap <- features_df %>%
        filter((UQ(as.name(master_protein_col)) %in% crap_proteins)|
                 grepl("cRAP", !!as.symbol(protein_col), ignore.case=FALSE)) %>%
        pull(UQ(as.name(protein_col))) %>%
        strsplit("; ") %>%
        unlist()
      associated_crap <- associated_crap[!grepl("cRAP", associated_crap)]
    
    }

    features_df <- features_df %>% filter(!UQ(as.name(master_protein_col)) %in% crap_proteins,
                                          !grepl("cRAP", !!as.symbol(protein_col), ignore.case=FALSE))
    print_summaries(features_df, master_protein_col, "Excluding features matching a cRAP protein")
    
    if(filter_associated_crap){
      cat(sprintf("Identified an additional %s proteins as 'cRAP associated'\n", length(associated_crap)))

      if(length(associated_crap)>0){
        # remove isoforms
        associated_crap_no_isoform <- unique(sapply(strsplit(associated_crap, "-"), "[[", 1))
        associated_crap_regex <- paste(associated_crap_no_isoform, collapse="|")
        features_df <- features_df %>% filter(!grepl(associated_crap_regex, UQ(as.name(protein_col))))
        print_summaries(features_df, master_protein_col, "Excluding features associated with a cRAP protein")
      }
    }
    
  }

  if(unique_master){
    features_df <- features_df %>% filter(Number.of.Protein.Groups==1)
    print_summaries(
      features_df, master_protein_col, "Excluding features without a unique master protein")
  }
  
  
  if(silac|TMT & level=="peptide"){
    features_df <- features_df %>% filter(Quan.Info=="")
    print_summaries(features_df, master_protein_col, "Excluding features without quantification")
  }
  return(features_df)
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
    abundance_columns <- abundance_columns[grep(sprintf('%s.*Count', quant_name), abundance_columns, invert=TRUE)]
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
  
  e_cols <- colnames(exprsCsv)
  p_rows <- rownames(pdataCsv)
  if(length(union(p_rows, e_cols))!=length(e_cols)){
    stop(sprintf(
      paste0('Samples_inf rownames should match abundance columns in data. ',
             'sample_inf rownames=%s,  Abundance columns=%s'),
      paste(p_rows, collapse=','),
      paste(e_cols, collapse=',')))
  }
  
  exprsCsv <- exprsCsv[,rownames(pdataCsv)]
  
  res <- MSnSet(as.matrix(sapply(exprsCsv, as.numeric)), fdataCsv, pdataCsv)
  
  summariseMissing(res)
  
  cat(sprintf("\n%s peptides do not have any quantification values\n\n", sum(rowSums(is.na(exprs(res)))==ncol(exprs(res)))))
  res <- res[rowSums(is.na(exprs(res)))!=ncol(exprs(res)),] # exclude peptides without any quantification
  return(res)
}

plotMissing <- function(obj, missing_col=cbPalette[3], present_col='black', verbose=TRUE, ...){
  tmp_obj <- obj
  exprs(tmp_obj)[!is.na(exprs(tmp_obj))] <- 1
  exprs(tmp_obj)[is.na(exprs(tmp_obj))] <- 0
  
  missing <- exprs(tmp_obj)

  missing <- missing[rowSums(missing==0)>0,,drop=FALSE] # identify features with missing values

  if(verbose){
    cat(sprintf("Out of %s total features, %s (%s%%) have missing values\n",
                length(rownames(exprs(tmp_obj))), length(rownames(missing)),
                round(100*length(rownames(missing))/length(rownames(exprs(tmp_obj))),3)))
    if(nrow(missing)>1){
      print(table(rowSums(missing==0)))
    }
  }

  if(nrow(missing)>1){
    missing_twice <- missing[rowSums(missing==0)>1,]
    
    if(verbose){
      cat(sprintf("And %s features have more than one missing value\n", length(rownames(missing_twice))))
    }
    
    colnames(missing) <- pData(tmp_obj)$Sample_name
    
    heatmap.2(missing, col = c(missing_col, present_col),
              scale = "none", dendrogram = "none",
              trace = "none", keysize = 0.5, key = FALSE,Colv = F, labRow=F,
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

  pep_agg <- combineFeatures(obj, groupBy=gb, fun=function(x) myAggFunction(x, fun))
  pData(pep_agg) <- pData(obj)
  return(pep_agg)
}


agg_to_peptides <- function(obj, gb=NULL, seq_col="Sequence", fun="median"){
  
  # remove any "CV.Abundance" columns. Otherwise, combineFeatures will throw an error
  fData(obj) <- fData(obj)[,grep("CV.*", colnames(fData(obj)), invert=TRUE)]
  
  if(missing(gb)){
    gb <- fData(obj)[[seq_col]]
  }
  
  pep_agg <- combineFeatures(obj, groupBy=gb, fun=fun, na.rm=TRUE)
  pData(pep_agg) <- pData(obj)
  return(pep_agg)
}


agg_to_protein <- function(obj, protein_col="Master.Protein.Accessions", gb=NULL, fun="median"){
  
  if(missing(gb)){
    gb <- fData(obj)[[protein_col]]
  }
  
  fData(obj) <- fData(obj)[,grep("CV.*", colnames(fData(obj)), invert=TRUE)]
  
  # we'll use our own aggregation function here to obtain the median over all the peptides
  pep_agg <- combineFeatures(obj, groupBy=gb, fun=fun, na.rm=TRUE)
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

# The functions below are used to replace sequential missing values with zeros
replace_missing_not_at_random_row <- function(x, min_sequential_missing){
  
  missing_rle <- rle(is.na(x))
  
  sequential_blocks <- missing_rle$values
  sequential_blocks[missing_rle$lengths<min_sequential_missing] <- FALSE
  
  replace_with_zero <- rep(sequential_blocks, missing_rle$lengths)
  
  out <- x
  out[replace_with_zero]<-0
  
  return(out)
}

replace_missing_not_at_random <- function(obj, max_total_missing=10, min_sequential_missing=4){
  
  missing_n <- rowSums(is.na(exprs(obj)))
  cat(sprintf("Removing %s features with more than %s missing values\n",
              sum(missing_n>max_total_missing), max_total_missing))
  obj <- obj[missing_n<=max_total_missing,]
  
  exprs(obj) <- exprs(obj) %>%
    apply(1, function(x) replace_missing_not_at_random_row(x, min_sequential_missing)) %>%
    t()
  
  zero_n <- rowSums(exprs(obj)==0, na.rm=TRUE)
  fData(obj)$zero_n <- zero_n
  cat("Number of sequential values imputed with zero\n")
  print(table(zero_n))
  # assuming that zeros are from the step above only...
  cat(sprintf("Sequential missing values imputed with zero in %s features\n", sum(zero_n>=min_sequential_missing)))
  
  
  return(obj)
}

filterPSM <- function(psm,
                      master_protein_col="Master.Protein.Accessions", # master protein column name
                      SN_threshold=5, # Signal:noise threshold to retain PSM quantification values
                      intensity_filter=2^2.25, # Minimum intensity value, e.g 2^2.25
                      interference_threshold=50, # Maximum interference
                      int_column="Isolation.Interference.in.Percent",
                      sn_column="Average.Reporter.SN"
                      
){
  ####################################################
  # Remove PSMs with Interference above threshold
  ####################################################
  if(!is.null(interference_threshold)){
    cat("Removing PSMs with high interference (e.g co-isolation)\n")

    psm_int <- psm[!is.na(fData(psm)[[int_column]]),]
    psm_int <- psm_int[fData(psm_int)[[int_column]]<=interference_threshold,]
    
    cat(sprintf("PSMs with co-isolation interference under threshold(%s): %s/%s\n",
                round(interference_threshold, 2), length(rownames(psm_int)), length(rownames(psm))))
    
    cat(sprintf("Proteins retained = %s/%s\n",
                length(unique(fData(psm_int)[[master_protein_col]])),
                length(unique(fData(psm)[[master_protein_col]]))))
  }
  else{
    psm_int <- psm
  }

  
  ####################################################
  # Remove PSMs with Signal:Noise (SN) below threshold
  ####################################################
  if(!is.null(SN_threshold)){
    cat("Removing PSMs with low Signal:noise ratio\n")
    psm_int_sn <- psm_int[fData(psm_int)[[sn_column]]>=SN_threshold,]
    cat(sprintf("PSMs with SN over threshold(%s): %s/%s\n",
                round(SN_threshold, 2), length(rownames(psm_int_sn)), length(rownames(psm_int))))
    
    cat(sprintf("Proteins retained = %s/%s\n",
                length(unique(fData(psm_int_sn)[[master_protein_col]])),
                length(unique(fData(psm_int)[[master_protein_col]]))))
  }
  else{
    psm_int_sn <- psm_int
  }
  ###################################################
  # Remove low intensity estimates
  ###################################################
  if(!is.null(intensity_filter)){
    cat("Replacing low intensity values with NA\n")
    psm_int_sn_filter_low <- psm_int_sn
    
    #plotLabelQuant(psm_int_sn_filter_low, log=TRUE, print=FALSE)$p2 +
    #  geom_vline(xintercept=intensity_filter) %>% print()
    
    p <- plotLabelQuant(psm_int_sn_filter_low, log=TRUE, print=FALSE)
    print(p$p1)
    print(p$p2 + geom_vline(xintercept=log2(intensity_filter)))
    
    exprs(psm_int_sn_filter_low)[exprs(psm_int_sn_filter_low)<intensity_filter] <- NA
  }
  else{
    psm_int_sn_filter_low <- psm_int_sn
  }
  
  p <- plotLabelQuant(psm_int_sn_filter_low, log=TRUE, print=FALSE)$p2 +
    geom_vline(xintercept=log2(intensity_filter))
  
  print(p)
  
  invisible(psm_int_sn_filter_low)
}



# Given a MSnSet with peptide level abundance, remove peptides which
# are assigned to proteins with too few peptides assigned in total
restrictPepPerProtein <- function(peptide_obj, min_peptides, plot=TRUE){
  pep2protein <- fData(peptide_obj) %>%
    select(Master.Protein.Accessions)
  
  n_pep_per_prot <- tidy(peptide_obj, addPheno=TRUE) %>%
    merge(pep2protein, by.x="protein", by.y="row.names") %>%
    filter(is.finite(value)) %>%
    group_by(Master.Protein.Accessions, sample) %>%
    tally()
  
  if(plot){
    p <- ggplot(n_pep_per_prot, aes(log2(n))) +
      geom_histogram() +
      my_theme
    print(p)
  }
  
  retain_mask <- n_pep_per_prot %>%
    merge(fData(peptide_obj)[,c('Sequence', 'Master.Protein.Accessions'), drop=FALSE],
          by='Master.Protein.Accessions') %>%
    mutate(retain=n>=min_peptides) %>%
    select(Sequence, sample, retain) %>%
    spread(key=sample, value=retain) %>%
    tibble::column_to_rownames('Sequence')
  
  retain_mask[is.na(retain_mask)] <- FALSE
  retain_mask <- retain_mask[rownames(peptide_obj),colnames(peptide_obj)]
  
  masked_exprs <- exprs(peptide_obj) * retain_mask
  masked_exprs[masked_exprs==0] <- NA
  
  exprs(peptide_obj) <- as.matrix(masked_exprs)
  
  retain_proteins <- n_pep_per_prot %>%
    filter(n>=min_peptides) %>%
    pull(Master.Protein.Accessions)
  
  peptide_obj <- peptide_obj[fData(peptide_obj)$Master.Protein.Accessions %in% retain_proteins,]
  invisible(peptide_obj)
}


