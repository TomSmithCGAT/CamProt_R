# Utility functions useful across all analyses
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(VennDiagram))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

my_theme <- theme_bw() + theme(
  text=element_text(size=20,  family="serif"),
  panel.grid=element_blank(), aspect.ratio=1,
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


### parsePTMScores ###
# Function to parse the PD PTM score probabilities column
# and add columns for downstream filtering
# obj = data.frame with PD output
# threshold = score threshold
# ptm_col = colname for PTM scores
# prob_split = regex to split PTM probabilities
# collapse_delimiter = delimiter between multiple values in output columns
# verbose = print log of events

# to do:
# add logging?

parsePTMScores <- function(obj, threshold=95, ptm_col="PhosphoRS..Best.Site.Probabilities",
                           prob_split='; |: ', collapse_delimiter=";", verbose=TRUE){
  
  if(class(obj)!="data.frame"){
    stop("first argument must be a data.frame")
  }
  
  # initiate vectors with empty string
  # where the PTM(s) pass the threshold, these will be updated
  filtered_ptm_desc <- filtered_ptm_res <- filtered_ptm_pos <- filtered_ptm <- filtered_ptm_score <- rep("", nrow(df))
  
  split_probabities <- strsplit(obj[[ptm_col]], split=prob_split)
  
  log <- list("Total peptides"=0, "Retained peptides"=0, "Filtered peptides"=0,
              "..but some sites above threshold"=0,
              "Total sites"=0, "Retained sites"=0, "Filtered sites"=0,
              "Too many isoforms"=0)

  for (i in seq_along(split_probabities)){
    
    peptide_ptm_scores <- split_probabities[[i]]
    
    log$`Total peptides` <- log$`Total peptides` + 1
    
    if (peptide_ptm_scores[[1]] == "Too many isoforms"){
      log$`Too many isoforms` <- log$`Too many isoforms` + 1
      log$`Filtered peptides` <- log$`Filtered peptides` + 1
      next()
    }

    log$`Total sites` <- log$`Total sites` + length(peptide_ptm_scores)/2
    
    scores <- peptide_ptm_scores[seq(2, length(peptide_ptm_scores), 2)]
    
    # if any score is below threshold, disregard all putative ptm sites
    if (any(as.numeric(scores)<threshold)){
      log$`Filtered sites` <- log$`Filtered sites` + length(peptide_ptm_scores)/2
      log$`Filtered peptides` <- log$`Filtered peptides` + 1
      if (any(as.numeric(scores)>=threshold)){
        log$`..but some sites above threshold` <- log$`..but some sites above threshold` + 1  
      }
      # if we want to handle this differently, can implement an alternative approach here
      # and move the rest of the code below into an else clause
      next()
    }
    
    log$`Retained sites` <- log$`Retained sites` + length(peptide_ptm_scores)/2
    log$`Retained peptides` <- log$`Retained peptides` + 1
    
    ptms <- peptide_ptm_scores[seq(1, length(peptide_ptm_scores), 2)] # extract the PTMs info
    split_ptms <- unlist(strsplit(ptms, split =  '\\(|\\)')) # split to remove parantheses
    modifications <- split_ptms[seq(2, length(split_ptms), 2)] # extract modifications, e.g "phospho"
    positions <- split_ptms[seq(1, length(split_ptms), 2)] # extract the positions, e.g "S6"
    residues <- substr(positions, 1, 1) # extract first element, e.g S
    positions <- sub('.', '', positions) # remove first element and leave position, e.g 6
    
    # paste together the value, separated by option(collapse_delimiter) and update vectors which will become columns
    filtered_ptm_desc[[i]] <- paste(peptide_ptm_scores, collapse=collapse_delimiter)
    filtered_ptm_res[[i]] <- paste(residues, collapse=collapse_delimiter)
    filtered_ptm_pos[[i]] <- paste(positions, collapse=collapse_delimiter)
    filtered_ptm[[i]] <- paste(modifications, collapse=collapse_delimiter)
    filtered_ptm_score[[i]] <- paste(scores, collapse=collapse_delimiter)
    
  }
  
  # add columns
  obj['filtered_PTM_desc'] = filtered_ptm_desc
  obj['filtered_res'] = filtered_ptm_res
  obj['filtered_pos'] = filtered_ptm_pos
  obj['filtered_ptm'] = filtered_ptm
  obj['filtered_score'] = filtered_ptm_score
  
  if(verbose){
    for(event in names(log)){
      print(sprintf("%s: %i", event, log[[event]]))
    }
  }
  
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

