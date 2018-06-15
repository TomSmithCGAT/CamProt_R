# Utility functions useful across all analyses
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(VennDiagram))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

my_theme <- theme_bw() + theme(
  text=element_text(size=20,  family="serif"), panel.grid=element_blank(), aspect.ratio=1,
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
    #stat_summary(geom="point", fun.y=mean) +
    scale_colour_discrete(guide=F) +
    xlab("")+ ylab(ylb) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
    facet_wrap(~lab)
  
  p2 <- p + geom_violin() + geom_boxplot(notch=TRUE, outlier.size=NA, width=0.1)
  suppressMessages(print(p2))
  suppressMessages(print(p2 %+% final_df + aes(x=Var2)))
  
  p3 <- p +
    aes(fill=lab) +
    #geom_quasirandom(colour="black", pch=21, size=2, stroke = 0.2) +
    scale_fill_discrete(guide=F) +
    #stat_summary(geom="point", fun.y = median, colour="black", size=3) +
    stat_summary(geom="pointrange", colour="black", size=0.5)#width=0.3, lwd=1, size=1)
  #fun.ymin = function(z) {quantile(z,0.25)},
  #fun.ymax = function(z) {quantile(z,0.75)})
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

