# emulate standard pROLOC colours

suppressMessages(library("pRoloc"))
suppressMessages(library("reshape2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(require(gridExtra))

setStockcol(NULL) ## reset first
setStockcol(paste0(getStockcol(), 70))

getClassColours <- function(){
  cols = getStockcol()[c(1:5,12,7:11,13:20)]
  return(cols)
  
}

organelle_order <- c("CYTOSOL", "Cytosol", "PROTEASOME", "PROTEASOME 19S", "PROTEASOME 20S", "RIBOSOME", "RIBOSOME 40S", "RIBOSOME 60S",
                     "er", "ER", "GOLGI", "GA", "LYSOSOME", "PM", "PEROXISOME", "MITOCHONDRION", "MITOCHONDRIA", "Mitochondria",
                     "Mitochondrion", "NUCLEUS/CHROMATIN", "NUCLEUS", "Nuclear", "NUCLEUS-CHROMATIN", "CHROMATIN",
                     "unknown", "missing")

organelle2colour <- list("CYTOSOL"="#E41A1C",
                         "Cytosol"="#E41A1C",
                         "PROTEASOME"="#984EA3",
                         "PROTEASOME 19S"="#704da2",
                         "PROTEASOME 20S"="#a24d70",
                         "RIBOSOME"="#9999FF",
                         "RIBOSOME 40S"="#9999FF",
                         "RIBOSOME 60S"="#000099",
                         "ER"="#FF7F00",
                         "er"="#FF7F00",
                         "GOLGI"="#377EB8",
                         "GA"="#377EB8",
                         "LYSOSOME"="#F781BF",
                         "PM"="#00CED1",
                         "PEROXISOME"="#A65628",
                         "MITOCHONDRION"="#FFD700",
                         "MITOCHONDRIA"="#FFD700",
                         "Mitochondria"="#FFD700",
                         "Mitochondrion"="#FFD700",
                         "NUCLEUS/CHROMATIN"="#9ACD32",
                         "NUCLEUS-CHROMATIN"="#238B45",
                         "NUCLEUS"="#9ACD32",
                         "Nuclear"="#9ACD32",
                         "CHROMATIN"="#238B45",
                         "unknown"="grey50",
                         "missing"="grey50")

# colour-blind friendly pallete when to be used when not plotting organelles
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

getColours <-function(organelles){
  unlist(organelle2colour[names(organelle2colour) %in% organelles], use.names = FALSE)
}

my_theme <- theme_bw() + theme(text=element_text(size=20), aspect.ratio=1, panel.grid=element_blank())

# ------------------------------------------------------------------------------------------------------------
# Function	: PlotProteinProfiles 
# ------------------------------------------------------------------------------------------------------------

PlotProteinProfiles <- function(res){
  
  exprs_df <- exprs(res)
  exprs_df <- melt(exprs_df)
  
  p <- ggplot(exprs_df, aes(Var2, value, colour=Var1, group=Var1)) +
    my_theme + geom_line() +
    theme(aspect.ratio=1, axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  +
    xlab("") + ylab("norm. abundance") + 
    scale_colour_discrete(name="", na.value="grey")
  
  return(p)
}

# ------------------------------------------------------------------------------------------------------------
# Function	: PlotMarkerProfiles 
# ------------------------------------------------------------------------------------------------------------

PlotMarkerProfiles <- function(res, fcol="master_protein",
                               keep_markers=c("CYTOSOL"),
                               unknown=F, plot_all=T,
                               individual_plots=F, foi=NULL,
                               organelle_order=NULL, alpha=0.25,
                               foi_alpha=1){
  
  
  exprs_df <- exprs(res)
  f_df <- fData(res)
  f_df$markers <- f_df[[fcol]]
  exprs_df <- melt(exprs_df)
  exprs_df <- merge(exprs_df, f_df, by.x="Var1", by.y="row.names")
  
  if(unknown){
    exprs_trim_df <- exprs_df[exprs_df$markers %in% c(keep_markers, "unknown"),]
    exprs_trim_df$unknown <- exprs_trim_df$markers=="unknown"
    exprs_trim_df$markers[exprs_trim_df$markers=="unknown"] <- NA
  }
  else{
    exprs_trim_df <- exprs_df[exprs_df$markers %in% keep_markers,]
  }
  
  if(!missing(organelle_order)){
    exprs_trim_df$markers <- factor(exprs_trim_df$markers, levels=organelle_order)
  }
  
  p <- ggplot(exprs_trim_df, aes(factor(Var2), value, colour=markers)) +
    my_theme +
    theme(aspect.ratio=1, axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  +
    xlab("") + ylab("norm. abundance") + 
    scale_colour_discrete(name="", na.value="grey")
  
  if(plot_all){
    p <- p + aes(group=Var1)
    
    if(unknown){
      p <- p + geom_line(aes(alpha=unknown)) +
        scale_alpha_manual(values=c(1, 0.01), guide=F)
    }
    else{
      p <- p + geom_line(alpha=alpha)
    }
    if(individual_plots){
      p <- p + facet_wrap(~Var1) + scale_colour_discrete(guide=F) +
        geom_line(alpha=1) +
        theme(text=element_text(size=10), panel.grid=element_blank(),
              axis.text=element_blank())
    }
  }
  else{
    if(unknown){
      p <- p + stat_summary(aes(group=markers), geom="line", fun.y=mean, size=1) +
        scale_alpha_manual(values=c(1, 0.01), guide=F)
    }
    else{
      p <- p + stat_summary(aes(group=markers), geom="line", fun.y=mean, size=1)
    }
  }
  
  if(!missing(foi)){
    exprs_df <- data.frame(exprs_df)
    features_df <- exprs_df[exprs_df$Var1 %in% foi,]
    #if(plot_all){
    p <- p + geom_line(data=features_df,
                       aes(factor(Var2), value, group=Var1),
                       colour="black", alpha=foi_alpha)
    #}
    #else{
    #p <- p + stat_summary(data=features_df, aes(Var2, value), group="1",
    #                      geom="line", fun.y=mean, colour="black", alpha=1)
    #}
    
  }
  
  p <- p + guides(colour=guide_legend(override.aes = list(alpha = 1, size=1)))
  
  return(p)
}

PlotMarkerHeatmap <- function(res, m_col="markers", organelle_order=NULL){
  
  exprs_df <- exprs(res)
  f_df <- fData(res)
  f_df$markers <- f_df[[m_col]]
  exprs_df <- exprs_df %>% fitler(markers!="unknown")
  long_exprs_df <- melt(exprs_df)
  colnames(long_exprs_df) <- c("feature", "fraction", "intensity")
  long_exprs_df <- merge(long_exprs_df, f_df, by.x="feature", by.y="row.names")
  long_exprs_df <- long_exprs_df %>% group_by(markers, fraction) %>% summarise(intensity=mean(intensity))
  
  if(missing(organelle_order)){
    
    average_profile_matrix <- exprs_df %>% spread(key=Var2, value=intensity) %>% data.frame()
    rownames(average_profile_matrix) <- average_profile_matrix$markers
    average_profile_matrix <- average_profile_matrix %>% dplyr::select(-markers)
    
    clust <- hclust(dist(as.matrix(average_profile_matrix)))
    plot(clust)
    my_order <- rownames(average_profile_matrix)[clust$order]
    
    exprs_df$markers <- factor(exprs_df$markers, my_order)
  }
  
  else{
    exprs_df$markers <- factor(exprs_df$markers, organelle_order)
  }
  
  p <- ggplot(exprs_df, aes(Var2, markers, fill=intensity)) +
    geom_tile() +
    scale_fill_continuous(low="white", high=cbPalette[3], limits=c(0,NA), name="Normalised\nabundance") +
    my_theme +
    xlab("Fraction") +
    ylab("") +
    theme(axis.text.x=element_text(size=12, angle=45, vjust=1, hjust=1), axis.text.y=element_text(size=15))
  
  return(p)
  
}
# ------------------------------------------------------------------------------------------------------------
# Function	: makeQSepDistance
# ------------------------------------------------------------------------------------------------------------

makeQSepDistance <- function(res_with_markers, plot_cluster=F, ...){
  qsep_dist <- qsep(QSep(res_with_markers, ...))
  clustering <- hclust(as.dist(qsep_dist))
  
  if(plot_cluster){
    plot(clustering, hang = -1)
  }
  
  qsep_df <- melt(qsep_dist)
  qsep_df$Var1 <- factor(qsep_df$Var1, levels = clustering$labels[clustering$order])
  qsep_df$Var2 <- factor(qsep_df$Var2, levels = clustering$labels[clustering$order])
  
  return(qsep_df)
}

# ------------------------------------------------------------------------------------------------------------
# Function	: makeQSepDistancePlot 
# ------------------------------------------------------------------------------------------------------------
makeQSepDistancePlot <- function(qsep_df){
  
  p <- ggplot(qsep_df, aes(Var1, Var2, fill=value)) +
    geom_tile() + my_theme +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), aspect.ratio = 1) +
    theme(axis.text=element_text(size=10)) + xlab("") + ylab("") +
    scale_fill_continuous(low="white", high="steelblue", name="QSep\ndistance")
  
  return(p)
  
}

# ------------------------------------------------------------------------------------------------------------
# Function	: makeQSepDistancePlot 
# ------------------------------------------------------------------------------------------------------------

plotCompareDistances <- function(qsep_df1, qsep_df2, exp1="All", exp2="TMT2",
                                 qsep_df3=F, exp3=F, qsep_df4=F, exp4=F,
                                 membrane_organelles=FALSE, facet_out=F, plot_den=F, organelle_order=NULL, remove=F){
  
  if(membrane_organelles==FALSE){
    membrane_organelles <- c("ER", "PEROXISOME", "PM", "MITOCHONDRIA", "GOLGI", "LYSOSOME")
  } 
  
  qsep_df1$samples <- exp1
  qsep_df2$samples <- exp2
  qsep_cat_df <- rbind(qsep_df1, qsep_df2)
  
  if(!qsep_df3==FALSE){
    qsep_df3$samples <- exp3
    qsep_cat_df <- rbind(qsep_cat_df, qsep_df3)
  }
  
  if(!qsep_df4==FALSE){
    qsep_df4$samples <- exp4
    qsep_cat_df <- rbind(qsep_cat_df, qsep_df4)
  }
  
  if(!missing(organelle_order)){
    qsep_cat_df$Var1 <- factor(qsep_cat_df$Var1, levels=organelle_order)
    qsep_cat_df$Var2 <- factor(qsep_cat_df$Var2, levels=organelle_order)
  }
  
  if(!missing(organelle_order)){
    qsep_cat_df <- qsep_cat_df %>% filter(Var1!=remove, Var2!=remove)
  }
  
  plots <-function(qsep_cat_df, colour_shape_same=T, x, y, facet_out=F){
    
    p <- makeQSepDistancePlot(qsep_cat_df)
    p_facet <- p + facet_wrap(~samples) +
      theme(axis.text=element_text(size=8), text=element_text(size=15))
    print(p_facet)
    if(facet_out!=FALSE){
      ggsave(p_facet, filename=facet_out)
    }
    p <- ggplot(qsep_cat_df, aes(value, colour=samples)) + geom_density() + my_theme +
      xlab("Distance")
    
    if(plot_den){
      print(p)
    }
  
    qsep_cat_wide_df <- dcast(Var1+Var2~samples, value.var="value", data=qsep_cat_df)
    qsep_cat_wide_df$label <- qsep_cat_wide_df$Var2
    qsep_cat_wide_df$label[abs(qsep_cat_wide_df[[exp2]] - qsep_cat_wide_df[[exp1]])<2] <- ""
    
    qsep_cat_wide_df <- qsep_cat_wide_df[qsep_cat_wide_df$Var1 != qsep_cat_wide_df$Var2,]
    
    p <- ggplot(qsep_cat_wide_df, aes_string(x, y)) +
      geom_point(size=3) + my_theme + geom_abline(slope=1, colour="red", linetype=2)
    
    if(colour_shape_same){
      p <- p + aes(group=interaction(Var1, Var2), colour=Var1, shape=Var1) +
        scale_shape_manual(values=c(1:20), name="") +
        scale_colour_discrete(name="")
      
    }
    else{
      p <- p + aes(group=interaction(Var1, Var2), colour=Var1, shape=Var2) +
        scale_shape_discrete(name="marker set 1") +
        scale_colour_discrete(name="marker set 2")
      
    }
    
    print(p)
    
  }
  
  plots(qsep_cat_df, x=exp1, y=exp2, facet_out=facet_out)
  plots(qsep_cat_df[(qsep_cat_df$Var1 %in% membrane_organelles & qsep_cat_df$Var2 %in% membrane_organelles),],
        colour_shape_same=FALSE, x=exp1, y=exp2)
  
  if(!qsep_df3==FALSE){
    plots(qsep_cat_df, x=exp2, y=exp3)
    plots(qsep_cat_df[(qsep_cat_df$Var1 %in% membrane_organelles & qsep_cat_df$Var2 %in% membrane_organelles),],
          colour_shape_same=FALSE, x=exp2, y=exp3)
    
    plots(qsep_cat_df, x=exp1, y=exp3)
    plots(qsep_cat_df[(qsep_cat_df$Var1 %in% membrane_organelles & qsep_cat_df$Var2 %in% membrane_organelles),],
          colour_shape_same=FALSE, x=exp1, y=exp3)
    
  }
  
}


# ------------------------------------------------------------------------------------------------------------
# Function	: plotMinDistance
# ------------------------------------------------------------------------------------------------------------
plotMinDistance <- function(qsep1, qsep2, order_organelle=NULL){

  qsep1$type="qsep1"
  qsep2$type="qsep2"
    
  qsep_min_comparison <- rbind(qsep1, qsep2) %>% filter(Var1!=Var2) %>% group_by(type, Var1) %>%
  summarise(min_distance=min(value)) %>%
  spread(key=type, value=min_distance) %>%
  mutate(Var1=factor(Var1, levels=sort(as.character(Var1)))) %>%
  data.frame()

  if(!missing(order_organelle)){
    qsep_min_comparison$Var1 <- factor(qsep_min_comparison$Var1, levels=order_organelle)
  }
  p <- qsep_min_comparison %>%
    ggplot(aes(qsep1, qsep2)) +
    geom_abline(slope=1, linetype=2, colour="grey50") +
    geom_point(aes(colour=Var1), size=6) +
    my_theme
  
  return(p)
}

# ------------------------------------------------------------------------------------------------------------
# Function	: plotPCA
# ------------------------------------------------------------------------------------------------------------
makePCA_proj <- function(...){
  make_proj("PCA", ...)
}

make_proj <- function(method="PCA", res, fcol="markers", dims=c(1,2), ...){
  PCA_matrix <- plot2D(res, method=method, fcol=fcol, plot=FALSE, dims=dims, ...)
  axes <- colnames(PCA_matrix)
  colnames(PCA_matrix) <- c("X", "Y")
  PCA_df <- PCA_matrix %>% merge(data.frame(fData(res)), by="row.names")
  PCA_df$markers <- PCA_df[[fcol]]
  PCA_df$unknown <- PCA_df$markers=="unknown"
  #PCA_df$markers[PCA_df$unknown] <- NA
  PCA_df <- PCA_df[order(-PCA_df$unknown),]
  return(list("PCA_df"=PCA_df, "axes"=axes))
}

LOPITPCAPlotter <- function(...){
  LOPITPlotter(...)
}


LOPITPlotter <- function(
  PCA_df, axes, xlims=FALSE, ylims=FALSE, foi=FALSE, add_foi_names=FALSE,
  just_markers=FALSE, m_colours=FALSE, re_order_markers=FALSE, marker_levels=NULL,
  title=FALSE, fcol=FALSE, point_size=FALSE, foi_size=1, foi_alpha=1, foi_shape=8,
  foi_colour="black", foi_fill="black", return_proj=FALSE, cex=1, add_density=FALSE,
  unknown_name="unknown", unknown_colour="grey70", show_unknown=TRUE){
  
  if(!missing(fcol)){
    # remake marker column
    PCA_df$markers <- PCA_df[[fcol]]
    PCA_df$unknown <- PCA_df$markers=="unknown"
    PCA_df$markers[PCA_df$unknown] <- unknown_name
    PCA_df <- PCA_df[order(-PCA_df$unknown),]
  }
  else{
    PCA_df$unknown <- PCA_df$markers=="unknown"
    PCA_df$markers[PCA_df$unknown] <- unknown_name
  }
  
  if(re_order_markers!=FALSE){
    marker_levels <- setdiff(marker_levels, unknown_name)
    PCA_df$markers <- factor(
      PCA_df$markers, levels=c(marker_levels, unknown_name))
  }
  else{
    PCA_df$markers <- factor(
      PCA_df$markers,
      levels=c(setdiff(unique(PCA_df$markers), unknown_name), unknown_name))
  }
  
  plots = NULL
  if(m_colours==FALSE){
    m_colours <- getStockcol()[c(1:5,12,7:11,13:20)][
      1:length(setdiff(unique(PCA_df$markers), unknown_name))]
  }
  p <- ggplot(PCA_df) +
    geom_point(aes(X, Y, colour=markers, alpha=unknown), size=cex) +
    scale_alpha_manual(values=c(0.5,0.1), guide=F) +

    my_theme + xlab(axes[1]) + ylab(axes[2])

  if(!missing(point_size)){
    p <- p + aes_string(size=point_size) +
      scale_size_area(guide=F, max_size=2)
    print(p)
  }
  
  if(show_unknown){
    p <- p + scale_colour_manual(
      values=c(m_colours,  unknown_colour), name="")}
  else{
    p <- p + scale_colour_manual(
      values=c(m_colours,  unknown_colour), name="",
      breaks = setdiff(levels(PCA_df$markers), unknown_name))}
  
  if(!missing(title)){
    p <- p + ggtitle(title)
  }
  
  if(xlims){
    p <- p + xlim(xlims)
  }
  
  if(ylims){
    p <- p + ylim(ylims)
  }
  
  if(add_density){
    p <- p + geom_density_2d(aes(x=X, y=Y), colour="black", alpha=0.5)
  }
  
  plots[['p']] <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size=2)))
  
  if(just_markers){
    PCA_df2 <-PCA_df[!is.na(PCA_df$markers),]
    p2  <- p %+% PCA_df2
    plots[['p2']] <- p2 + guides(colour = guide_legend(override.aes = list(alpha = 1, size=2)))
    
  }
  
  if(!missing(foi)){
    features_df <- PCA_df[PCA_df[['Row.names']] %in% foi,]
    p <- p + geom_point(data=features_df, aes(X, Y), shape=foi_shape, size=foi_size, alpha=foi_alpha,
                        colour=foi_colour, fill=foi_fill, stroke=0.4)
    
    if(!missing(point_size)){
      p <- p + aes_string(size=point_size)
    }
    
    if(add_foi_names){
      p <- p + geom_text(data=features_df, aes(X, Y+0.25, label=Row.names))
    }
    
    plots[['p_foi']] <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size=2)))
    
    if(just_markers){
      features_df2 <- PCA_df2[PCA_df2[['Row.names']] %in% foi,]
      p2 <- p2 + geom_point(data=features_df2, aes(X,Y), shape=8)
      
      if(!missing(point_size)){
        p2 <- p2 + aes_string(size=point_size)
      }
      
      plots[['p2_foi']] <- p2 + guides(colour = guide_legend(override.aes = list(alpha = 1, size=2)))
    }
  }
  
  return(plots)
}

plotPCA <- function(res, fcol, dims=c(1,2), ...){
  
    PCA <- makePCA_proj(res, fcol, dims)
    LOPITPCAPlotter(PCA$PCA_df, PCA$axes, ...)
}

plotHexbin <- function(res, fcol, dims=c(1,2), ...){
  
  PCA <- make_proj("PCA", res, fcol, dims)
  LOPITPCAPlotter_hexbin(PCA$PCA_df, PCA$axes, ...)
}

LOPITPCAPlotter_hexbin <- function(PCA_df, axes){
  p <- ggplot(PCA_df) +
    geom_hex(aes(X, Y)) +
    scale_alpha_manual(values=c(0.5,0.1), guide=F) +
    scale_fill_continuous(low="grey90", high=cbPalette[6]) +
    my_theme + xlab(axes[1]) + ylab(axes[2])
  
  invisible(p)  
}

# ------------------------------------------------------------------------------------------------------------
# Function	: checkParams
# ------------------------------------------------------------------------------------------------------------

# Use times = 100 for final results but takes a long time (~40 mins for each function call!)
checkParams <- function(res_with_markers, ml_method=svmOptimisation,
                        desc="", ml_desc="svm", fcol="master_protein", times=100, weight=F, plot_output=F,
                        ...){
  
  exprs(res_with_markers)[is.na(exprs(res_with_markers))] <- 0
  
  sink("tmp") # we don't want the output
  if(weight==F){
    params <- ml_method(res_with_markers, fcol=fcol, times=times,
                        verbose=T, ...)
  }
  else{
    params <- ml_method(res_with_markers, fcol=fcol, times=times,
                        verbose=T, class.weights=classWeights(res_with_markers, fcol=fcol), ...)
  }
  sink()
  
  if(plot_output){
    print(plot(params))
    print(levelPlot(params))
  }
  
  all_F1_matrix <- do.call("rbind", params@f1Matrices)
  all_F1_df <- data.frame(melt(all_F1_matrix))
  colnames(all_F1_df) <- c("sigma", "cost", "F1")
  F1_df <-data.frame(getF1Scores(params))
  F1_df$desc <- desc
  F1_df$weight <- weight
  F1_df$ml <- ml_desc
  
  # code below copied from Ollie
  f1scoreperO <- matrix(0, times, dim(params@cmMatrices[[1]])[1])
  
  for(i in 1:times){
    
    conf <- params@cmMatrices[[i]]
    
    f1perO <- MLInterfaces::F1(conf, naAs0 = TRUE)
    
    f1scoreperO[i, ] <- f1perO
  }
  
  OrganelleF1 <- as.data.frame(f1scoreperO)
  colnames(OrganelleF1) <- getMarkerClasses(res_with_markers, fcol=fcol)
  OrganelleF1$desc <- desc
  OrganelleF1$weight <- weight
  
  return(list("F1_df"=F1_df, "OrganelleF1"=OrganelleF1))
  
}

# input is concatenate F1_dfs from multiple calls to checkParams
makeF1plots <- function(F1_df, facet_by_ml=TRUE){
  
  plots <- NULL
  
  p1 <- ggplot(F1_df, aes(desc, F1, colour=desc)) +
    scale_colour_discrete(guide=FALSE) +
    my_theme +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), aspect.ratio = 1) +
    xlab("")
  
  if(facet_by_ml){
    p1 <- p1 + facet_wrap(~ml)
  }
  
  p1a <- p1 + stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.data = mean_se, geom = "errorbar")
  plots[['p1a']] <- p1a
  plots[['p1b']] <- p1 + geom_boxplot(notch=TRUE)
  

  p2 <- ggplot(F1_df, aes(desc, F1, colour=ml)) +
    my_theme +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), aspect.ratio = 1) +
    xlab("")
  
  p2a <- p2 + stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.data = mean_se, geom = "errorbar")
  
  plots[['p2a']] <- p2a
  plots[['p2b']] <- p2b <- p2 + geom_boxplot(notch=TRUE)
  
  return(plots)
}

# ------------------------------------------------------------------------------------------------------------
# Function	: plotBothLOPITs
# ------------------------------------------------------------------------------------------------------------

plotBothLOPITs <- function(LOPIT1, LOPIT2, foi=F, xlims=c(-5.5, 8), ylims=c(-6,6)){
  num_iden <- length(intersect(foi, row.names(exprs(LOPIT1))))
  num_iden2 <- length(intersect(foi, row.names(exprs(LOPIT2))))
  cat(sprintf(
    "Out of %s proteins of interest, %s are observed in LOPIT1 and %s are observed in LOPIT2\n",
    length(foi), num_iden, num_iden2))
  
  plotPCA(LOPIT1, xlims=xlims, ylims=ylims, foi=foi)
  plotPCA(LOPIT2, fcol="Accession", foi=foi)
}

# ------------------------------------------------------------------------------------------------------------
# Function	: plot collection of PCs
# ------------------------------------------------------------------------------------------------------------

plotPC1to_4 <- function(res, just_markers=F, ...){
  p1 <- plotPCA(res, just_markers=just_markers, ...)
  p2 <- plotPCA(res, dims=c(3,4), just_markers=just_markers, ...)
  if(just_markers){
    grid.arrange(p1$p2, p2$p2, ncol=1)}
  else{
    grid.arrange(p1$p, p2$p, ncol=1)
  }
}

plotPC1to_8 <- function(res, just_markers=F, ...){
  p1 <- plotPCA(res, just_markers=just_markers, ...)
  p2 <- plotPCA(res, dims=c(3,4), just_markers=just_markers, ...)
  p3 <- plotPCA(res, dims=c(5,6), just_markers=just_markers, ...)
  p4 <- plotPCA(res, dims=c(7,8), just_markers=just_markers, ...)
  if(just_markers){
    grid.arrange(p1$p2, p2$p2, p3$p2, p4$p2, ncol=2)}
  else{
    grid.arrange(p1$p, p2$p, p3$p, p4$p, ncol=2)
  }
}

# ----------------------------------------
# Function	: plot per-organelle F1 scores
# ----------------------------------------

plotOrganelleF1 <- function(melted_F1_organelle){
  p <- ggplot(melted_F1_organelle, aes(variable, value, colour=desc)) + 
    geom_boxplot() +
    my_theme +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +
    scale_x_discrete(name="") +
    xlab("") + ylab("F1") +
    scale_colour_discrete(name="")
  
  print(p)
}


catLOPITs <- function(lopit_a, lopit_b, fdata_merge="master_protein"){
  intersect_proteins <- intersect(rownames(exprs(lopit_a)), rownames(exprs(lopit_b)))
  
  subsetToCommonExprs <- function(lopit, common_proteins){
    exprs_data <- exprs(lopit)
    exprs_data <- exprs_data[rownames(exprs_data) %in% common_proteins,]
    exprs_data <- exprs_data[common_proteins, ]
    return(exprs_data)
  }
  
  exprs_df1 <- subsetToCommonExprs(lopit_a, intersect_proteins)
  colnames(exprs_df1) <- paste0("1_", colnames(exprs_df1))
  exprs_df2 <- subsetToCommonExprs(lopit_b, intersect_proteins)
  colnames(exprs_df2) <- paste0("2_", colnames(exprs_df2))
  
  print(sum(rownames(exprs_df1) != rownames(exprs_df2))) # manual check for mismatch in protein index
  
  cat_exprs_df <- cbind(exprs_df1, exprs_df2)
  
  pdata_df1 <- data.frame(pData(lopit_a))
  rownames(pdata_df1) <- paste0("1_", rownames(pdata_df1))
  pdata_df2 <- data.frame(pData(lopit_b))
  rownames(pdata_df2) <- paste0("2_", rownames(pdata_df2))
  
  cat_pData <- rbind(pdata_df1, pdata_df2)
  
  subsetToCommonfData <- function(lopit_norm, common_proteins){
    fdata_df <- data.frame(fData(lopit_norm))
    fdata_df <- fdata_df[rownames(fdata_df) %in% common_proteins,]
    fdata_df <- fdata_df[common_proteins, ]
    return(fdata_df)
  }
  
  fdata_df1 <- subsetToCommonfData(lopit_a, intersect_proteins)
  fdata_df2 <- subsetToCommonfData(lopit_b, intersect_proteins)
  cat_fdata <- merge(fdata_df1, fdata_df2, by=fdata_merge)
  rownames(cat_fdata) <- cat_fdata[[fdata_merge]]
  
  cat_LOPIT <- MSnSet(cat_exprs_df, cat_fdata, cat_pData)

  return(cat_LOPIT)
  
}


getPCALoadings <- function(obj, ...){
  .pca <- prcomp(exprs(obj), ...)
  loadings <- .pca$rotation
  return(loadings)
}  

plotLoadings <- function(loadings){
  
  p <- melt(loadings) %>%
    ggplot(aes(factor(Var1, levels=unique(Var1)), value)) +
    geom_bar(stat="identity") +
    facet_wrap(~Var2) +
    my_theme +
    theme(axis.text.x=element_text(size=10, angle=45, vjust=1, hjust=1))
  
  return(p)
  
}

##' The function plots marker consensus profiles obtained from mrkConsProfile
##'
##' @title Plot marker consenses profiles.
##' @param object A \code{matrix} containing marker consensus profiles as output
##'   from \code{\link{mrkConsProfiles}}.
##' @param order Order for markers (Optional).
##' @param plot A \code{logical} defining whether the heatmap should be plotted.
##'   Default is \code{TRUE}.
##' @return Invisibly returns \code{\link[ggplot2]{ggplot}} object.
##' @author Tom Smith
##' @examples
##' library("pRolocdata")
##' data(E14TG2aS1)
##' hc <- mrkHClust(E14TG2aS1, plot=FALSE)
##' mm <- getMarkerClasses(E14TG2aS1)
##' m_order <- levels(factor(mm))[order.dendrogram(hc)]
##' fmat <- mrkConsProfiles(E14TG2aS1)
##' plotConsProfiles(fmat, order=m_order)
plotConsProfiles <- function(object, order=NULL, plot=TRUE){
  
  fmatlong <- cbind(expand.grid("feature" = rownames(object),
                                "sample" = colnames(object),
                                stringsAsFactors=FALSE),
                    "intensity" = as.vector(object))
  
  if (!is.null(order))
    fmatlong$feature <- factor(fmatlong$feature, order)
  
  fmatlong$sample <- factor(fmatlong$sample, colnames(object))
  
  p <- ggplot(fmatlong, aes(sample, feature, fill = intensity)) +
    geom_tile() +
    scale_fill_continuous(low="white", high="#56B4E9", limits=c(0,NA), name="Intensity") +
    theme_bw() +
    xlab("") +
    ylab("") +
    theme(panel.grid=element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(),
          axis.ticks=element_blank(),
          axis.text.x=element_text(angle=45, vjust=1, hjust=1),
          aspect.ratio=1) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))
  
  if (plot)
    print(p)
  
  invisible(p)
}

plotConsProfileAll <- function(obj, fcol="markers", label_all_x=FALSE, foi=NULL, foi_name=NULL){
  
  if(!is.null(foi)){
    if(is.null(foi_name)){
      stop("need to provide foi_name")
    }
    fData(obj)[rownames(obj) %in% foi, fcol] <- foi_name
    fData(obj)[[fcol]] <- factor(fData(obj)[[fcol]])
  }
  
  x <- markerMSnSet(obj, fcol = fcol)
  fmat_all <- exprs(x)
  
  order <- NULL
  hc <- mrkHClust(obj, fcol, plot=FALSE)
  mm <- getMarkerClasses(x, fcol)
  m_order <- levels(factor(mm))[order.dendrogram(hc)]
  
  for (mc in m_order){
    fmat <- exprs(x[fData(x)[[fcol]]==mc,])
    if(nrow(fmat)>1){
      #umm <- levels(factor(rownames(fmat)))
      x.dist <- dist(fmat)
      hc <- hclust(x.dist)
      #hc <- hclust(as.dist(1-abs(cor(t(fmat)))))
      hc <- as.dendrogram(hc)
      order <- append(order, rownames(fmat)[order.dendrogram(hc)])
    }
    else{
      order <- append(order, rownames(fmat))
    }
  }
  
  fmatlong <- cbind(expand.grid("feature" = rownames(fmat_all),
                                "sample" = colnames(fmat_all),
                                stringsAsFactors=FALSE),
                    "intensity" = as.vector(fmat_all))
  mm <- getMarkers(x, fcol, verbose=FALSE)
  fmatlong$marker_class <- mm[match(fmatlong$feature, names(mm))]
  fmatlong$marker_class <- factor(fmatlong$marker_class, levels=rev(m_order))
  fmatlong$feature <- factor(fmatlong$feature, levels=order)
  fmatlong$sample <- factor(fmatlong$sample, colnames(fmat))
  
  if(label_all_x){
    labels <- as.character(unique(fmatlong$sample))
    }
  else{
    labels <- as.character(unique(fmatlong$sample)[seq(0, length(unique(fmatlong$sample)), by=2)])
  }
  
  p <- ggplot(fmatlong, aes(sample, feature, fill = intensity)) +
    geom_tile() +
    scale_fill_continuous(low="white", high="#56B4E9", limits=c(0,NA), name="Intensity") +
    theme_bw() +
    xlab("") +
    ylab("") +
    theme(panel.grid=element_blank(),
          panel.border=element_rect(colour="grey70"),
          axis.line.x = element_line(),
          axis.ticks=element_blank(),
          axis.text.x=element_text(angle=45, vjust=1, hjust=1),
          axis.text.y=element_blank(),
          strip.text.y=element_text(size=10, angle=180, hjust=1),
          strip.background=element_blank(),
          panel.spacing = unit(0, "lines"),
          aspect.ratio=0.1) +
    scale_x_discrete(expand=c(0,0), breaks=labels) +
    scale_y_discrete(expand=c(0,0)) +
    facet_grid(marker_class~., scales="free", switch="y")
  
  return(p)
  
}



getPhenoMarkers <- function(obj, fcol){
  m <- getMarkers(obj, fcol, verbose=FALSE)
  pm <- m[grep("Phenotype ", m)]
  pm <- table(pm)
  invisible(pm)
}

getPhenoConvergence <- function(obj){
  last_n_pc <- 0
  last_pm <- NULL
  rows <- NULL
  for (iteration in fvarLabels(obj)[grep(".pd", fvarLabels(obj))]){
    pm <- getPhenoMarkers(obj, iteration)
    total_pm <- sum(pm) 
    n_pc <- length(pm)

    if (n_pc == last_n_pc){
      distance <- sqrt(sum((as.matrix(pm)-as.matrix(last_pm))^2))
    }
    else{
      distance <- NA
    }
    
    if(iteration==".pd"){
      iteration_n <- "final"
    }
    else{
      iteration_n <- gsub(".pd", "", iteration)
    }
    rows[[iteration]] <- c(iteration_n, n_pc, total_pm, distance)
    
    last_pm <- pm
    last_n_pc <- n_pc
  }
  
  pheno_convergence <- data.frame(do.call("rbind", rows), stringsAsFactors=FALSE)
  colnames(pheno_convergence) <- c("iteration", "phenotypes", "proteins_in_phenotypes", "distance")
  pheno_convergence <- pheno_convergence %>% mutate_all(as.numeric)
  
  invisible(pheno_convergence)
}


plotPhenoConvergence <- function(obj){
  p <- ggplot(obj, aes(iteration, phenotypes)) + geom_point() +
    my_theme
  
  p2 <- p + aes(y=proteins_in_phenotypes)
  
  p3 <- p + aes(y=distance)
  
  invisible(list("p1"=p, "p2"=p2, "p3"=p3))
}

plotTempPhenoDisco <- function(infile){
  load(infile)
  top_iteration <- rev(fvarLabels(object)[grep(".pd", fvarLabels(object))])[1]
  pheno_convergence <- getPhenoConvergence(object)
  pheno_convergence_plots <- plotPhenoConvergence(pheno_convergence)
  print(pheno_convergence_plots)
  
  fData(object)$.pd_cleaned <- fData(object)[[top_iteration]]
  fData(object)$.pd_cleaned[!grepl("Phenotype|GOLGI", fData(object)[[top_iteration]])] <- "unknown"
  print(plotPCA(prot_res_imputed_w_m, "markers"))
  print(plotPCA(object, ".pd_cleaned"))
  pca_plot <- plotPCA(object, top_iteration)
  print(pca_plot)

  print(plotPCA(prot_res_imputed_w_m, "markers", dims=c(3,4)))
  print(plotPCA(object, ".pd_cleaned", dims=c(3,4)))
  
  cons_plot <- plotConsProfileAll(object, top_iteration)
  print(cons_plot)
  
  invisible(list("convergence"=pheno_convergence_plots, "pca_plot"=pca_plot, "cons_plot"=cons_plot,
                 "object"=object, "top_iteration"=top_iteration))
}

runTAGMMAP <- function(obj, min_prob_threshold=0.9, max_outlier_prob=0.1){
  
  print(table(fData(obj)$markers))
  params <- tagmMapTrain(obj)
  
  plotEllipse(obj, params)
  addLegend(object = obj, where = "bottomleft", cex = 1)
  
  plotEllipse(obj, params, dims=c(3,4))
  addLegend(object = obj, where = "topright", cex = 1)
  
  res_map <- tagmPredict(obj, params)
  
  print(table(fData(res_map)$tagm.map.allocation, fData(res_map)[["markers"]]))
  print(table(fData(res_map)$tagm.map.allocation, fData(res_map)$tagm.map.probability>min_prob_threshold))
  print(table(fData(res_map)$tagm.map.allocation, fData(res_map)$tagm.map.outlier>max_outlier_prob))
  
  fData(res_map)$tagm.map.allocation_filtered <- fData(res_map)$tagm.map.allocation
  fData(res_map)$tagm.map.allocation_filtered[fData(res_map)$tagm.map.probability<min_prob_threshold] <- "unknown"
  fData(res_map)$tagm.map.allocation_filtered[fData(res_map)$tagm.map.outlier>max_outlier_prob] <- "unknown"
  
  p <- fData(res_map) %>% data.frame() %>% ggplot() +
    aes(tagm.map.probability) + geom_histogram() + my_theme +
    facet_wrap(~tagm.map.allocation, scales="free") +
    theme(axis.text.x=element_text(angle=90)) +
    xlim(0,NA)
  
  print(p)
  
  ptsze <- (exp(fData(res_map)$tagm.map.probability) - 1)/2
  ptsze <- fData(res_map)$tagm.map.probability
  
  plot2D(res_map,
         fcol = "tagm.map.allocation",
         cex = ptsze,
         main = "protein pointer scaled with posterior localisation probability")
  
  addLegend(object = res_map, where = "bottomright", cex = 1)
  
  print(plotPCA(res_map, "tagm.map.allocation_filtered"))
  print(plotPCA(res_map, "tagm.map.allocation_filtered", dims=c(3,4)))
  
  print(PlotMarkerProfiles(res_map, "markers", keep_markers=getMarkerClasses(res_map)))
  print(PlotMarkerProfiles(res_map, "tagm.map.allocation", keep_markers=getMarkerClasses(res_map)))
  print(PlotMarkerProfiles(res_map, "tagm.map.allocation_filtered", keep_markers=getMarkerClasses(res_map)))
  
  invisible(res_map)
}

combineLOPITs <- function(lopit1, lopit2){
  
  tmp_exprs <- rbind(exprs(lopit1), exprs(lopit2))
  
  tmp_fData1 <- fData(lopit1)
  tmp_fData2 <- fData(lopit2)
  
  for(col in setdiff(colnames(tmp_fData2), colnames(tmp_fData1))){
    tmp_fData1[[col]] <- NA
  }
  
  for(col in setdiff(colnames(tmp_fData1), colnames(tmp_fData2))){
    tmp_fData2[[col]] <- NA
  }
  
  tmp_fData <- rbind(tmp_fData1, tmp_fData2)
  
  mod_prot_plus_mod_pep <- MSnSet(tmp_exprs, tmp_fData, pData(lopit1))
  
  invisible(mod_prot_plus_mod_pep)
}

combineUnmodProtModPep <- function(unmod_prot, mod_pep, mod_to_unknown=FALSE){
  
  combined_data <- combineLOPITs(unmod_prot, mod_pep)
  fData(combined_data)$mod <- c(rep("Unmodified", length(rownames(unmod_prot))),
                                rep("Modified", length(rownames(mod_pep))))
  
  for(col in fvarLabels(combined_data)[grep("^CV", fvarLabels(combined_data))])
    
    if(mod_to_unknown){
      # reset the marker column to unknown so that all modified features can be classified  
      fData(combined_data)$markers[fData(combined_data)$mod=="Modified"] <- "unknown"
    }
  
  invisible(combined_data)
}
