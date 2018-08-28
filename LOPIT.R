# emulate standard pROLOC colours

suppressMessages(library("pRoloc"))
suppressMessages(library("reshape2"))
suppressMessages(library("ggplot2"))
suppressMessages(require(gridExtra))

setStockcol(NULL) ## reset first
setStockcol(paste0(getStockcol(), 70))

getClassColours <- function(){
  cols = getStockcol()[c(1:5,12,7:11,13:20)]
  return(cols)
  
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
  print(p)
}

# ------------------------------------------------------------------------------------------------------------
# Function	: PlotMarkerProfiles 
# ------------------------------------------------------------------------------------------------------------

PlotMarkerProfiles <- function(res, fcol="master_protein", keep_markers=c("CYTOSOL"), unknown=F, plot_all=T,
                               individual_plots=F, foi=NULL){
  
  
  exprs_df <- exprs(res)
  f_df <- fData(res)
  exprs_df <- melt(exprs_df)
  exprs_df <- merge(exprs_df, f_df['markers'], by.x="Var1", by.y="row.names")
  
  if(unknown){
    exprs_trim_df <- exprs_df[exprs_df$markers %in% c(keep_markers, "unknown"),]
    exprs_trim_df$unknown <- exprs_trim_df$markers=="unknown"
    exprs_trim_df$markers[exprs_trim_df$markers=="unknown"] <- NA
  }
  else{
    exprs_trim_df <- exprs_df[exprs_df$markers %in% keep_markers,]
  }
  
  p <- ggplot(exprs_trim_df, aes(Var2, value, colour=markers)) +
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
      p <- p + geom_line(alpha=0.25)
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
      p <- p + stat_summary(aes(group=markers), geom="line", fun.y=mean) +
        scale_alpha_manual(values=c(1, 0.01), guide=F)
    }
    else{
      p <- p + stat_summary(aes(group=markers), geom="line", fun.y=mean)
    }
  }
  
  if(!missing(foi)){
    exprs_df <- data.frame(exprs_df)
    features_df <- exprs_df[exprs_df$Var1 %in% foi,]
    p <- p + geom_line(data=features_df, aes(Var2, value), colour="black", alpha=0.3)
    
  }
  
  #print(p)
  return(p)
}

# ------------------------------------------------------------------------------------------------------------
# Function	: makeQSepDistance
# ------------------------------------------------------------------------------------------------------------

makeQSepDistance <- function(res_with_markers, plot_cluster=F){
  qsep_dist <- qsep(QSep(res_with_markers))
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
    scale_fill_continuous(name="Distance")
  
  return(p)
  
}

# ------------------------------------------------------------------------------------------------------------
# Function	: makeQSepDistancePlot 
# ------------------------------------------------------------------------------------------------------------

plotCompareDistances <- function(qsep_df1, qsep_df2, exp1="All", exp2="TMT2", qsep_df3=F, exp3=F,
                                 membrane_organelles=FALSE, facet_out=F, plot_den=F){
  
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
# Function	: plotPCA
# ------------------------------------------------------------------------------------------------------------
makePCA_proj <- function(res, fcol="markers", dims=c(1,2)){
  PCA_matrix <- plot2D(res, fcol=fcol, plot=FALSE, dims=dims)
  axes <- colnames(PCA_matrix)
  colnames(PCA_matrix) <- c("X", "Y")
  PCA_df <- merge(data.frame(PCA_matrix), data.frame(fData(res)), by="row.names")
  PCA_df$markers <- PCA_df[[fcol]]
  PCA_df$unknown <- PCA_df$markers=="unknown"
  PCA_df$markers[PCA_df$unknown] <- NA
  PCA_df <- PCA_df[order(-PCA_df$unknown),]
  return(list("PCA_df"=PCA_df, "axes"=axes))
}

LOPITPCAPlotter <- function(PCA_df, axes, xlims=FALSE, ylims=FALSE, foi=FALSE, add_foi_names=FALSE,
                    just_markers=FALSE, m_colours=FALSE, re_order_markers=FALSE, marker_levels=NULL,
                    title=FALSE, fcol=FALSE, point_size=FALSE, foi_size=1, foi_alpha=1, foi_shape=8,
                    foi_colour="black", foi_fill="black", return_proj=FALSE){
  
  if(!missing(fcol)){
    # remake marker column
    PCA_df$markers <- PCA_df[[fcol]]
    PCA_df$markers <- PCA_df[[fcol]]
    PCA_df$unknown <- PCA_df$markers=="unknown"
    PCA_df$markers[PCA_df$unknown] <- NA
    PCA_df <- PCA_df[order(-PCA_df$unknown),]
  }
  
  if(re_order_markers!=FALSE){
    PCA_df$markers <- factor(PCA_df$markers, levels=marker_levels)
  }
  
  plots = NULL
  if(m_colours==FALSE){
    m_colours <- getStockcol()[c(1:5,12,7:11,13:20)]
  }
  p <- ggplot(PCA_df) +
    geom_point(aes(X, Y, colour=markers, alpha=unknown)) +
    scale_alpha_manual(values=c(0.5,0.1), guide=F) +
    #scale_colour_manual(values=rep("#E41A1C70", 12), na.value="grey95") +
    scale_colour_manual(values=m_colours, na.value="grey80", name="") +
    my_theme + xlab(axes[1]) + ylab(axes[2])
    
  if(!missing(point_size)){
    p <- p + aes_string(size=point_size) +
      scale_size_area(guide=F, max_size=2)
    print(p)
  }
  
  if(!missing(title)){
    p <- p + ggtitle(title)
  }
  
  if(xlims){
    p <- p + xlim(xlims)
  }
  
  if(ylims){
    p <- p + ylim(ylims)
  }
  
  plots[['p']] <- p + guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
  if(just_markers){
    PCA_df2 <-PCA_df[!is.na(PCA_df$markers),]
    p2  <- p %+% PCA_df2
    plots[['p2']] <- p2 + guides(colour = guide_legend(override.aes = list(alpha = 1)))
    
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
    
    plots[['p_foi']] <- p + guides(colour = guide_legend(override.aes = list(alpha = 1)))
    
    if(just_markers){
      features_df2 <- PCA_df2[PCA_df2[['Row.names']] %in% foi,]
      p2 <- p2 + geom_point(data=features_df2, aes(X,y), shape=8)
      
      if(!missing(point_size)){
        p2 <- p2 + aes_string(size=point_size)
      }
      
      plots[['p2_foi']] <- p2 + guides(colour = guide_legend(override.aes = list(alpha = 1)))
    }
  }
  
  return(plots)
}

plotPCA <- function(res, fcol, dims=c(1,2), ...){
  
    PCA <- makePCA_proj(res, fcol, dims)
    LOPITPCAPlotter(PCA$PCA_df, PCA$axes, ...)
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
    params <- ml_method(res_with_markers, fcol="markers", times=times,
                        verbose=T, ...)
  }
  else{
    params <- ml_method(res_with_markers, fcol="markers", times=times,
                        verbose=T, class.weights=classWeights(res_with_markers), ...)
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
  colnames(OrganelleF1) <- getMarkerClasses(res_with_markers)
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

plotPC1to_4 <- function(res, just_markers=F){
  p1 <- plotPCA(res, just_markers=just_markers)
  p2 <- plotPCA(res, dims=c(3,4), just_markers=just_markers)
  if(just_markers){
    grid.arrange(p1$p2, p2$p2, ncol=1)}
  else{
    grid.arrange(p1$p, p2$p, ncol=1)
  }
}

plotPC1to_8 <- function(res, just_markers=F){
  p1 <- plotPCA(res, just_markers=just_markers)
  p2 <- plotPCA(res, dims=c(3,4), just_markers=just_markers)
  p3 <- plotPCA(res, dims=c(5,6), just_markers=just_markers)
  p4 <- plotPCA(res, dims=c(7,8), just_markers=just_markers)
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


