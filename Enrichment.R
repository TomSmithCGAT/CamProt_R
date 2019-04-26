suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(goseq))
suppressMessages(library(Hmisc))
suppressMessages(library(mygene))
suppressMessages(library(data.table))

source("Utility.R")

# ------------------------------------------------------------------------------------------------------------
# Function	: myProtMapper 
# Aim		: To use the function 'queryMany' from Bioconductor package mygene as fast and most up-to-date
# Input 
# 	: ids = a character list of ids which can be uniprot, ensembl gene, gene symbol,etc
#       : id.type = what type of ids have you provided in the 'ids' list. Default = "uniprot"
#       : outlist = list of ids you want as an output. Default = c("interpro","ensembl.gene","go")
#       : modify = Logical, Default = T; Would you like to modify fields such as interpro, enseml, go to make them more human readable.
# Output  	: A dataframe with required ids and input ids 
# ------------------------------------------------------------------------------------------------------------------------

myProtMapper <- function(ids,id.type="uniprot",out.fields=c("interpro.short_desc","ensembl.gene","go.MF.id","go.CC.id","go.BP.id"),modify=T){
  
  # Get the mapping
  qm = queryMany(ids,scopes=id.type,fields=out.fields)
  
  # Returning variable
  ret.qm = NULL
  
  # Resolve the mappings to make them human readable
  if(modify == T){
    qm$go.all = NULL
    
    # Interpro mappings
    if(is.element("interpro",colnames(qm))){
      qm$domains = sapply(qm$interpro,function(x) paste(unlist(x),collapse=";"))
    }
    else{
      print("No Interpro domains")
    }
    
    # KEGG
    if(is.element("pathway.kegg",colnames(qm))){
      qm$kegg = sapply(qm$pathway.kegg,function(x) paste(unlist(x),collapse=";"))
    }
    else{
      print("No KEGG")
    }
    
    # GO mappings
    if(!is.na(grep("go",colnames(qm)))){
      
      # Grep all the go columns 'go.CC','go.MF','go.BP'
      f = grep("go",colnames(qm), value=T)
      
      qm$go.all = apply(qm[,f], MARGIN=1, FUN = function(x) paste0(as.character(unique(unlist(x))), collapse=";"))
      qm$go.all = gsub("^;","",gsub(";;",";",qm$go.all))
      qm$go.count = lengths(strsplit(qm$go.all,";"))
    }
    else{
      print("No GO terms")
    }
    
    # Ensembl.gene mappings
    if(is.element("ensembl",colnames(qm))){
      qm$ens = sapply(qm$ensembl,function(x) paste(unlist(x),collapse=";"))
    }
    else{
      print("No Ensembl Ids")
    }
    
    # Return mapped structure with tidy columns
    ret.qm = qm
  }
  else{
    ret.qm = qm
  }
  
  return(data.frame(ret.qm))
}



# ------------------------------------------------------------------------------------------------------------
# Function	: myProtMapper 
# Aim		: To use the function 'queryMany' from Bioconductor package mygene as fast and most up-to-date
# Input 
# 	: ids = a character list of ids which can be uniprot, ensembl gene, gene symbol,etc
#       : id.type = what type of ids have you provided in the 'ids' list. Default = "uniprot"
#       : outlist = list of ids you want as an output. Default = c("interpro","ensembl.gene","go")
#       : modify = Logical, Default = T; Would you like to modify fields such as interpro, enseml, go to make them more human readable.
# Output  	: A dataframe with required ids and input ids 
# --
cleanQM <- function(qm){
  duplicate_query <- names(table(qm$query)[table(qm$query)>1])
  duplicate_ens <- names(table(qm$ens)[table(qm$ens)>1])
  multi_ens <- qm$ens[grep(";", qm$ens)]
  
  qm <- qm[!qm$query %in% duplicate_query,]
  qm <- qm[!qm$ens %in% duplicate_ens,]
  qm <- qm[!qm$ens %in% multi_ens,]
  
  return(qm)
}

# ------------------------------------------------------------------------------------------------------------
# Function	: plotGOTerms
# ------------------------------------------------------------------------------------------------------------
# 15.9.17 Added 'term.col' as parameter as in Interpro enrichment, there is no term, just 'category'

plotGOTerms <- function(GO_df, len_foreground=NULL, len_background=NULL,
                        BH_filter=0.01, enrichment_filter=2, numObs_filter=50,
                        switch_axes=F, plot_top=10, term.col="term", cut_at=30,
                        order_by_enrichment=T, enrichment_col=NULL){
  

  GO_df <- GO_df %>% filter(!is.na(GO_df$term))
  # GO_df = dataframe from call to goseq
  # len_foreground = total number of proteins in foreground
  # len_background = total number of proteins in background
  
  GO_df$BH <-  p.adjust(GO_df$over_represented_pvalue, method="BH")
  GO_df$BH[GO_df$BH==0] <- 1E-16
  
  if(missing(enrichment_col)){
    if(missing(len_foreground) | missing(len_background)){
      stop("need to provide lengths of foregound and background")
    }
    GO_df$enrichment <- log(((GO_df$numDEInCat/GO_df$numInCat) / (len_foreground/len_background)),2)
  }
  else{
    GO_df$enrichment <- GO_df[[enrichment_col]]
  }
  
  print(head(GO_df))
  
  GO_filtered_df <- GO_df[GO_df$BH < BH_filter,]
  GO_filtered_df <- GO_filtered_df[GO_filtered_df$enrichment > enrichment_filter,]
  GO_filtered_df <- GO_filtered_df[GO_filtered_df$numDEInCat > numObs_filter,]
  
  
  
  if(switch_axes){ # determine how long to allow GO term desc. before truncating
    if(order_by_enrichment){
      GO_filtered_df <- GO_filtered_df[order(GO_filtered_df$enrichment),]}
    else{
      GO_filtered_df <- GO_filtered_df[order(-GO_filtered_df$BH),]}
    if(missing(cut_at)){
      cut_at <- 40
    }
  }
  else{
    if(order_by_enrichment){
      GO_filtered_df <- GO_filtered_df[order(-GO_filtered_df$enrichment),]}
    else{
      GO_filtered_df <- GO_filtered_df[order(GO_filtered_df$BH),]}
    if(missing(cut_at)){
      cut_at <- 25
    }
  }
  
  GO_filtered_df$term <- factor(GO_filtered_df[,term.col], levels=unique(GO_filtered_df[,term.col]))
  
  
  plots = NULL

  if(length(GO_filtered_df$ontology)<=plot_top){
    GO_filtered_df <- GO_filtered_df[order(GO_filtered_df$ontology),]
    p <- ggplot(GO_filtered_df, aes(interaction(term, ontology), enrichment, fill=log(BH,10))) +
      geom_bar(stat="identity") + 
      xlab("") + ylab("Over-representation (Log2)") +
      scale_fill_continuous(name="BH adj. p-value\n(Log 10)\n", low=cbPalette[3], high="grey30") +
      my_theme +
      theme(
        text=element_text(size=15),
        plot.title=element_text(hjust=0.5)) +
      scale_x_discrete(labels=paste0(sapply(as.character(GO_filtered_df[[term.col]]), FUN=function(x) shortenTerm(x, cut_at)),
                                     " (", GO_filtered_df$ontology, ")"))
    
    if(switch_axes){
      p <- p + coord_flip()
    }
    else{
      p <- p + theme(axis.text.x=element_text(size=12, angle=45, vjust=1, hjust=1))
    }
    
    plots[["all"]] <- p
  }
  else{
    
    onto_class2name = c("Biological Process", "Cellular Compartment", "Molecular Function")
    names(onto_class2name) <- c("BP", "CC", "MF")
    
    for(onto_class in unique(GO_filtered_df$ontology)){
      if(is.na(onto_class)){
        next
      }

      tmp_df <- GO_filtered_df[GO_filtered_df$ontology == onto_class,]
      
      if(switch_axes){
        tmp_df <- tail(tmp_df,min(plot_top, length(tmp_df[,1])))
      }
      else{
        tmp_df <- head(tmp_df,min(plot_top, length(tmp_df[,1])))
      }

      tmp_df$short_term <- sapply(as.character(tmp_df$term), FUN=function(x) shortenTerm(x, cut_at))
      
      p <- ggplot(tmp_df, aes(term, enrichment, fill=log(BH,10))) +
        geom_bar(stat="identity") + 
        xlab("") + ylab("Over-representation (Log2)") +
        scale_fill_continuous(name="BH adj. p-value\n(Log 10)\n", low=cbPalette[3], high="grey30") +
        ggtitle(onto_class2name[onto_class]) +
        my_theme +
        theme(
          text=element_text(size=15),
          plot.title=element_text(hjust=0.5)) +
        scale_x_discrete(labels=tmp_df$short_term)
      
      if(switch_axes){
        p <- p + coord_flip() + theme(aspect.ratio=2, axis.text.y=element_text(size=15), axis.text.x=element_text(size=15))
      }
      else{
        p <- p + theme(axis.text.x=element_text(size=12, angle=45, vjust=1, hjust=1))
      }
      plots[[onto_class]] <- p
    }
  }
  return(list(p=plots,g=GO_df))
}

# ------------------------------------------------------------------------------------------------------------
# Function	: plotInterPro
# ------------------------------------------------------------------------------------------------------------
# 15.9.17 Added 'term.col' as parameter as in Interpro enrichment, there is no term, just 'category'

plotInterPro <- function(GO_df, len_foreground, len_background,
                         BH_filter=0.01, enrichment_filter=2, numObs_filter=20,
                         plot_top=10, term.col="category"){
  
  # GO_df = dataframe from call to goseq
  # len_foreground = total number of proteins in foreground
  # len_background = total number of proteins in background
  
  GO_df$BH <-  p.adjust(GO_df$over_represented_pvalue, method="BH")
  GO_df$enrichment <- log((GO_df$numDEInCat/GO_df$numInCat) / (GO_df$numInCat/len_background),2)
  GO_df <- GO_df[order(GO_df$enrichment,decreasing=T),]
  
  GO_filtered_df <- GO_df[GO_df$BH < BH_filter,]
  GO_filtered_df <- GO_filtered_df[GO_filtered_df$enrichment > enrichment_filter,]
  GO_filtered_df <- GO_filtered_df[GO_filtered_df$numDEInCat > numObs_filter,]
  GO_filtered_df <- head(GO_filtered_df,min(plot_top, length(GO_filtered_df[,1])))
  
  
  # ggplot2 function being called
  p <- ggplot(GO_filtered_df, aes(GO_filtered_df[,term.col], enrichment, fill=log(BH,10))) +
    geom_bar(stat="identity") + 
    xlab("") + ylab("Over-representation (Log2)") +
    scale_fill_continuous(name="BH adj. p-value\n(Log 10)\n") +
    my_theme +
    theme(
      text=element_text(size=15),
      plot.title=element_text(hjust=0.5)) +
    scale_x_discrete(labels=GO_filtered_df[,term.col],limits=GO_filtered_df[,term.col])+ coord_flip()
  
  return(list(p=p,g=GO_df))
}


plotKEGGterms <- function(kegg_over_rep_obj, kegg_terms=NULL, facet_by=NULL, BH_filter=0.1){
  
  kegg_over_rep_obj$short_cat <- sapply(strsplit(kegg_over_rep_obj$category, " - Homo sapiens"), "[[", 1)
  
  if(missing(kegg_terms)){
    tmp_df <- kegg_over_rep_obj[kegg_over_rep_obj$BH<BH_filter,]  
  } else{
    tmp_df <- kegg_over_rep_obj[kegg_over_rep_obj$category %in% kegg_terms,]  
    print(head(tmp_df))
  }
  
  tmp_df <- tmp_df[order(tmp_df$adj_over_rep),]
  tmp_df$short_cat <- factor(tmp_df$short_cat, levels=unique(tmp_df$short_cat))
  
  if(!missing(facet_by)){
    tmp_df$facet = tmp_df[[facet_by]]
  }
  
  p <- ggplot(tmp_df, aes(short_cat, adj_over_rep, fill=-log10(BH))) +
    geom_bar(stat="identity") +
    my_theme +
    theme() +
    coord_flip() +
    scale_fill_continuous(name="-log10 BH-adjusted p-value",
                          guide = guide_colorbar(title.position="top"),
                          high=cbPalette[3], low="grey30") +
    theme(text=element_text(size=15),
          legend.position="top",
          aspect.ratio=1.5) +
    ylab("Log2-fold Over-representation") +
    xlab("")
  
  if(!missing(facet_by)){
    p <- p + facet_grid(facet~., scales="free")
  }
  
  print(p)
  return(p)
}


makePWF <- function(df, sig_col, bias_col, identifier_col){
  bias <- df[[bias_col]]
  names(bias) <- df[[identifier_col]]
  
  sig_genes <- df[[sig_col]]==TRUE
  
  names(sig_genes) <- df[[identifier_col]]
  print(sum(sig_genes))
  
  pwf.counts=nullp(sig_genes, bias.data=bias, plot.fit = F)
  return(pwf.counts)
}


makePWFPlot <- function(pwf.counts, bias.data_col='bias.data', bins=20, xlab, spline=T){
  
  pwf.counts <- pwf.counts[order(pwf.counts[[bias.data_col]]),]
  pwf.counts$DEgenes <- as.numeric(pwf.counts$DEgenes)
  
  if(bins>0){
    print("not missing bins")
    pwf.counts$bin <- as.numeric(as.character(cut2(pwf.counts[[bias.data_col]], g = bins, levels.mean = T)))
  }
  else{
    pwf.counts$bin <- as.numeric(as.character(pwf.counts[[bias.data_col]]))
  }

  p <- ggplot(pwf.counts, aes(bin, DEgenes)) +
    stat_summary(fun.y=mean, geom="point")+ 
    my_theme +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
    xlab(xlab) +
    ylab("Prob. identified")
  
  if(spline){
    p <- p + geom_line(aes(x=bias.data, y=pwf), colour="green", linetype=2)
  }
  
  return(p)
  
}

getEnrichedGO <- function(pwf, gene2cat=FALSE, shorten_term=TRUE, ...){
  
  if(!missing(gene2cat)){
    GO.abundance<-goseq(pwf,gene2cat=gene2cat, ...)
  }
  else{
    GO.abundance<-goseq(pwf,"hg38","ensGene", ...)
  }
  
  GO.abundance$BH <-  p.adjust(GO.abundance$over_represented_pvalue, method="BH")
  if(shorten_term){GO.abundance$short_term <- substr(GO.abundance$term, 1, 30)}
  return(GO.abundance)
}

# ------------------------------------------------------------------------------------------------------------
# Function  : 'addAdjustedOverRep' A crude function to add an adjusted estimate of the over-representation of a term
# Input 
#           : obj = A data frame with the results from goseq. As generated by GetEnrichedGO
#           : pwf = a PWF from goseq. As generated by makePWF
#           : gene2cat = A dataframe mapping features to categories
#           : term_col = The name of the column with the "terms" of interest in gene2cat (e.g GO terms)
#           : target_id = THe name of the features in gene2cat, e.g "gene_id"
# Output  : The input obj + a column with estimated adjusted over-representation for each term ($adj_over_rep)
# ------------------------------------------------------------------------------------------------------------------------

addAdjustedOverRep <- function(obj, pwf, gene2cat, term_col="transcript_biotype", target_col="transcript_id"){
  len_fore <- sum(pwf$DEgenes)
  len_back <- length(pwf$DEgenes)

  obj$adj_over_rep <- apply(obj[,c("numDEInCat", "numInCat", "category")], MARGIN=1, function(x){
    term_features <- gene2cat[gene2cat[[term_col]]==x[["category"]], target_col]
    term_weight <- mean(pwf[rownames(pwf) %in% term_features, "pwf"])
    non_term_weight <- mean(pwf[!rownames(pwf) %in% term_features, "pwf"])
    as.numeric(x[['numDEInCat']])/as.numeric(x[['numInCat']]) / (term_weight/non_term_weight) / (len_fore/len_back)})
  
  return(obj)
}

# ------------------------------------------------------------------------------------------------------------
# Function	: plotGOTerms
# ------------------------------------------------------------------------------------------------------------
plotTerms <- function(terms_df, pwf, gene2cat,
                      term_col="transcript_biotype",
                      target_col="transcript_id",
                      BH_filter=0.01,
                      overrep_filter=1, numObs_filter=50,
                      switch_axes=F, plot_top=10){
  terms_df$BH <-  p.adjust(terms_df$over_represented_pvalue, method="BH")
  terms_df$BH[terms_df$BH==0] <- 1E-16

  terms_filtered_df <- terms_df[terms_df$BH <= BH_filter,]
  terms_filtered_df <- terms_filtered_df[terms_filtered_df$numDEInCat > numObs_filter,]
  
  terms_filtered_df <- addAdjustedOverRep(terms_filtered_df, pwf, gene2cat, term_col, target_col)
  terms_filtered_df <- terms_filtered_df[terms_filtered_df$adj_over_rep > overrep_filter,]
  
  
  if(switch_axes){
    terms_filtered_df <- terms_filtered_df[order(terms_filtered_df$adj_over_rep),]
  }
  else{
    terms_filtered_df <- terms_filtered_df[order(-terms_filtered_df$adj_over_rep),]
  }
  
  p <- terms_filtered_df %>% head(plot_top) %>%
    ggplot(aes(category, log(adj_over_rep,2), fill=log(BH,10))) +
    geom_bar(stat="identity") + 
    xlab("") + ylab("Over-representation (Log2)") +
    scale_fill_continuous(name="BH adj. p-value\n(Log 10)\n", low=cbPalette[3], high="grey30") +
    my_theme +
    theme(
      text=element_text(size=15),
      plot.title=element_text(hjust=0.5))
    
  if(switch_axes){
    p <- p + coord_flip()
  }
  else{
    p <- p + theme(axis.text.x=element_text(size=12, angle=45, vjust=1, hjust=1))
  }
  
  return(list(plot=p,g=terms_filtered_df))
}
# ------------------------------------------------------------------------------------------------------------
# Function  : 'makeGene2Cat' to produce a 1:1 mapping of uniprot/ensembl/symbols to GO/Interpro terms. 
#              Will be used as input into the 'goseq' function in the gene2cat slot
# Input 
#           : dat = dataframe with ids and go/interpro terms (obtained from myProtMapper)
#           : from.id =  ids you want to map 'from'. Default = "uniprot"
#           : to.id =  ids you want to map to c("interpro","ensembl.gene","go")
#           : splt = symbol you want to split by if there are multiple ids
# Output  : A two column dataframe with Uniprot ids in the first and Go/Interpro in the second
# ------------------------------------------------------------------------------------------------------------------------

makeGene2Cat <- function(dat,from.id,to.id,splt){
  
  cat.frame = dat[,c(from.id,to.id)]
  d.dt = data.table(cat.frame,key=colnames(cat.frame[,from.id]))
  cat.out = data.frame(d.dt[, list(to.id = unlist(strsplit(get(to.id), splt))), by=from.id])
  cat.out = unique(cat.out)
  
  return(cat.out)
}

# ------------------------------------------------------------------------
# Function  : runGoseq
# Aim       : runs goseq on a list of genes
# Input     
#       : genelist = ids for genes of interest (default is Uniprot)
#       : bglist = background list against which enrichment is performed. Here, it is Geiger et al., 2012
#                  This is a dataframe that contains Unipt ids, gene symbols as well as bias data for the Geiger list of proteins
#       : bias  = column that contains bias information in the bglist object or a separate object with same nrows as bglist
#       : cat.oligo = mapping from uniprot ids to GO terms or Interpro terms to use in goseq
#       : signif = cut-off for calling significantly enriched genes after multiple-testing correction
# Output    : enriched list of Interpro domains
# ------------------------------------------------------------------------

runGoseq <- function(genelist,bglist,bias=NULL,cat.oligo=NULL,signif=0.05,method="Wallenius"){
  
  # setting up goseq object
  all.genes.comp = rep(0,nrow(bglist))
  names(all.genes.comp) = rownames(bglist)
  all.genes.comp[which(names(all.genes.comp) %in% unique(genelist))] = 1
  table(all.genes.comp)
  
  # Remove missing values
  comp.no.missing = all.genes.comp[which(!is.na(names(all.genes.comp)))]
  table(comp.no.missing)
  
  # Running the function to calculate weights - with user provided bias or default gene length bias
  pwf.comp = nullp(comp.no.missing,bias.data = bias)
  
  # Running goseq enrichment using default mapping within goseq or user provided mapper
  goseq.comp.cat = goseq(pwf.comp,gene2cat = cat.oligo, method=method)
  goseq.comp.cat$BH_over_represented_pvalue = p.adjust(goseq.comp.cat$over_represented_pvalue,method = "BH")
  enriched.goseq.comp.cat = goseq.comp.cat[which(goseq.comp.cat$BH_over_represented_pvalue <= signif),]
  
  return(list(goseq.comp.cat,enriched.goseq.comp.cat))
  
}


# ------------------------------------------------------------------------
# Function  : mymakespline
# Aim       : 
# Input     
#       : 
# Output    :
# ------------------------------------------------------------------------

mymakespline=function (x, y, newX=NULL, nKnots = 6, lower_bound=10^-3){
  # TS:this is a direct copy of the goseq function. Seems this function only works for
  # some input data when we generate x,y ourselves!
  
  #Should not be used outside of goseq package.  Refer to the help pages for pcls in the mgcv package for more general 
  #contstrained spline fitting.
  #We handle montonocily decreasing problems by reformulating them as monotonicly increasing by reflecting about "infinity"
  #Compare the first 10% to the last 10%
  ww=order(x)
  size=ceiling(length(y)/10) 
  low=sum(y[ww][1:size]) 
  hi=sum(y[ww][(length(y)-size):length(y)]) 
  #Is the trend decreasing
  if(hi<=low){
    #Reform as a montonicly increasing fit by reflecting about infitinity
    reflectionFactor=10^10
    x=reflectionFactor-x
  }
  #The number of knots to use when generating the spline
  nKnots <- round(nKnots)
  if (is.null(newX)) 
    newX <- x
  #We need to force the fit to go through (0,0), mono.con constructs a requirment that the LOWEST x VALUE be greater than 0
  #so to force the fit through (0,0) we add a fake lowest x value of 0 by adding in a dummy point. 
  #We use an inequality constraint rather than equality to only skew the fit enough that the spline is always >0
  #This won't work if we have a monotonically decreasing function (and we don't need to it), so check that
  if(hi>low){
    x=c(0,x)
    y=c(0,y)
  }
  
  f.ug <- gam(y ~ s(x, k = nKnots, bs = "cr"),family=binomial())
  dat <- data.frame(x = x, y = y)
  sm <- smoothCon(s(x, k = nKnots, bs = "cr"), dat, knots = NULL)[[1]]
  if (length(sm$xp) < 6) 
    warning("Few than 6 nKnots were specified!\n Possible inaccurate fitting!\n")
  F <- mono.con(sm$xp,TRUE)
  # The lower and upper parameters here don't seem to build the right constraints 100% of the time, so add them manually
  # ND-6/14- removing the upper limit constraint because it wasn't satisfied by initial parameters.
  # Hope it doesn't break the code.
  F$A=rbind(F$A,c(1,rep(0,ncol(F$A)-1))) #,c(rep(0,ncol(F$A)-1)),-1))
  F$b=c(F$b,0) #,1) ND - should this last one be -1?
  G <- list(X = sm$X, C = matrix(0, 0, 0), sp = f.ug$sp, p = sm$xp,
            y = y, w = y * 0 + 1, Ain =F$A, bin = F$b, S = sm$S, off = 0)
  # do some error checking as pcls can return NaNs sometimes if this isn't true
  if( any((G$Ain%*%G$p - G$bin) < 0 )){
    show(G$p)
    stop("Constraints for spline fit not satisfied by initial parameters")
  }
  #Do the actual fitting with penalized contstrained least squares
  p <- pcls(G)
  
  # Do some error checking for the rare case that p returns "Na"
  if(any(is.na(p))){
    warning("The splines fit failed! Returning a flat weight distribution. This may compromise your GO analysis, please check what happened.")
    return(rep(mean(y),length(newX)))
  }
  
  #Now what we want is the value of the spline at each data point x,
  fv <- Predict.matrix(sm, data.frame(x = newX)) %*% p
  fv <- as.vector(fv)
  
  #If pcls still fails for some reason and returns negative values (or values that are so low that they'll have an effective relative weight of Inf
  #then we need to add a little bit to every weight to make things non-negative, the choice of minimum weight is somewhat arbitrary and serves only to
  #ensure that we don't have positive weights that will be infinitely prefered as a ratio.
  if(min(fv)<lower_bound)
    fv=fv-min(fv)+lower_bound
  return(fv)
}





