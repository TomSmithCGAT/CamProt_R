#script.dir <- dirname(sys.frame(1)$ofile)
source("Utility.R")

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(VennDiagram))
suppressMessages(library(gplots))
suppressMessages(library(gridExtra))


# ----------------------------------------------------------------------------------------------------------------------
# Function	: parse_peptides
# Aim		: To aggregate peptides into master-proteins and carry out QC on dataset. 
# Input :
#    infile = file containing output from mascot server
#    silac = input is from a SILAC experiment
# Output  	: A dataframe with peptides associated to a unique master.protein, which is not cRAP or cRAP-associated. 
# ----------------------------------------------------------------------------------------------------------------------
parse_peptides <- function(infile, silac=FALSE, TMT=FALSE){
  
  print_n_pep <- function(peptide_df, message){
    cat(sprintf("%s\t%s\n", length(rownames(peptide_df)), message))
  }
  print_n_prot <- function(peptide_df){
    cat(sprintf("These peptides are associated with %s master proteins\n",
                length(unique(peptide_df$master_protein))))}
  
  print_summaries <- function(peptide_df, message){
    print_n_pep(peptide_df, message)
    print_n_prot(peptide_df)
  }
  
  peptide_df <- read.delim(infile,  header=T)
  cat("Tally of peptides at each stage:\n")
  
  print_summaries(peptide_df, "All peptides with a PSM")
  
  peptide_df <- peptide_df[peptide_df$master_protein!="",]
  print_summaries(peptide_df, "Excluding peptides without a master protein")
  
  peptide_df <- peptide_df[peptide_df$unique==1,]
  print_summaries(peptide_df, "Excluding peptides without a unique master protein")
  
  peptide_df <- peptide_df[peptide_df$crap_protein==0,]
  print_summaries(peptide_df, "Excluding peptides matching a cRAP protein")

  peptide_df <- peptide_df[peptide_df$associated_crap_protein==0,]
  print_summaries(peptide_df, "Excluding peptides associated with a cRAP protein")
  
  if (!TMT){
    peptide_df <- peptide_df[peptide_df$PSM.Ambiguity!="Rejected",]
    print_summaries(peptide_df, "Excluding peptides rejected due to PSM ambiguity")
  }
  
  if(silac|TMT){
    peptide_df <- peptide_df[peptide_df$Quan.Info!="NoQuanValues",]    
    print_summaries(peptide_df, "Excluding peptides without quantification")
  }
  return(peptide_df)
}


# ----------------------------------------------------------------------------------------------------------------------
# Function	: aggregate_modifications 
# Aim		: Remove duplicated (filename + Sequence + modification) rows since
#         these are just the same peptide but diff. modifications
# Input :
# 	: peptide_df = Peptide-level data.frame as obtained from parse_peptides above
# Output  	: A data.frame with duplicated rows removed
# ----------------------------------------------------------------------------------------------------------------------
aggregate_modifications <- function(peptide_df){
  return(peptide_df %>% distinct(Sequence, master_protein, filename, .keep_all=TRUE) %>% data.frame())}
    
# ----------------------------------------------------------------------------------------------------------------------
# Function	: aggregate_modifications 
# Aim		: Remove duplicated (filename + master_protein) rows since
#         these are just the same master_protein but different proteins
# Input :
# 	: peptide_df = Peptide-deduped-level data.frame as obtained from aggregate_modifications  above
# Output  	: A data.frame with duplicated rows removed
# ----------------------------------------------------------------------------------------------------------------------
aggregate_peptides <- function(peptide_agg_df){
  return(peptide_agg_df %>% distinct(master_protein, filename, .keep_all=TRUE) %>% data.frame())}


# ----------------------------------------------------------------------------------------------------------------------
# Function	: plot_pep_counts
# Aim		: Plot the peptide counts per protein, per sample
# Input :
# 	: peptide_df = Peptide-level data.frame as obtained from aggregate_modifications above
# Output  	: Plots
# ----------------------------------------------------------------------------------------------------------------------
plot_pep_counts <- function(peptide_df){
  
  pep_counts <- peptide_df %>%
    group_by(filename, master_protein) %>%
    dplyr::count() %>%
    separate(col = filename, into=c("species", "treatment", "rep"), sep="_", remove=FALSE)
  
  single_pep_hits <- pep_counts %>% filter(n==1) %>% ungroup() %>% dplyr::select(master_protein, treatment, rep)
  
  p1 <- ggplot(pep_counts, aes(log(n,2))) +
    geom_histogram() +
    xlab("Peptide counts (log2)") +
    ylab("Count") +
    my_theme +
    facet_grid(treatment~rep, scales='free')
  
  p2 <- ggplot(pep_counts, aes(x=rep, fill=n==1)) +
    geom_bar(position='fill') +
    xlab("Replicate") +
    ylab("Fraction of proteins") +
    scale_fill_manual(name="Single-hit-\nwonders", values=c("grey90", "grey10")) +
    my_theme +
    facet_wrap(~treatment)
  
  return(list("p1"=p1, "p2"=p2, "single_hits"=single_pep_hits))
}


# ----------------------------------------------------------------------------------------------------------------------
# Function	: tidy_protein
# Aim		: To extract information from filename to master_protein. 
# Input :
# 	: inputdata = file containing output from mascot server
#       : order = type of information in filename
# Output  	: A dataframe containing master_protein, treatment and replicate. 
# ----------------------------------------------------------------------------------------------------------------------

tidy <- function(inputdata, order=c("species", "treatment", "rep")){
  
  filtered_data <- data.frame(inputdata)
  filtered_data <- filtered_data[,c("master_protein", "filename")]
  filtered_data <- filtered_data %>% separate(col = filename, into=order, sep="_")
  filtered_data <- filtered_data[, c("master_protein", "treatment", "rep")]
  
}


get_prot_counts <- function(filtered_data){
  
  protein_count <- filtered_data %>%
    group_by(master_protein, treatment) %>%
    dplyr::count() %>%
    group_by(treatment, n) %>%
    dplyr::count()
  
  return(protein_count)
}

# ----------------------------------------------------------------------------------------------------------------------
# Function	: plotProtCounts
# Aim		: Plot the # times each protein was observed, split by NC vs. CL
# Input :
# 	: inputdata = filtered_data as output from tidy above
# Output  	: plot
# ----------------------------------------------------------------------------------------------------------------------

plotProtCounts <- function(filtered_data){
  
  protein_count <- get_prot_counts(filtered_data)
    
  p <- ggplot(protein_count, aes(n, nn)) + geom_bar(stat='identity') + 
    facet_wrap(~treatment) + my_theme +
    xlab("# times protein seen") +
    ylab("Protein count")
  
  return(p)
}

# ----------------------------------------------------------------------------------------------------------------------
# Function	: get_rep_counts
# Aim		: Get the counts of proteins for each possible combination of # times seen in NC or CL reps
# Input :
# 	: data = dataframe from tidy function
# Output  	:  data.frame of unique peptides and observations in each treatment
# ----------------------------------------------------------------------------------------------------------------------

get_rep_counts <- function(tidy_df){
  
  protein_counts <- tidy_df %>%
    group_by(master_protein, treatment) %>%
    count() %>%
    ungroup() %>%
    spread(key=treatment, value=n) %>%
    replace_na(list(CL=0)) 

}

# ----------------------------------------------------------------------------------------------------------------------
# Function	: compiler
# Aim		: To aggregate the data into NC and CL groups   
# Input :
# 	: data = dataframe from tidy function
# Output  	:  data.frame of unique peptides and observations in each treatment
# ----------------------------------------------------------------------------------------------------------------------

compiler <- function(data){
  
    dedup_data <- data[!duplicated(data),]
    dedup_data <- dedup_data %>% group_by(master_protein, treatment) %>%
    summarise('obs'=length(master_protein)) %>%
    spread(key = treatment, value = obs)
    dedup_data$NC[is.na(dedup_data$NC)] <- 0
    dedup_data$CL[is.na(dedup_data$CL)] <- 0
    
    
    return(dedup_data)
}

# ----------------------------------------------------------------------------------------------------------------------
# Function	: aggregation
# Aim		: To aggregate the data into number of times seen in CL 
# Input :
# 	: data = dataframe from compiler function
# Output  	:  vectors seen 1,2,3,4,5 times 
# ----------------------------------------------------------------------------------------------------------------------

aggregation <- function(compiler.data){
  
  a <- NULL
  b <- NULL
  c <- NULL
  d <- NULL
  e <- NULL
  for (CL in unique(compiler.data$CL)){
    tmp_df <- compiler.data[compiler.data$CL==CL & compiler.data$NC==0,]
    print(CL)
    print(length(unique(tmp_df$master_protein)))
    if(tmp_df$CL == 1){
      a <- append(a, unique(tmp_df$master_protein))
    }else if(tmp_df$CL == 2){
      b <- append(b, unique(tmp_df$master_protein))
    }else if(tmp_df$CL == 3){
      c <- append(c, unique(tmp_df$master_protein))
    }else if(tmp_df$CL == 4){
      d <- append(d, unique(tmp_df$master_protein))
    }else{
      e <- append(e, unique(tmp_df$master_protein))
    }
  }
  
frame <- list(a=a,b=b,c=c,d=d,e=e)
frame %>% lapply(length) %>% unlist %>% max -> mymax
frame %>% lapply(function(x) c(x, rep(NA, mymax-length(x))))
return(frame)
}


# ----------------------------------------------------------------------------------------------------------------------
# Function	: vennCL
# Aim		: To visualize intersect between CL and NC  
# Input :
# 	: inputdata = dataframe from tidy function
# Output  	:  venn diagram of CL versus NC
# ----------------------------------------------------------------------------------------------------------------------

vennCL <- function(data){
  
  grid.newpage()
  
  CL <- data[which(data$treatment == "CL"),]
  NC <- data[which(data$treatment == "NC"),]
  temp <- draw.pairwise.venn(length(unique(CL$master_protein)),
                     length(unique(NC$master_protein)), 
                     length(intersect(CL$master_protein, NC$master_protein)),
                     category = c("CL", "NC"), 	
                     fill = c("hotpink", "forestgreen"), alpha = rep(0.4, 2), 
                      lty = "blank", 
                      fontface = "bold",
                      fontfamily = "sans", 
                      cat.fontface = "bold",
                      cat.default.pos = "outer",
                      cat.fontfamily = "sans")
  grid.draw(temp)
}

# ----------------------------------------------------------------------------------------------------------------------
# Function	: histCL
# Aim		: Plot histogram of CL master_proteins in 1:5 replicates  
# Input : 
# 	: inputdata = data.frame from function tidy
#   : main = title of histogram
# Output  	:  histogram for replicates of CL fraction
# ----------------------------------------------------------------------------------------------------------------------

histCL <- function(data){

  data <- data[which(data$treatment == "CL"),]
  ONE <- NULL
  TWO <- NULL
  THREE <- NULL
  FOUR <- NULL
  FIVE <- NULL
  for (rep in unique(data$rep)){
    tempdf <- data[data$rep==rep,]
    print(rep)
    print(length(unique(tempdf$master_protein)))
    if(rep =="1"){
      ONE <- unique(as.character(append(ONE, unique(tempdf$master_protein))))}
    if(rep =="2"){
      TWO <- unique(as.character(append(TWO, unique(tempdf$master_protein))))}
    if(rep =="3"){
      THREE <- unique(as.character(append(THREE, unique(tempdf$master_protein))))}
    if(rep =="4"){
      FOUR <- unique(as.character(append(FOUR, unique(tempdf$master_protein))))}
    else{
      FIVE <- unique(as.character(append(FIVE, unique(tempdf$master_protein))))}
  }
  
  vect <- c(ONE, TWO, THREE, FOUR, FIVE)
  print(hist(table(vect)))

}

# ----------------------------------------------------------------------------------------------------------------------
# Function	: RBP
# Aim		:   Quantify number of RBP identified by both methods
# Input :
# 	: data.go = dataframe with GO term description
#   : data.uniprot = dataframe with uniprot accession numbers
# Output  	:  Print number of RBP identified by Beckmann and Manasa
# ----------------------------------------------------------------------------------------------------------------------

RBP <- function(data.go, data.uniprot){
    toMatch <- c("mRNA", "splic", "RNA binding", "RNA", "RNP", "translation", "ribosom", "nuclease", "exosome")
    RBP_Beckmann <- grep(paste(toMatch,collapse="|"), unique(data.go$GO), value = F)
    
    data <- as.data.frame(unique(data.uniprot$master_protein))
    colnames(data) <- "MP"
    suppressMessages(quer <- myProtMapper(data$MP))
    suppressMessages(uniprot_at.go <- makeGene2Cat(quer,"query","go.all",";"))
    RBP_Manasa <- uniprot_at.go %>% filter(uniprot_at.go$to.id =="GO:0003723")
    
    cat(sprintf("\n%s\tRBPs from Beckmann method\n", length(RBP_Beckmann)))
    cat(sprintf("%s\tRBPs from Manasa method\n", length(RBP_Manasa$query)))
    return(RBP_Manasa)
    
}
  
# ----------------------------------------------------------------------------------------------------------------------
# Function	: five.venn
# Aim		: generate a fivee.veen diagram between 5 setss   
# Input :
# 	: a-e = sets to be analyzed (must be an atomic vector)
#   : set_names = labels for the venn diagram 
# Output  	:  Print number of RBP identified by Beckmann and Manasa
# ----------------------------------------------------------------------------------------------------------------------


five.venn <- function(a, b, c, d, e, 
                      set1_name,
                      set2_name,
                      set3_name, 
                      set4_name,
                      set5_name){

  area1 = length(a)
  area2 = length(b)
  area3 = length(c)
  area4 = length(d)
  area5 = length(e)
  n12 = (intersect(a, b))
  n13 = (intersect((a),(c)))
  n14 = (intersect( (a),  (d)))
  n15 = (intersect( (a),  (e)))
  n23 = (intersect( (b),  (c)))
  n24 = (intersect( (b),  (d)))
  n25 = intersect(b, e)
  n34 = (intersect( (c),  (d)))
  n35 = (intersect( (c),  (e)))
  n45 = (intersect( (d),  (e)))
  n123 = (intersect(n12, c))
  n124 = intersect(n12, d)
  n125 = intersect(n12, e)
  n134 = intersect(n13, d)
  n135 = intersect(n13, e)
  n145 = intersect(n14, e)
  n234 = intersect(n23, d)
  n235 = intersect(n23, e)
  n245 = intersect(n24, e)
  n345 = intersect(n34, e)
  n1234 = intersect(n123, d)
  n1235 = intersect(n123, e)
  n1245 = intersect(n124, e)
  n1345 = intersect(n134, e)
  n2345 = intersect(n234, e)
  n12345 = intersect(n123, n45)
  
  
  draw.quintuple.venn(area1, area2, area3, area4, area5, length(n12), length(n13), length(n14), length(n15),
                      length(n23), length(n24), length(n25), length(n34), length(n35), length(n45), length(n123), 
                      length(n124), length(n125), length(n134), length(n135), length(n145), length(n234), 
                      length(n235), length(n245), length(n345), length(n1234), length(n1235),
                      length(n1245), length(n1345), length(n2345), length(n12345), 
                      category = c(set1_name, set2_name, set3_name, set4_name, set5_name),
                      fill = c("hotpink", "forestgreen", "yellow3", "darkslategray2", "gray"), 
                      alpha = rep(0.4, 5), 
                      lty = "blank", 
                      fontface = "bold",
                      fontfamily = "sans", 
                      cat.fontface = "bold",
                      cat.default.pos = "outer",
                      cat.fontfamily = "sans")
}

# ----------------------------------------------------------------------------------------------------------------------
# Function	: plot_five_venn
# Aim		: Overlap the CL replicatesion
# Input :
#   : tidy.ec: protein identification data.frame as generated by func:tidy_protein
# Output  	:  Venn plot
# ----------------------------------------------------------------------------------------------------------------------

plot_five_venn <- function(tidy_df){
  split_reps <- tidy_df %>% filter(treatment=="CL") %>%
    group_by(master_protein, rep)
  
  split_reps <- split(split_reps, split_reps$rep)
  
  grid.newpage()
  fullven <- five.venn(
    split_reps$`1`$master_protein, 
    split_reps$`2`$master_protein,
    split_reps$`3`$master_protein,
    split_reps$`4`$master_protein,
    split_reps$`5`$master_protein, 
    set1_name= "rep1",
    set2_name="rep2",
    set3_name="rep3", 
    set4_name="rep4",
    set5_name="rep5")
  
grid.draw(fullven) 

}
    
# ----------------------------------------------------------------------------------------------------------------------
# Function	: RBP
# Aim		:   
# Input :
# 	: foreground = vector of master_proteins assesed for enrichment
#   : background = vector of whole proteome with master_proteins
#   : id.type = annotation of master_proteins (uniprot, ensembl, tair, etc.)
#   : method: method to calculate enrichment ("BH" or "Hypergeometric")
#   : uniprot.go: go_mappings
# Output  	:  GO enrichment plots 
# NOTE: HAVE INCLUDED SHORTEN_TERM = FALSE IN getEnrichedGO function as original code based on Hs... 
# not sure about the implications of this decision
# ----------------------------------------------------------------------------------------------------------------------

compare.go <- function(foreground, background, id.type, method, uniprot.go){
  
  DE = (background$master_protein %in% foreground$master_protein)
  ABUN = log(background$abundance+1, 2)
  df = data.frame("ID"=background$master_protein,"DE"=DE,"ABUN"=ABUN)
  
  pwf <- makePWF(df, "DE", "ABUN", "ID")
  makePWFPlot(pwf, xlab = "protein_abundance")
  go_all <- getEnrichedGO(pwf, gene2cat=uniprot.go, shorten_term=FALSE,  method) 
  
  plotGOTerms(go_all, sum(as.numeric(DE)),
              length(df$ABUN), numObs_filter=10, 
              enrichment_filter=0, switch_axes=T, BH_filter=0.01)
}

# ----------------------------------------------------------------------------------------------------------------------
# Function	: myheatmap
# Aim		:   produce binary heatmap across the samples
# Input :
# 	: whole = dataframe of unique tidy output
#   : split_res = dataframe CL split by rep
#   : split_res.nc = dataframe NC split by rep
# Output  	: heatmap 10 channels
# ----------------------------------------------------------------------------------------------------------------------


myheatmap = function(whole, split_reps, split_reps.nc){

whole = whole
split_reps = split_reps
split_reps.nc = split_reps.nc
  
  
wholeR <- unique(whole$master_protein)
  
cl1 <- split_reps$`1`[, "master_protein"]
cl2 <- split_reps$`2`[, "master_protein"]
cl3 <- split_reps$`3`[, "master_protein"]
cl4 <- split_reps$`4`[, "master_protein"]
cl5 <- split_reps$`5`[, "master_protein"]

cl1 <- as.numeric(wholeR %in% cl1$master_protein)
cl2 <- as.numeric(wholeR %in% cl2$master_protein)
cl3 <- as.numeric(wholeR %in% cl3$master_protein)
cl4 <- as.numeric(wholeR %in% cl4$master_protein)
cl5 <- as.numeric(wholeR %in% cl5$master_protein)

nc1 <- split_reps.nc$`1`[, "master_protein"]
nc2 <- split_reps.nc$`2`[, "master_protein"]
nc3 <- split_reps.nc$`3`[, "master_protein"]
nc4 <- split_reps.nc$`4`[, "master_protein"]
nc5 <- split_reps.nc$`5`[, "master_protein"]

nc1 <- as.numeric(wholeR %in% nc1$master_protein)
nc2 <- as.numeric(wholeR %in% nc2$master_protein)
nc3 <- as.numeric(wholeR %in% nc3$master_protein)
nc4 <- as.numeric(wholeR %in% nc4$master_protein)
nc5 <- as.numeric(wholeR %in% nc5$master_protein)

full.df <- cbind(cl1, cl2, cl3, cl4, cl5, nc1, nc2, nc3, nc4, nc5)
breaks=seq(0, 1, by=0.5)
heatmap.2(full.df, breaks = breaks, tracecol = "black")
  
}

# wrapper for the above heatmap function but staight from "tidy" proteins identification level data
toms_my_heatmap <- function(tidy_df){

  split_reps_cl <- tidy_df %>% filter(treatment=="CL") %>%
    group_by(master_protein, rep)
                  
  split_reps_cl <- split(split_reps_cl, split_reps_cl$rep)
  
  split_reps_nc <- tidy_df %>% filter(treatment=="NC") %>%
    group_by(master_protein, rep)
  
  split_reps_nc <- split(split_reps_nc, split_reps_nc$rep)
  
  myheatmap(tidy_df, split_reps_cl, split_reps_nc)

  }

# ----------------------------------------------------------------------------------------------------------------------
# Function	:  myenrichment
# Aim		:   assess enrichment across the number of replicates
# Input :
#   : compiled.data = dataframe from compiler function
#   : background = all proteins identified in oops
#   : id.type = uniprot, tair, etc.
#   : sig_value = value for enrichment significance 
# Output  	: line plot of RBP enrichment 
# ----------------------------------------------------------------------------------------------------------------------

myenrichment = function(compiled.data, background, id.type, id, sig_value){

compiled.cl <- compiled.data[which(compiled.data$NC == 0),]

split.compiled <- split(compiled.cl, compiled.cl$CL)

one.cl    <- background$master_protein %in% split.compiled$`1`$master_protein
two.cl    <- background$master_protein %in% split.compiled$`2`$master_protein
three.cl  <- background$master_protein %in% split.compiled$`3`$master_protein
four.cl   <- background$master_protein %in% split.compiled$`4`$master_protein
five.cl   <- background$master_protein %in% split.compiled$`5`$master_protein

rbp.data = myProtMapper(unique(background$master_protein), id.type)
uniprot.go <- makeGene2Cat(rbp.data, "query", "go.all", ";")
rbpdata <- uniprot.go[which(uniprot.go$to.id == id), ]
RBP <- background$master_protein %in% rbpdata$query

full.df.2 <- cbind(background, RBP, one.cl, two.cl, three.cl, four.cl, five.cl)
#enrichment of cl proteins versus whole RBP, or need to consider only confirmed RBP in oops?

one.p <- print(runFishersEnrichment(full.df.2, "RBP", "one.cl"))
two.p <- print(runFishersEnrichment(full.df.2, "RBP", "two.cl"))
three.p <- print(runFishersEnrichment(full.df.2, "RBP", "three.cl"))
four.p <- print(runFishersEnrichment(full.df.2, "RBP", "four.cl"))
five.p <- print(runFishersEnrichment(full.df.2, "RBP", "five.cl"))

values <- c(one.p$p.value, two.p$p.value, three.p$p.value, four.p$p.value, five.p$p.value)
new.values <- log(values[1:5], 10)

r <- {plot(new.values, type = "l", ylab = "log{10}.p.values", col = "red")
  abline(h = log(sig_value), lty = 2)}

print(r)
cat(sprintf("\tNeed to identify master_protein at least %s replicates for a significant enrichment in RBPs\n", min(which(values < sig_value))))
}


# ----------------------------------------------------------------------------------------------------------------------
# Function	:  lineplots 
# Aim		:  line plot per number of times master_protein is seen 
# Input :  
#  compiled.data = output df from compiler function
# Output  	: line plot of proteins identified 
# ----------------------------------------------------------------------------------------------------------------------

lineplots = function(compiled.data){
  
cl.counts <- compiled.data[which(compiled.data$NC == 0),]
sep.cl <- split(cl.counts, cl.counts$CL)
sep <- split(compiled.data, compiled.data$CL) #nc != 0

lengths <- c(length(sep.cl$`1`$master_protein), 
             length(sep.cl$`2`$master_protein),
             length(sep.cl$`3`$master_protein),
             length(sep.cl$`4`$master_protein),
             length(sep.cl$`5`$master_protein))
lengths.with.nc <- c(length(sep$`1`$master_protein),
                     length(sep$`2`$master_protein),
                     length(sep$`3`$master_protein),
                     length(sep$`4`$master_protein),
                     length(sep$`5`$master_protein))
indexes <- c(1,2,3,4,5)
both <- cbind(indexes, lengths, lengths.with.nc)
colnames(both) <- c("indexes", "no nc", "with nc")
both <- as.data.frame(both)
{plot(both$indexes, both$`no nc`, type="o", col="black", pch="o", lty=1, ylim=c(0,700), 
      xlab = "number of times seen", ylab = "unique master proteins", main = "proteins identified")
  legend(3.5,700,legend=c("no NC proteins","with"), col=c("black","orange"), lty=c(1,1), ncol=1)
  points(both$indexes, both$`with nc`, col="orange", pch="o")
  lines(both$indexes, both$`with nc`, col="orange",lty=1)}
}

##
plotKEGG <- function(GO_df, len_foreground, len_background,
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
    xlab("") + ylab("Enrichment (Log2)") +
    scale_fill_continuous(name="BH adj. p-value\n(Log 10)\n") +
    my_theme +
    theme(
      text=element_text(size=15),
      plot.title=element_text(hjust=0.5)) +
    scale_x_discrete(labels=GO_filtered_df[,term.col],limits=GO_filtered_df[,term.col])+ coord_flip()
  
  return(list(p=p,g=GO_df))
}

#########TS: Doc. this!
runFTest <- function(not_in_foi, in_foi, len_back, len_foi_back){
  f_test <- fisher.test(matrix(c(in_foi, not_in_foi, len_foi_back, len_back-len_foi_back), nrow=2))
  return(c(f_test$p.value, f_test$estimate))
}
#########TS: Doc. this!


#########TS: Doc. this!
plotFOIfraction <- function(tidy_df, selected, foi, foi_name="RBP"){

  foi_identified <- tidy_df %>%
    filter(as.character(master_protein) %in% selected$uniprot_ac) %>%
    group_by(master_protein, treatment) %>%
    count() %>%
    spread(key=treatment, value=n) %>%
    replace_na(list(CL=0, NC=0)) %>%
    mutate(is_foi=master_protein %in% foi) %>%
    group_by(CL, NC, is_foi) %>%
    count() %>%
    spread(key=is_foi, value=n) %>%
    replace_na(list("FALSE"=0, "TRUE"=0)) %>%
    arrange(desc(CL), desc(NC))
  
  print(foi_identified)

  foi_identified$fraction <-  foi_identified[['TRUE']] / (foi_identified[['FALSE']] + foi_identified[['TRUE']])
  
  background <- unique(selected$uniprot_ac)
  
  foi_in_back <- intersect(background, foi)
  
  len_back <- length(background)
  len_foi_back <- length(foi_in_back)
  
  foi_identified[c("p", "log.odds")] <- t(apply(foi_identified, MARGIN=1, FUN=function(x) runFTest(x[['FALSE']], x[['TRUE']], len_back, len_foi_back)))
  foi_identified$label <- paste(foi_identified[['TRUE']], (foi_identified[['FALSE']]+foi_identified[['TRUE']]), sep=" / ")
  
  p <- ggplot(foi_identified, aes(CL, NC, fill=100*fraction)) +
    geom_tile() +
    scale_x_continuous(breaks=0:100) +
    scale_y_continuous(breaks=0:100) +
    scale_fill_gradient(name=sprintf("%s (%%) \n", foi_name), low="grey99", high="skyblue3") +
    geom_text(aes(label=label)) +
    my_theme
  
  p2 <- p + aes(fill=-(log10(p))) +
    scale_fill_gradient2(name=sprintf("%s over-rep\np-value\n(-log10)\n", foi_name), low="orangered2", high="skyblue3", mid="grey99", midpoint=1)
  
  p3 <- p + aes(fill=log.odds) +
    scale_fill_gradient2(name=sprintf("%s over-rep\nlog.odds\n", foi_name), low="orangered2", high="skyblue3", mid="grey99", midpoint=1)
  
  print(p)
  print(p2)
  print(p3)
  
  return(foi_identified)
}
#########TS: Doc. this!

myMakePWF <- function(df, sig_col, bias_col, identifier_col, cut_g=50){
  #print("making pwf")
  pwf <- df[,c(sig_col, bias_col)]
  #print(table(df[[sig_col]]))
  rownames(pwf) <- df[[identifier_col]]
  colnames(pwf) <- c("DEgenes", "bias.data")
  
  pwf$bias_col <- cut2(as.numeric(as.character(pwf$bias.data)), g=cut_g, levels.mean=T)
  if((sum(as.numeric(pwf$DEgenes))>40 & length(as.numeric(pwf$DEgenes))<0.8 | T==T)){
    spline <- makespline(pwf$bias_col, as.numeric(pwf$DEgenes))
  }
  else{
    spline=mean(as.numeric(df[[sig_col]]))
  }
  pwf$pwf <- spline
  
  return(pwf)
}

#########TS: Doc. this!
plotAbundanceAdjustedOverRep <- function(input_df, selected, go.full,
                                         RBP_GO_TERM="GO:0003723", min_obs=100, pwf_function=makePWF,
                                         NC_levels=None, CL_levels=None){
  
  foi_identified <- input_df %>%
    filter(master_protein %in% selected$uniprot_ac) %>%
    group_by(master_protein, treatment) %>%
    count() %>%
    spread(key=treatment, value=n) %>%
    replace_na(list(CL=0, NC=0))
  
  thresholds_selected <- selected
  thresholds_selected$log_abundance <- log2(thresholds_selected$abundance+1)
  go.rbp.only <- go.full %>% filter(GO.ID==RBP_GO_TERM)
  
  rows <- NULL
  plist <- NULL
  ix = 1
  
  if(missing(NC_levels)){
    NC_levels <- unique(foi_identified$NC)
  } 
  
  if(missing(CL_levels)){
    CL_levels <- unique(foi_identified$CL)
  } 

    all_identified <- foi_identified  %>% pull(master_protein) %>% unique()
  for(CL_n in CL_levels){
    for(NC_n in NC_levels){
      
      if (is.na(NC_n)){
        identified_proteins<- foi_identified %>% filter(CL==CL_n) %>% pull(master_protein) %>% unique()
      }
      
      else if(is.na(CL_n)){
        identified_proteins<- foi_identified %>% filter(NC==NC_n) %>% pull(master_protein) %>% unique()
      }
      
      else{
        identified_proteins<- foi_identified %>% filter(CL==CL_n, NC==NC_n) %>% pull(master_protein) %>% unique()
      }
      if(length(identified_proteins)>=min_obs){
        tmp_thresholds_selected <- thresholds_selected %>% filter(!uniprot_ac %in% setdiff(all_identified, identified_proteins))
        tmp_thresholds_selected$oops <- tmp_thresholds_selected$uniprot_ac %in% identified_proteins

        pwf <- pwf_function(tmp_thresholds_selected, sig_col="oops", bias_col="log_abundance", identifier_col="uniprot_ac")
        pwf_plot <- makePWFPlot(pwf, xlab=sprintf("log2 abundance %s:%s", CL_n, NC_n))
        plist[[ix]] <- pwf_plot + theme(text=element_text(size=5))
        goseq_results <- suppressMessages(goseq(pwf, gene2cat=go.full))
        print(head(goseq_results, 30))
        rows[[ix]] <- c(CL_n, NC_n,
                        sprintf("%s / %s", goseq_results %>% filter(category==RBP_GO_TERM) %>% pull(numDEInCat), length(identified_proteins)),
                        goseq_results %>% filter(category==RBP_GO_TERM) %>% pull(over_represented_pvalue))
        ix <- ix + 1
      }
    }
  }
  
  thresholds_RBP_over_rep <- data.frame(do.call("rbind", rows))
  colnames(thresholds_RBP_over_rep)=c("CL", "NC", "label", "p")
  thresholds_RBP_over_rep$p <- as.numeric(as.character(thresholds_RBP_over_rep$p))
  thresholds_RBP_over_rep$CL <- as.numeric(as.character(thresholds_RBP_over_rep$CL))
  thresholds_RBP_over_rep$NC <- as.numeric(as.character(thresholds_RBP_over_rep$NC))
  
  print(do.call("grid.arrange", c(plist, ncol=4)))
  print(thresholds_RBP_over_rep)
  p <- ggplot(data.frame(thresholds_RBP_over_rep), aes(CL, NC, fill=-log10(p))) +
    geom_tile() +
    scale_x_continuous(breaks=0:100) +
    scale_y_continuous(breaks=0:100) +
    geom_text(aes(label=label)) +
    scale_fill_gradient2(name="RBP over-rep\np-value\n(-log10)\n", low="orangered2", high="skyblue3", mid="grey99", midpoint=1, limits=c(0, NA)) +
    my_theme
  
  return(list("p"=p, "results"=thresholds_RBP_over_rep))


}
#########TS: Doc. this!
