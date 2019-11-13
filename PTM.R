suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
require("GenomicRanges")
require("Biostrings")

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

# TS: Note, this function was initially copied entirely from Phospho.R and
# hard-coded references to phosphorylation removed
parsePTMScores <- function(obj, threshold=95, ptm_col="PhosphoRS.Best.Site.Probabilities",
                           prob_split='; |: ', collapse_delimiter=";", verbose=TRUE){
  
  if(class(obj)!="data.frame"){
    stop("first argument must be a data.frame")
  }
  
  cat(sprintf("Removed %s peptides where the ptm_col value == `Inconclusive data`\n",
              sum(obj[[ptm_col]]=="Inconclusive data")))
  obj <- obj[obj[[ptm_col]]!="Inconclusive data",]
  
  
  # initiate vectors with empty string
  # where the PTM(s) pass the threshold, these will be updated
  filtered_ptm_desc <- filtered_ptm_res <- filtered_ptm_pos <- filtered_ptm <- filtered_ptm_score <- rep("", nrow(obj))
  
  split_probabities <- strsplit(obj[[ptm_col]], split=prob_split)
  
  log <- list("Total peptides"=0, "Total detected PTMpeptides"=0,
              "Peptides passing filter"=0, "Peptides failing filter"=0,
              "BiPTM/multiPTM peptides where some sites fail filter"=0,
              "Total detected sites"=0, "Sites passing filter"=0, "Sites failing filter"=0,
              "monoPTM passing filter"=0, "biPTM passing filter"=0, "multiPTM passing filter"=0,
              "Too many isoforms"=0)
  
  
  for (i in seq_along(split_probabities)){
    
    log[["Total peptides"]] <- log[["Total peptides"]] +1
    
    peptide_ptm_scores <- split_probabities[[i]]
    if(length(peptide_ptm_scores)==0){next()} # no PTM detected
    if(is.na(peptide_ptm_scores[[1]])){ next() } # no PTM detected
    
    log[["Total detected PTMpeptides"]] <- log[["Total detected PTMpeptides"]] + 1
    
    if (peptide_ptm_scores[[1]] == "Too many isoforms"){
      log[["Too many isoforms"]] <- log[["Too many isoforms"]] + 1
      log[["Peptides failing filter"]] <- log[["Peptides failing filter"]] + 1
      next()
    }
    
    log[["Total detected sites"]] <- log[["Total detected sites"]] + length(peptide_ptm_scores)/2
    
    scores <- peptide_ptm_scores[seq(2, length(peptide_ptm_scores), 2)]
    
    # if any score is below threshold, disregard all putative ptm sites
    if (any(as.numeric(scores)<threshold)){
      log[["Sites failing filter"]] <- log[["Sites failing filter"]] + length(peptide_ptm_scores)/2
      log[["Peptides failing filter"]] <- log[["Peptides failing filter"]] + 1
      if (any(as.numeric(scores)>=threshold)){
        log[["BiPTM/multiPTM peptides where some sites fail filter"]] <- log[[
          "BiPTM/multiPTM peptides where some sites fail filter"]] + 1  
      }
      # if we want to handle this differently, can implement an alternative approach here
      # and move the rest of the code below into an else clause
      next()
    }
    
    log[["Sites passing filter"]] <- log[["Sites passing filter"]] + length(peptide_ptm_scores)/2
    log[["Peptides passing filter"]] <- log[["Peptides passing filter"]] + 1
    
    if(length(scores)==1){
      log[["monoPTM passing filter"]] <- log[["monoPTM passing filter"]] + 1
    }
    else if(length(scores)==2){
      log[["biPTM passing filter"]] <- log[["biPTM passing filter"]] + 1
    }
    else{
      log[["multiPTM passing filter"]] <- log[["multiPTM passing filter"]] + 1
    }
    
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
      cat(sprintf("%s: %i\n", event, log[[event]]))
    }
  }
  
  return(obj)
}


# TS: Note, this function was initially copied entirely from Phospho.R and
# hard-coded references to phosphorylation removed. It has since been modified 
# to use a proteome fasta to identify the position of the PTM with respect to the protein
addPTMPositions <- function(obj, proteome_fasta, master_protein_col="Master.Protein.Accessions"){
  
  proteome <- readAAStringSet(proteome_fasta)
  names(proteome) <- sapply(base::strsplit(names(proteome), split="\\|") , "[[", 2)
  
  combine_peptide_ptm_positions <- function(proteome, protein, sequence, filtered_pos){
    # Given a master protein(s), sequence and positions of PTM AA with respect to peptide sequence (filtered_pos),
    # return the PTM positions with respect to protein sequence
    
    if(filtered_pos==""){
      return(NA)
    }
    return_string <- NULL
    
    if(length(proteins)>1){
      return("")
    }
    
    peptide_starts <- start(matchPattern(sequence, proteome[[protein]]))

    for(p_start in peptide_starts){
      position_string <- NULL    
      for(ptm_p in strsplit(filtered_pos, split=";")[[1]]){
        position_string[[ptm_p]] <- p_start + as.numeric(ptm_p) - 1
        
      }
      return_string[[as.character(p_start)]] <- paste0(position_string, collapse="|")
    }
    
    return(paste0(return_string, collapse=";"))
  }
  
  obj$ptm_position <- apply(
    obj, MARGIN=1, function(x) combine_peptide_ptm_positions(
      proteome, x[[master_protein_col]], x[["Sequence"]], x[["filtered_pos"]]))
  
  return(obj)
}

### Get the Sequence around the PTMsite. Will return NA if peptide maps to multiple proteins or has multiple PTMs ###
# TS: Note, this function was initially copied entirely from Phospho.R and
# hard-coded references to phosphorylation removed
getSequence <-function(proteome, Master.Protein.Accessions, ptm_position, pad=7){
  
  ptm_position <- suppressWarnings(as.numeric(as.character(ptm_position)))
  
  if(is.na(ptm_position)){
    return(NA)
  }
  
  if(grepl("; ", ptm_position)){
    return(NA)
  }

  if(!Master.Protein.Accessions %in% names(proteome)){
    return(NA)
  }
  
  protein_length <- length(proteome[[Master.Protein.Accessions]])
  
  if(ptm_position>protein_length){
    warning(sprintf("PTM positions is outside protein sequence! Returning NA. %s: [-%s], PTM: %s",
                    Master.Protein.Accessions, protein_length, ptm_position))
    return(NA)
  }
  
  start_pad <- end_pad <- ""
  
  start <- ptm_position - pad
  if(start < 0){
    start_pad <- paste0(rep("_", start*-1), collapse="")
    start <- 0
  }
  
  
  
  end <- ptm_position + pad
  if (end > protein_length){
    end_pad <- paste0(rep("_", (end - protein_length)), collapse="")
    end <- protein_length
  }
  
  mod_position <- pad + 1
  
  sequence <- as.character(proteome[[Master.Protein.Accessions]][start:end])
  sequence <- paste0(start_pad, sequence, end_pad)
  
  sequence <- paste(base::substr(sequence, 1, pad),
                    tolower(base::substr(sequence, pad+1, pad+1)),
                    base::substr(sequence, pad+2, pad+pad+1), sep="")
  
  return(sequence)
  
}

# TS: Note, this function was initially copied entirely from Phospho.R and
# hard-coded references to phosphorylation removed
addSiteSequence <- function(obj, proteome_fasta){
  require("GenomicRanges")
  require("Biostrings")
  
  proteome <- readAAStringSet(proteome_fasta)
  names(proteome) <- sapply(base::strsplit(names(proteome), split="\\|") , "[[", 2)
  
  obj <- obj %>% rowwise() %>%
    mutate(site_seq = getSequence(proteome, Master.Protein.Accessions, ptm_position)) %>%
    data.frame()
  
  return(obj)
}




# ----------------------------------------------------------------------------------------------------------------------
# Function	: parsePSMAndAggregate
# Aim		: parse PSMs and prepare for LOPIT analysis. 
#
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
#   plots_dir: [OPTIONAL] Where to save plots. If empty, no plots are saved
#   modifications: [OPTIONAL] Expected modifications to filter against, e.g c("Succinyl", "Malonyl", "Glutaryl"). 
#   SN_threshold: [REQUIRED; DEFAULT=5] Signal:noise threshold to retain PSM quantification values
#   subset_cols_ix: [OPTIONAL] Subset the TMT dataset to these columns, e.g 1:10
#   intensity_filter: [REQUIRED; DEFAULT=2^2.25], # Minimum intensity value
#   interference_threshold: [REQUIRED; DEFAULT=50], # Maximum allowed interference (e.g co-isolation)
#   total_data: [OPTIONAL] Center-normalise using this "total" data rather than the observed intensities
#   max_missing: [REQUIRED; DEFAULT=3] Maximum number of allowed missing values to be imputed (using KNN, k=10)
#   mod_col: [REQUIRED; DEFAULT="Modifications"] column to identify the modifications for grouping PSM into pep seq + mod. 
#            Note that if you have previous parsed the input with `parsePTMScores` this will add a `ptm_position` and
#            `filtered_ptm` column which is a more suitable value here since this will ensure all peptides with the same 
#            sequence and PTMs are aggregated
#   agg_to_unique_pep: [OPTIONAL] set to TRUE (for total, e.g unmodified) peptides to aggregate to peptides, e.g ignoring their
#                      modifications
#   agg_to_prot: [OPTIONAL] set to TRUE (for total, e.g unmodified) peptides to aggregate to protein
#
# Output: A list with three values:
#         datasets - The aggregated datasets. All are MsnSet objects
#         descriptions - Descriptions of the aggregated datasets
#         plots - The plots generated during parsing
# ----------------------------------------------------------------------------------------------------------------------


parsePSMAndAggregate <- function(
  raw_psm, # the PSM level dataframe
  sample_infile, # file mapping TMT tags to sample names
  plots_dir, # Where to save plots
  master_protein_col="Master.Protein.Accessions", # master protein column name
  modifications, # Expected modifications to filter against, e.g c("Succinyl", "Malonyl", "Glutaryl")
  SN_threshold=5, # Signal:noise threshold to retain PSM quantification values
  subset_cols_ix, # Subset the TMT dataset to these columns, e.g 1:10
  intensity_filter=2^2.25, # Minimum intensity value, e.g 2^2.25
  interference_threshold=50, # Maximum interference
  total_data, # Center-normalise using this "total" data rather than the observed intensities
  max_missing=3, # Maximum number of allowed missing values
  sequential_missing_max=6, # minimum (total) number of missing values to allow for sequential zero imputation
  sequential_missing_min=4, # minimum number of sequential missing values for zero imputation
  mod_col="Modifications", # column to identify the modifications for grouping PSM into pep seq + mod
  agg_to_unique_pep=FALSE, # set to TRUE (for total, e.g unmodified) peptides to agg. to pep. ignoring mods
  agg_to_prot=FALSE, # set to TRUE (for total, e.g unmodified) peptides to agg. to protein
  MS_level="MS3"
){
  
  printModTally <- function(data, modifications){
    for(mod in modifications){
      data[[mod]] <- as.numeric(grepl(mod, data$Modifications))
    }
    
    print(data %>% group_by(.dots=lapply(modifications, as.symbol)) %>% tally())
  }
  
  plots <- list()
  datasets <- list()
  descriptions <- list()
  
  if(hasArg(modifications)){
    printModTally(raw_psm, modifications)
  }
  
  ###################################################
  # filter out peptides without desired modifications
  ###################################################
  if(hasArg(modifications)){
    cat("Removing PSMs without required modifications\n")
    modifications_grep <- paste(modifications, collapse="|")
    raw_psm_only_mod <- raw_psm %>% filter(grepl(modifications_grep, Modifications))
    
    cat(sprintf("\nAll PSMs: %s, PSMs with modifications: %s\n", nrow(raw_psm), nrow(raw_psm_only_mod)))
  }
  else{
    raw_psm_only_mod <- raw_psm
  }
  
  # Make MsnSet
  #if(MS_level=="MS2"){
  #  raw_res <- makeMSNSet(raw_psm_only_mod, sample_infile, level="PSM", ab_col_ix=2) 
  #}
  #else if(MS_level=="MS3"){
  #  raw_res <- makeMSNSet(raw_psm_only_mod, sample_infile, level="PSM", ab_col_ix=3) 
  #}
  raw_res <- makeMSNSet(raw_psm_only_mod, sample_infile, level="PSM", ab_col_ix=2) 
  
  if(hasArg(modifications)){
    printModTally(fData(raw_res), modifications)
  }
  
  ####################################################
  # Remove PSMs with Interference above threshold
  ####################################################
  cat("Removing PSMs with high interference (e.g co-isolation)\n")
  raw_res_int <- raw_res[fData(raw_res)$Isolation.Interference.in.Percent<=interference_threshold,]
  cat(sprintf("PSMs with co-isolation interference under threshold(%s): %s/%s\n",
              round(interference_threshold, 2), length(rownames(raw_res_int)), length(rownames(raw_res))))
  
  cat(sprintf("Proteins retained = %s/%s\n",
              length(unique(fData(raw_res_int)[[master_protein_col]])),
              length(unique(fData(raw_res)[[master_protein_col]]))))
  
  ####################################################
  # Remove PSMs with Signal:Noise (SN) below threshold
  ####################################################
  if(MS_level=="MS3"){
    cat("Removing PSMs with low Signal:noise ratio\n")
    raw_res_sn <- raw_res_int[fData(raw_res_int)$Average.Reporter.SN>=SN_threshold,]
    cat(sprintf("PSMs with SN over threshold(%s): %s/%s\n",
                round(SN_threshold, 2), length(rownames(raw_res_sn)), length(rownames(raw_res_int))))
    
    cat(sprintf("Proteins retained = %s/%s\n",
                length(unique(fData(raw_res_sn)[[master_protein_col]])),
                length(unique(fData(raw_res_int)[[master_protein_col]]))))
  }
  else{
    raw_res_sn <- raw_res_int
  }
  
  ###################################################
  # plot distribution of intensities per tag
  ###################################################
  p <- plotLabelQuant(raw_res_sn, log=TRUE, print=FALSE)$p2 +
    geom_vline(xintercept=2.25)
  
  plots[["pep_itensities_sn"]] <- p
  
  if(hasArg(plots_dir)){
    ggsave(file.path(plots_dir, "pep_intensities_sn.png"))
  }
  
  
  ###################################################
  # plot missing values at PSM level
  ###################################################
  if(hasArg(plots_dir)){
    png(file.path(plots_dir, "missing_peptide.png"))
    plotMissing(raw_res_sn, verbose=FALSE)
    dev.off()
  }
  
  plotMissing(raw_res_sn)
  
  
  ###################################################
  # Remove low intensity estimates
  ###################################################
  cat("Replacing low intensity values with NA\n")
  # If we want to subset the tags, e.g one contains total lysate, subset here
  if(hasArg(subset_cols_ix)){
    raw_res_sn_filter_low <- raw_res_sn[,subset_cols_ix]
  }
  else{
    raw_res_sn_filter_low <- raw_res_sn
  }
  
  exprs(raw_res_sn_filter_low)[exprs(raw_res_sn_filter_low)<intensity_filter] <- NA
  
  plotMissing(raw_res_sn_filter_low)
  
  if(hasArg(plots_dir)){
    png(file.path(plots_dir, "missing_peptide_low_intensity_filtered.png"))
    plotMissing(raw_res_sn_filter_low, verbose=FALSE)
    dev.off()
  }
  
  p <- plotLabelQuant(raw_res_sn_filter_low, log=TRUE, print=FALSE)
  plots[["pep_itensities_sn_intensity_filter_box"]] <- p$p1
  plots[["pep_itensities_sn_intensity_filter_density"]] <- p$p2
  
  datasets[["sn_intensity_filter"]] <- raw_res_sn_filter_low
  descriptions[["sn_intensity_filter"]] <- "PSM level, Signal:noise filtered, low intensity filtered"
  
  ###################################################
  # Normalise 
  ###################################################
  cat("Center-median normalising data\n")
  if(hasArg(total_data)){
    cat("Using medians from total data to center normalise\n")
    raw_res_sn_filter_low_cm_norm <- centerNormalise(raw_res_sn_filter_low, getMedians(total_data))
  }
  else{
    raw_res_sn_filter_low_cm_norm <- centerNormalise(raw_res_sn_filter_low)
  }
  
  p <- plotLabelQuant(raw_res_sn_filter_low_cm_norm, log=TRUE, print=FALSE)
  plots[["pep_itensities_sn_intensity_filter_sm_normalised_box"]] <- p$p1
  plots[["pep_itensities_sn_intensity_filter_sm_normalised_density"]] <- p$p2
  
  datasets[["sn_intensity_filter_sm_normalised"]] <- raw_res_sn_filter_low_cm_norm
  descriptions[["sn_intensity_filter_sm_normalised"]] <- "
  PSM level, Signal:noise filtered, low intensity filtered, center-median normalised"
  
  ######################################################
  # Aggregate to unique peptide sequence + modifications
  ######################################################
  cat("Aggregating to unique peptide sequence + modifications\n")
  
  agg_pep_mod <- agg_to_peptide_mod(raw_res_sn_filter_low_cm_norm, mod_col=mod_col)
  
  cat(sprintf("%s unique seq + mod peptides\n", length(rownames(agg_pep_mod))))
  
  plotMissing(agg_pep_mod)
  
  if(hasArg(plots_dir)){
    png(file.path(plots_dir, "missing_peptide_agg_seq_mod.png"))
    plotMissing(agg_pep_mod, verbose=FALSE)
    dev.off()
  }
  
  p <- plotLabelQuant(agg_pep_mod, log=TRUE, print=FALSE)
  plots[["pep_agg_seq_mod_intensities_box"]] <- p$p1
  plots[["pep_agg_seq_mod_intensities_density"]] <- p$p2
  
  datasets[["agg_pep_mod"]] <- agg_pep_mod
  descriptions[["agg_pep_mod"]] <- "Aggregated unique peptide sequence + modifications level"
  
  ######################################################
  # Sum normalise
  ######################################################
  cat("Sum normalisation\n")
  agg_pep_mod_sum_norm <- normalise(agg_pep_mod, "sum")
  
  # make a unique rowname. For now this is just the sequence and modification concatenated
  rownames(agg_pep_mod_sum_norm) <- paste0(fData(agg_pep_mod_sum_norm)$Sequence, "_",
                                           fData(agg_pep_mod_sum_norm)[[mod_col]])
  
  ######################################################
  # Impute missing values
  ######################################################
  
  cat(sprintf("Imputing between %s-%s sequential missing values with zero\n",
              sequential_missing_min, sequential_missing_max))
  agg_pep_mod_sum_norm_impute <- replace_missing_not_at_random(
    agg_pep_mod_sum_norm,
    max_total_missing=sequential_missing_max,
    min_sequential_missing=sequential_missing_min)
  
  cat(sprintf("Imputing up to %s missing values by KNN\n", max_missing))
  agg_pep_mod_sum_norm_impute <- ImputeMissing(
    agg_pep_mod_sum_norm_impute, missing=max_missing)
  
  datasets["agg_pep_mod_sum_norm_impute"] <- agg_pep_mod_sum_norm_impute
  descriptions[["agg_pep_mod_sum_norm_impute"]] <- "
  Aggregated unique peptide sequence + modifications level, sum normalised, imputed"
  
  cat(sprintf("At the end, we have retained %s peptides in %s proteins\n",
              length(rownames(agg_pep_mod_sum_norm_impute)),
              length(unique(fData(agg_pep_mod_sum_norm_impute)[[master_protein_col]]))))
  
  ######################################
  # Aggregate to unique peptide sequence, sum normalise and impute
  # This is for total peptides where we want aggregation to unique sequences,
  # e.g ignoring variable modifications 
  ######################################
  if(agg_to_unique_pep){
    cat("Aggregating to unique peptide sequence, then sum normalise, then impute\n")
    agg_pep <- agg_to_peptides(raw_res_sn_filter_low_cm_norm)
    
    cat(sprintf("%s unique seq\n", length(rownames(agg_pep))))
    
    plotMissing(agg_pep)
    
    if(hasArg(plots_dir)){
      png(file.path(plots_dir, "missing_peptide_agg_seq.png"))
      plotMissing(agg_pep, verbose=FALSE)
      dev.off()
    }
    
    p <- plotLabelQuant(agg_pep, log=TRUE, print=FALSE)
    plots[["pep_agg_seq_intensities_box"]] <- p$p1
    plots[["pep_agg_seq_intensities_density"]] <- p$p2
    
    datasets[["agg_pep"]] <- agg_pep
    descriptions[["agg_pep"]] <- "Aggregated unique peptide sequence level"
    
    cat("Sum normalisation on unique pep sequence level\n")
    agg_pep_sum_norm <- normalise(agg_pep, "sum")
    
    cat(sprintf("Imputing up to %s missing values by KNN at unique pep sequence level\n", max_missing))
    agg_pep_sum_norm_impute <- ImputeMissing(
      agg_pep_sum_norm, missing=max_missing)
    
    datasets[["agg_pep_sum_norm_impute"]] <- agg_pep_sum_norm_impute
    descriptions[["agg_pep_mod_sum_norm_impute"]] <- "Aggregated unique peptide sequence level, sum normalised, imputed"
    
    cat(sprintf("At the end, we have retained %s peptides in %s proteins\n",
                length(rownames(agg_pep_sum_norm_impute)),
                length(unique(fData(agg_pep_sum_norm_impute)[[master_protein_col]]))))
    
  }
  
  
  ######################################
  # Aggregate to unique protein, sum normalise and impute
  # This is useful for total, unmodified proteins only,
  ######################################
  if(agg_to_prot){
    
    if(!agg_to_unique_pep){
      stop("Need to also set `agg_to_unique_prot` to TRUE to aggregate to protein level")
    }
    
    cat("Aggregating to unique protein, then sum normalise, then impute\n")
    
    agg_prot <- agg_to_protein(agg_pep, master_protein_col)
    
    cat(sprintf("%s unique proteins\n", length(rownames(agg_prot))))
    
    plotMissing(agg_prot)

    if(hasArg(plots_dir)){
      png(file.path(plots_dir, "missing_protein.png"))
      plotMissing(agg_prot, verbose=FALSE)
      dev.off()
    }
    
    p <- plotLabelQuant(agg_prot, log=TRUE, print=FALSE)
    plots[["prot_agg_intensities_box"]] <- p$p1
    plots[["prot_agg_intensities_density"]] <- p$p2
    
    
    datasets[["agg_prot"]] <- agg_prot
    descriptions[["agg_prot"]] <- "Aggregated unique protein ID level"
    
    cat("Sum normalisation on unique protein ID level\n")
    agg_prot_sum_norm <- normalise(agg_prot, "sum")
    
    cat(sprintf("Imputing up to %s missing values by KNN at unique protein ID level\n", max_missing))
    agg_prot_sum_norm_impute <- ImputeMissing(
      agg_prot_sum_norm, missing=max_missing)
    
    datasets[["agg_prot_sum_norm_impute"]] <- agg_prot_sum_norm_impute
    descriptions[["agg_prot_sum_norm_impute"]] <- "Aggregated unique protein ID level, sum normalised, imputed"
    
    cat(sprintf("At the end, we have retained %s proteins\n", length(rownames(agg_prot_sum_norm_impute))))
    
  }
  
  if(hasArg(plots_dir)){
    for(p_name in names(plots)){
      ggsave(filename=file.path(plots_dir, sprintf("%s.png", p_name)), plot=plots[[p_name]])
    }
  }
  
  if(hasArg(modifications)){
    printModTally(fData(agg_pep_mod_sum_norm_impute), modifications)
  }
  
  
  invisible(list("datasets"=datasets, "descriptions"=descriptions, "plots"=plots))
}
  