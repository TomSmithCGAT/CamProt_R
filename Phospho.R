suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

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


### combine_peptide_phospho_positions ###

addPhosphoPositions <- function(obj){
  
  combine_peptide_phospho_positions <- function(peptide_start, filtered_pos){
    # Given a string with start positions (peptide_start) and positions of phosphorylated AA with
    # respect to peptide sequence (phospho_positions), return the phospho positions with respect to protein sequence
    
    if(filtered_pos==""){
      return(NA)
    }
    return_string <- NULL
    #print(peptide_start)
    for(p_start in strsplit(peptide_start, split=";")[[1]]){
      position_string <- NULL    
      
      for(phospho_p in strsplit(filtered_pos, split=";")[[1]]){
        position_string[[phospho_p]] <- as.numeric(p_start) + as.numeric(phospho_p)
        
      }
      return_string[[p_start]] <- paste0(position_string, collapse="|")
    }
    return(paste0(return_string, collapse=";"))
  }
  
  phospho_positions <- strsplit(obj$Positions.in.Master.Proteins, split="; ")
  
  start_array <- rep("", length(phospho_positions))
  protein_array <- rep("", length(phospho_positions))
  
  for(ix in seq_along(phospho_positions)){
    phospho_position <- phospho_positions[[ix]]
    start_values <- NULL
    protein_values <- NULL
    protein <- NA
    for(site in phospho_position){
      site_split <- strsplit(site, " ")[[1]]
      if (length(site_split)==1){
        start_values[[site]] <- gsub("\\[|\\]", "", strsplit(site_split[[1]], "-")[[1]][1])
      }
      
      else if (length(site_split)==2){
        protein <- site_split[[1]]
        start_values[[site]] <- gsub("\\[|\\]", "", strsplit(site_split[[2]], "-")[[1]][1])
      }
      
      else{
        stop(sprintf("Unexpected input. Position should contain a maximum of two space-separated elements: %s", site))
      }
      protein_values[[site]] <- protein
    }
    start_array[[ix]] <- paste0(start_values, collapse=";")
    protein_array[[ix]] <- paste0(protein_values, collapse=";")
  }
  
  obj$peptide_start <- start_array
  obj$phospho_position_protein  <- protein_array
  
  obj$phospho_position <- apply(
    obj, MARGIN=1, function(x) combine_peptide_phospho_positions(x[["peptide_start"]], x[["filtered_pos"]]))
  
  return(obj)
}
