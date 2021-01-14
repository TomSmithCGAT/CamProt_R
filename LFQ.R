library(tidyverse)

# Given a list of infiles, read the files, retain the seqeunce and protein ID columns and bind the rows
readInfiles <- function(infiles){
  data <- lapply(infiles, function(x){
    read.delim(x, stringsAsFactors=FALSE) %>% select('Sequence', 'Protein.Accessions', 'Master.Protein.Accessions') %>%
      mutate(infile=x)
  }) %>% do.call(what='rbind')
  
  return(data)
}

# Uses an approximate parsimony approach (e.g least number of proteins required to account for all observed sequences)
# to determine a single map from sequence to protein for all samples
getSequenceToProtein <- function(data, prot_col='Protein.Accessions'){
  
  seq_protein <- data %>% separate_rows(Protein.Accessions, sep='; ') %>% select(Sequence, prot_col) %>% distinct()
  protein_to_sequences <- seq_protein %>% group_by_at(prot_col) %>% summarise(Sequences=list(Sequence))
  protein_to_sequences_map <- protein_to_sequences$Sequences
  names(protein_to_sequences_map) <- protein_to_sequences[[prot_col]]
  
  protein_counts <- protein_to_sequences_map %>% names() %>% sapply(function(x) length(protein_to_sequences_map[[x]]))
  
  sequences_left_to_account_for <- unique(seq_protein$Sequence)
  sequences_accounted_for <- NULL
  n_seq_per_protein_values <- rev(names(table(protein_counts)))
  
  s2p <- NULL
  for(n in n_seq_per_protein_values){
    
    if(length(sequences_left_to_account_for)==0){
      break()
    }
    
    proteins <- protein_counts[protein_counts==n]
    
    mini_s2p <- rep(names(proteins), each=n)
    names(mini_s2p) <- unlist(protein_to_sequences_map[names(proteins)], use.names=FALSE)
    mini_s2p <- mini_s2p[!names(mini_s2p) %in% sequences_accounted_for]
    
    rep_seq <- intersect(names(mini_s2p), sequences_accounted_for)
    
    if(length(rep_seq>0)){
      stop(sprintf('Something has gone wrong here, accounting for the same sequence more than once %s',
                   paste(rep_seq, collapse=',')))
    }
    
    sequences_accounted_for <- union(sequences_accounted_for, names(mini_s2p))
    sequences_left_to_account_for <- setdiff(sequences_left_to_account_for, names(mini_s2p))
    
    s2p <- c(s2p, mini_s2p)
  }
  
  new_seq_to_master <- data.frame(Sequence=names(s2p), Protein=s2p) %>%
    group_by(Sequence) %>%
    summarise('Updated.Master.Protein.Accessions'=paste(Protein, collapse='; '))
}

# Check if two protein IDs are the same allowing for:
# 1. multiple IDs (default separated by '; ')
# 2. mutiple IDs in any order
matchProteins <- function(proteins1, proteins2, sep1='; ', sep2='; '){
  if(proteins1==proteins2){
    return(TRUE)
  }
  
  proteins1 <- strsplit(proteins1, sep1)[[1]]
  proteins2 <- strsplit(proteins2, sep2)[[1]]
  
  if(length(proteins1)!=length(proteins2)) return(FALSE)
  else if(length(setdiff(proteins1, proteins2))==0) return(TRUE)
  else return(FALSE)
}

# With object containing original ('Master.Protein.Accessions') and updated ('Updated.Master.Protein.Accessions')
# master protein assingments and Sequences, compare the number of sequences per protein
compareSequencesPerProtein <- function(data_plus_new_master){
  
  master_protein_counts <- data_plus_new_master %>%
    filter(Master.Protein.Accessions!='') %>%
    select(Sequence, Master.Protein.Accessions) %>% distinct() %>% 
    group_by(Master.Protein.Accessions) %>%
    tally() %>%
    group_by(n) %>%
    tally() %>%
    mutate(type='Original')
  
  new_master_protein_counts <- data_plus_new_master %>%
    filter(Updated.Master.Protein.Accessions!='') %>%
    select(Sequence, Updated.Master.Protein.Accessions) %>% distinct() %>% 
    group_by(Updated.Master.Protein.Accessions) %>%
    tally() %>%
    group_by(n) %>%
    tally() %>%
    mutate(type='Updated')
  
  p <- new_master_protein_counts %>% rbind(master_protein_counts) %>%
    arrange(type) %>%
    mutate(type=factor(type, levels=c('Original', 'Updated'))) %>%
    group_by(type) %>%
    mutate(cum_nn =cumsum(nn)/sum(nn)) %>%
    ggplot(aes(log10(n), cum_nn, colour=type)) +
    geom_step() + 
    theme_bw() +
    theme(aspect.ratio=1) +
    scale_colour_discrete(name='Master protein assignment') +
    ylim(0, NA) +
    xlab('Sequences per protein (log10)') +
    ylab('Fraction of proteins')
  
  invisible(p)
}

# With object containing original ('Master.Protein.Accessions') and updated ('Updated.Master.Protein.Accessions')
# master protein assingments, identify the proportion of exact matches in protein IDs
compareIDs <- function(data_plus_new_master){
  id_match <- data_plus_new_master %>%
    rowwise() %>%
    mutate(same_id=ifelse(matchProteins(Master.Protein.Accessions, Updated.Master.Protein.Accessions), 'Same ID(s)', 'Different ID(s)')) %>%
    pull(same_id) %>%
    table()
  
  print(id_match)
  print(round(100*id_match/sum(id_match), 1))
}

# With object containing original ('Master.Protein.Accessions') and updated ('Updated.Master.Protein.Accessions')
# master protein assingments, identify the proportion with a single assigned master protein
compareSingleMaster <- function(data_plus_new_master){
  
  single_master <- data_plus_new_master %>% 
    mutate(original=ifelse(grepl('; ', Master.Protein.Accessions), 'Original:Multiple', 'Original:Single'),
           updated=ifelse(grepl('; ', Updated.Master.Protein.Accessions), 'Updated:Multiple', 'Updated:Single')) %>%
    group_by(original, updated) %>%
    count() %>%
    spread(key=updated, value=n) %>%
    tibble::column_to_rownames('original')
  
  print(single_master)
  print(round(100*single_master/sum(single_master), 2))
}

# Wrapper function to read in infiles, identify the parsimous proteins to account for the peptides and
# then compare the master protein assignments
getParsimonyPepToProt <- function(infiles, prot_col='Protein.Accessions', verbose=TRUE){
  data <- readInfiles(infiles)
  
  new_seq_to_master <- getSequenceToProtein(data, prot_col)
  
  data_plus_new_master <- data %>% merge(new_seq_to_master, by='Sequence')
  
  if(verbose){

    cat(sprintf('With original assignments, %s master proteins. With update, %s master proteins\n',
               length(unique(data_plus_new_master$Master.Protein.Accessions)),
               length(unique(data_plus_new_master$Updated.Master.Protein.Accessions))))
    
    print(compareSequencesPerProtein(data_plus_new_master))
    cat('Comparing Master Protein IDs\n')
    compareIDs(data_plus_new_master)
    compareSingleMaster(data_plus_new_master)
  }
  
  invisible(new_seq_to_master)
}


#From an MsnSet with peptide level 
removeLowSampleFeatures <- function(obj, min_samples, plot=TRUE){
  
  if(plot){
    p <- data.frame(table(rowSums(is.finite(exprs(obj))))) %>%
      mutate(Var1=factor(Var1, levels=rev(levels(Var1)))) %>%
      ggplot(aes(x=1, y=Freq, fill=Var1)) +
      geom_bar(stat='identity') +
      theme_bw() +
      scale_fill_discrete(name='# samples') +
      ylab('Frequency') +
      theme(aspect.ratio=4, axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
    
    print(p)
  }
  
  invisible(filterNA(obj, pNA=((ncol(obj)-min_samples)/ncol(obj))))
}

# From a tidied MSnSet, with a column of 'ratio' values,
# obtain the correction factor the ratios
getCorrectionFactor <- function(.data, group_col='Replicate', center=FALSE,
                                ratio_col='ratio', median_ratio_col='median_ratio'){
  
  cor_factors <- .data %>%
    group_by_at(vars(one_of(group_col))) %>%
    summarise(!!(median_ratio_col):=median(!!sym(ratio_col), na.rm=TRUE))
  
  if(center){ # center median ratios
    cor_factors$median_ratio <- cor_factors$median_ratio-mean(cor_factors$median_ratio)
  }
  
  invisible(cor_factors)
}


# From a tidied MSnSet, with a column of 'ratio' values
# and a vector of correction factors (see getCorrectionFactor()),
# correct the ratios
correctRatios <- function(.data, .cor_factor, group_col=Replicate,
                          ratio_col='ratio', median_ratio_col='median_ratio',
                          corrected_ratio_col='corrected_ratio'){
  return_data <- .data %>% merge(.cor_factor, by=group_col)
  return_data[[corrected_ratio_col]] <- return_data[[ratio_col]] - return_data[[median_ratio_col]] 
  invisible(return_data)
}
