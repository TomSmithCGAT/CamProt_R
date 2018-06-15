getRatio <- function(treatment, control, treatment_name, control_name){
  
  if (is.na(control)|is.na(treatment)){
    
    if (is.na(control) & is.na(treatment)){
      class = "Both Missing"
    }
    else if (is.na(control)){
      class = sprintf("%s Missing", control_name)
    }
    else {
      class = sprintf("%s Missing", treatment_name)
    }
    Ratio = NA
  }
  else{
    Ratio = treatment - control
    class = "Neither Missing"
  }
  return(c(Ratio, class))
}

my_aggregator <- function(values, func=mean){
  if(sum(is.na(values)==length(values))){
    return(NA)
  }
  else{
    return(func(values[is.finite(values)]))
  }
}


getMissingRatio <- function(treatment, control, treatment_name, control_name){
  
  if (is.na(control)|is.na(treatment)){
    Ratio = NA
    if (is.na(control) & is.na(treatment)){
      class = "Both Missing"
      pseudo_Ratio = NA
    }
    else if (is.na(control)){
      class = sprintf("%s Missing", control_name)
      pseudo_Ratio = log10((2^treatment) + 1)
    }
    else {
      class = sprintf("%s Missing", treatment_name)
      pseudo_Ratio = log10(1/((2^control) + 1))
    }
    
  }
  else{
    Ratio = log10(((2^treatment)+1) / ((2^control) + 1))
    pseudo_Ratio = log10(((2^treatment)+1) / ((2^control) + 1))
    class = "Neither Missing"
  }
  return(c(Ratio, pseudo_Ratio, class))
}
