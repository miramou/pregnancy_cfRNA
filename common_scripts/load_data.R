#Load and process data
load_process_qPCR = function(data, common_set, is_classification, filter_5695 = TRUE) {
  require(plyr)
  require(dplyr)

  if (filter_5695) {
    data = data %>%  filter(panel != '5695')
  }
  
  if(is_classification) {
    data = data[,colnames(data) %in% c(common_set, 'ga', 'dCD', 'delivery', 'patient')]
    data$Class = ifelse(data$delivery<37, 'PreTerm', 'Term')
    data$Class = factor(data$Class)
    data = data %>% filter(dCD <=0)
  } else {
    data = data[,colnames(data) %in% c(common_set, 'ga', 'dCD', 'delivery', 'patient', 'sample_id', 'cohort')]
  }
  return(data)
}
