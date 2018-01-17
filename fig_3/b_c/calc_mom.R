filter_bin_data = function(counts_all) {
  #Remove any post-partum samples
  counts_all = counts_all %>% filter(dCD<=0)
  
  #Bin for multiple of the median calculation
  counts_all$bin = counts_all$cohort
  counts_all$bin[counts_all$bin=='PT_UAB' | counts_all$bin =='T_UAB' ] = 'T_UAB'
  counts_all$bin[counts_all$bin == 'PT_Penn' | counts_all$bin == 'T_Penn'] = 'T_Penn'
  return(counts_all)
}

#Find median per trimester
calc_median_trimester = function(counts_all, trimester){
  if (trimester == 1) {
    median = counts_all %>% filter(ga<=12, delivery >= 37)  %>% group_by(bin) %>% select(-starts_with('dCD'), -starts_with('delivery')) %>% summarize_if(is.numeric, funs(median(., na.rm= TRUE)))
  } else if (trimester == 2) {
    median = counts_all %>% filter(ga<=24 & ga>12, delivery >= 37) %>% group_by(bin) %>% select(-starts_with('dCD'), -starts_with('delivery')) %>% summarize_if(is.numeric, funs(median(., na.rm= TRUE)))
  } else { #T3
    median = counts_all %>% filter(ga>24, delivery >= 37)%>% group_by(bin) %>% select(-starts_with('dCD'), -starts_with('delivery')) %>% summarize_if(is.numeric, funs(median(., na.rm= TRUE)))
  }
  return(median)
}

#Calculate MoM gene values per trimester
mom_calc = function(counts, trimester, median_counts_trimester, isClassification) {
  if (trimester == 1) {
    counts_t = counts %>% filter(ga<=12)
  } else if (trimester == 2) {
    counts_t = counts %>% filter(ga<=24 & ga>12)
  } else {
    counts_t =  counts %>% filter(ga>24)
  }
  
  mom = data.frame()
  for (i in unique(counts_t$bin)) {
    med = t(median_counts_trimester %>% filter(bin == i) %>% select(-starts_with("GAPDH")))
    med = as.numeric(med[-c(1, length(med))])
    mom_bin = sweep(counts_t %>% filter(bin == i) %>% select(-starts_with('sample_id'), -starts_with('dCD'), -starts_with('delivery'), -starts_with('ga'),
                                                             -starts_with('bin'), -starts_with('panel'), -starts_with('cohort'), -starts_with('patient')), 2, med, '/')
    
    mom_bin= cbind(data.frame(sample_id = counts_t %>% filter(bin==i) %>% select(sample_id),
                               cohort = counts_t %>% filter(bin==i) %>% select(cohort),
                               ga = counts_t %>% filter(bin==i) %>% select(ga),
                               delivery = counts_t %>% filter(bin==i) %>% select(delivery)), 
                   mom_bin)
    mom = rbind(mom, mom_bin)
  }
  return(mom)  
}

#Wrapper to calculate mom across trimesters
calc_mom_wrapper = function(counts_all) {
  counts_all = filter_bin_data(counts_all)
  
  median_t1 = calc_median_trimester(counts_all, 1)
  median_t2 = calc_median_trimester(counts_all, 2)
  median_t3 = calc_median_trimester(counts_all, 3)
  
  mom_counts = mom_calc(counts_all, 1, median_t1)
  mom_counts = rbind(mom_counts, mom_calc(counts_all, 2, median_t2))
  mom_counts = rbind(mom_counts, mom_calc(counts_all, 3, median_t3))
  
  return(mom_counts)
}

#Filters for matched Danish samples
danish_matched_sample = function(mom_counts) {
  mom_danish = mom_counts %>% filter(cohort == "Danish")
  mom_danish$patient = gsub("\\_[0-9]*$", "", mom_danish$sample_id)
  
  danish_matched_sample = mom_danish %>% filter(ga>24) %>% group_by(patient) %>% summarize_if(is.numeric, funs(mean(.,na.rm=TRUE)))
  danish_matched_sample$patient = gsub('$', '_3', danish_matched_sample$patient)

  danish_matched_sample$cohort = 'Danish'
  colnames(danish_matched_sample)[1] = 'sample_id'
  mom_counts = rbind(mom_counts %>% filter(cohort != 'Danish'), danish_matched_sample)
  return(mom_counts)
}

calc_mom_avg_wrapper = function(counts_all) {
  mom_counts = calc_mom_wrapper(counts_all)
  mom_counts = danish_matched_sample(mom_counts) %>% select(-starts_with('ACTB')) #housekeeping
  return(mom_counts)
}