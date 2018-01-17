calc_copy_num_denoised_wrapper = function(data_file, ntc_data_file, whichCohort, 
                                          sweepStart = 6, plasma_info) {
  
  require(reshape2)
  require(dplyr)
  
  samples = read.csv(data_file, header = TRUE, strip.white = TRUE)
  reps_denoised = techRepDenoised_count(samples, ntc_data_file, whichCohort, sweepStart)
  copy_num_per_mL = change_to_num_copies_per_mL(reps_denoised, sweepStart+1, plasma_info, whichCohort)
  
  return(copy_num_per_mL)
}

techRepDenoised_count = function(samples, ntc_data_file, whichCohort, sweepStart) {
  
  #Use experimental cutoff when available
  ntc = read.csv(ntc_data_file, header = TRUE, strip.white = TRUE)
  ntc_limit = ntc %>% summarise_if(is.numeric, funs(median(., na.rm=TRUE)))
  ntc_limit = melt(ntc_limit)
  
  if (whichCohort == 'Danish') {
    reps_long = melt(samples, id = c('sample_id', 'dilution', 'vol_input', 'patient', 'gestation_age', 'time_to_delivery'))
  } else {
    reps_long = melt(samples, id = c('sample_id', 'dilution', 'vol_input', 'cohort', 'sample_number'))
  }
  
  reps_denoised = data.frame()
  
  #Remove noise - either mean for NTC available or hard cutoff of 22.  
  for (gene in colnames(samples)[sweepStart:length(samples)]) {
    if (gene %in% ntc_limit$variable) {
      limit = ntc_limit$value[ntc_limit$variable==gene]
      if (limit == 999) {
        limit = 22
      }
    } else {
      limit = 22
    }
    reps_gene = reps_long %>% filter(variable == gene) %>% mutate(value= ifelse(value >= limit, 999, value))
    reps_denoised = rbind(reps_denoised, reps_gene)
  }
  
  if (whichCohort == 'Danish') {
    samples = dcast(reps_denoised,sample_id+dilution+vol_input+patient+gestation_age+time_to_delivery ~ variable)
  } else {
    samples = dcast(reps_denoised, sample_id+dilution+vol_input+cohort+sample_number~variable)
  }
  
  return(samples)
}

ercc_ct_num_copy = function(x) {
#See supplemental figures for details
  slope = -1.1228
  intercept = 22.6323
  return((x-intercept)/slope)
}

correct_for_mL_input = function(Ct_copy_num, plasma_info_file, Ct_info) {
  plasma_info = read.csv(plasma_info_file, header=TRUE, strip.white = TRUE)
  indices = match(Ct_info$sample_id, plasma_info$sample_id)
  plasma_mL = Ct_info$dilution*12/Ct_info$vol_input/plasma_info$plasma_mL[indices] #12 uL indicates original elution volume
  
  Ct_copy_num = sweep(Ct_copy_num, 1, plasma_mL, '*')
  
  return(Ct_copy_num)
}

change_to_num_copies_per_mL = function(Ct, sweepStart, plasma_info, whichCohort) {
  Ct_copy_num = ercc_ct_num_copy(Ct[,sweepStart:length(Ct)])
  Ct_copy_num = Ct_copy_num %>% mutate_if(is.numeric, funs(2^.))
  
  Ct_copy_num = correct_for_mL_input(Ct_copy_num, plasma_info, Ct[,1:sweepStart-1])
  Ct_copy_num = cbind(Ct[,1:sweepStart-1], Ct_copy_num)
  names(Ct_copy_num)[names(Ct_copy_num) == 'time_to_delivery'] = 'dCD'
  
  
  if (whichCohort == 'Danish') {
    Ct_copy_num = Ct_copy_num %>% filter(vol_input == 1, dilution == 5) %>% 
      select(-starts_with('dilution'), -starts_with('vol_input')) #10 uL input appears to saturate.
  } else if (whichCohort == 'UPenn_UAB') {
    Ct_copy_num = Ct_copy_num %>% filter(vol_input == 2, dilution == 5) %>% 
      select(-starts_with('dilution'), -starts_with('vol_input')) 
  }
  
  return(Ct_copy_num)
}
