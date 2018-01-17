avg_wrapper = function(counts, interval, isGA) {
  library(dplyr)
  
  pp = counts %>% filter(dCD>0)
  counts = smooth_intrapatient_avg(counts, interval, isGA)
  counts = interpatient_avg(counts, pp, interval, isGA)
  
  return(counts)
}

smooth_intrapatient_avg = function(counts, interval, isGA) {
  
  counts = counts %>% filter(dCD <= 0) #Remove post-partum samples
 
  #Take average for every interval weeks either for GA or dCD
  if (isGA) {
    weeks = seq(4, 44, by=interval) 
  } else {
    weeks = seq(floor(min(counts$dCD)), 1, by = interval) 
  }
  
  avg = data.frame()

  #Take intrapatient average either for GA or dCD
  for (idx in 2:length(weeks)) {
    if (isGA) {
      temp_avg = counts %>% group_by(patient) %>% filter(ga>weeks[idx-1] & ga<= weeks[idx]) %>% summarize_if(is.numeric, funs(mean(., na.rm=TRUE)))
    } else {
      temp_avg = counts %>% group_by(patient) %>% filter(dCD>weeks[idx-1] & dCD<= weeks[idx]) %>% summarize_if(is.numeric, funs(mean(., na.rm=TRUE)))
    }
    avg = rbind(avg, temp_avg)
  }
  return(avg)
}

interpatient_avg = function(counts, pp, interval, isGA) {
  library(reshape2)
  
  all_patient_avg = data.frame()
  all_patient_sd = data.frame()
  sem_n = c() #the n for every month's average to calculate SEM
  
  #Address post-partum separately
  pp_mean = pp %>% summarize_if(is.numeric, funs(mean(.,na.rm=TRUE)))
  pp_sd = pp %>% summarize_if(is.numeric, funs(sd(.,na.rm=TRUE)))
  pp_sd = pp_sd/sqrt(dim(pp)[1]) #sem
  
  #Take average for every interval weeks either for GA or dCD
  if (isGA) {
    weeks = seq(4, 44, by=interval)
  } else {
    weeks = seq(floor(min(counts$dCD)), 1, by = interval) 
  }
  
  for (idx in 2:length(weeks)) {
    if (isGA) {
      temp = counts %>% filter(ga>weeks[idx-1] & ga<= weeks[idx])
    } else {
      temp = counts %>% filter(dCD>weeks[idx-1] & dCD<= weeks[idx])
    }
    sem_n = c(sem_n, dim(temp)[1])
    all_patient_avg = rbind(all_patient_avg, temp %>% summarize_if(is.numeric, funs(mean(., na.rm=TRUE))))
    all_patient_sd =  rbind(all_patient_sd, temp  %>% summarize_if(is.numeric, funs(sd(., na.rm=TRUE))))
  }
  
  sem_n = 1/sqrt(sem_n)
  pp_mean$ga = round(pp_mean$ga) 
  pp_mean$dCD = round(pp_mean$dCD)
  
  all_patient_sd = sweep(all_patient_sd %>% select(-starts_with('ga'), -starts_with('dCD'), -starts_with('delivery')), 1, sem_n, '*')
  all_patient_avg = rbind(all_patient_avg, pp_mean)
  all_patient_sd = rbind(all_patient_sd, pp_sd %>% select(-starts_with('ga'), -starts_with('dCD'), -starts_with('delivery'), -starts_with('GAPDH')))
  
  #Add SEM column to Avg dataframe
  all_patient_avg = melt(all_patient_avg %>% select(-starts_with('GAPDH')), id = c('ga', 'dCD', 'delivery'))
  all_patient_sd = melt(all_patient_sd)
  all_patient_avg$sd = all_patient_sd$value
  
  return(all_patient_avg)
}

plot_panel = function(counts, panel, isGA, plot_name) {

  toPlot = counts[counts$variable %in% panel,]
  toPlot$variable = factor(toPlot$variable, levels = panel)

  row_n = ceiling(length(panel)/7)
  
  if (isGA) {
    plot = ggplot(data = toPlot, aes(x=ga, y= value)) +
      geom_line() +
      geom_point(size = 2.5) +
      geom_errorbar(aes(x=ga, ymin = value-sd, ymax = value+sd)) +
      theme_pub() +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold')) +
      coord_cartesian(xlim = c(0, 45)) +
      expand_limits(y=0) +
      labs(x='Gestational age (weeks)', y=expression(over('Estimated transcript count', 'mL plasma'))) +
      facet_wrap(~variable, nrow = row_n, scales = "free")
    
  } else {
    plot = ggplot(data = toPlot, aes(x=dCD, y= value)) +
      geom_line() +
      geom_point(size = 2.5) +
      geom_errorbar(aes(x=dCD, ymin = value-sd, ymax = value+sd)) +
      theme_pub() +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold')) +
      coord_cartesian(xlim = c(-40, 5)) +
      expand_limits(y=0) +
      labs(x='Time to delivery (weeks)', y=expression(over('Estimated transcript count', 'mL plasma'))) +
      facet_wrap(~variable, nrow = row_n, scales = "free")
  }
  pdf(plot_name, useDingbats = FALSE, width = 20, height = 10)
  print(plot)
  dev.off()
}

