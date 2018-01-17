#Calculate centered moving average

remove_na = function(dCt_cma, whichDim, cutoff) {
  to_remove = data.frame()
  for (vector in 1:dim(dCt_cma)[whichDim]) {
    if (whichDim == 1) {
      test= sum(is.na(dCt_cma[vector,]))
    } else {
      test = sum(is.na(dCt_cma[,vector]))
    }
    
    if (test > cutoff) {
      to_remove = rbind(to_remove, vector)
    }
  }
  
  if (whichDim == 1) {
    dCt_cma = dCt_cma[-to_remove[,1],]
  } else {
    dCt_cma = dCt_cma[,-to_remove[,1]]
  }
  
  return(dCt_cma)
}

centered_moving_avg = function(dCt_gene, order) {
  if (order%%2==0) {
    backstep = order/2 - 1
    frontstep = order/2
  } else {
    backstep = floor(order/2)
    frontstep = floor(order/2)
  }

  dCt_modified = vector()
  for (idx in (1+backstep):(length(dCt_gene)-frontstep)) {
    dCt_modified[idx] = mean(dCt_gene[(idx-backstep):(idx+frontstep)], na.rm=TRUE)
  }
  dCt_modified[1:backstep] = NA
  dCt_modified[(length(dCt_gene)-frontstep+1):length(dCt_gene)] = NA
  return (dCt_modified) 
}

cma = function(dCt_all) {
  dCt_all$patient = gsub("\\_[0-9]*$", "", dCt_all$sample_id)
  dCt_no_pp = dCt_all %>% filter(dCD <= 0)

  cma = dCt_no_pp %>% group_by(patient) %>% mutate_if(is.numeric, funs(centered_moving_avg(., order= 3)))
  cma = remove_na(cma, 1, 60)
  return(cma)
}