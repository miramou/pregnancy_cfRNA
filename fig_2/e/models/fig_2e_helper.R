estimate_delivery = function(preds, trimester) {
  test_est = data.frame()
  if (trimester == 2) {
    to_test = preds %>% filter(ga>12 & ga <=24, obs <=0)
  } else if (trimester == 3) {
    to_test = preds %>% filter(ga>24, obs <=0)
  } else {
    to_test = preds %>% filter(ga>12, obs<=0) #T2 & T3
  }
  
  for (p in unique(to_test$patient)) {
    #Get patient data
    pat_test = to_test %>% filter(patient == p)
    
    #Find min gestational age as anchoring pont
    earliest = min(pat_test$ga)
    
    #Find delta time from earliest point to any other sampling point
    pat_test$delta = pat_test$ga - earliest
    
    #Adjust prediction accordingly
    pat_test$adj_pred = pat_test$pred - pat_test$delta
    
    #Find median of all predictions
    pred = median(pat_test$adj_pred[order(pat_test$ga)]) 
    
    #Calc gestational age for avg prediction
    gest_age = mean(pat_test$ga[order(pat_test$ga)])
    sd_gest_age = sd(pat_test$ga[order(pat_test$ga)])
    
    #Put data together in temp dataframe
    temp = data.frame(patient = p, pred= pred, obs = pat_test$obs[pat_test$ga==earliest], 
                      est_dCD = pat_test$est_dCD[pat_test$ga==earliest], ga_test = gest_age, 
                      sd_test = sd_gest_age)
    
    #Add to existing dataframe
    test_est = rbind(test_est, temp)
  }
  
  #Difference between observed time to delivery and predicted for cfRNA
  test_est$delta = test_est$obs-test_est$pred
  
  #Difference between observed time to delivery and predicted for ultrasound
  test_est$us_delta = test_est$obs-test_est$est_dCD
  return(test_est)
}

get_stats = function(trimester_df, trimester_raw_data, type) {
  stats = c(dim(trimester_df %>% filter(condition == type, value < -2) %>% select(condition))[1]/(dim(trimester_raw_data)[1]), #Less than 2 weeks
                  dim(trimester_df %>% filter(condition == type, value >= -2 & value < -1) %>% select(condition))[1]/(dim(trimester_raw_data)[1]), #Between -1 to -2 weeks
                  dim(trimester_df %>% filter(condition == type, value <= 1 & value >-1) %>% select(condition))[1]/(dim(trimester_raw_data)[1]), #Between -1 to 1 weeks
                  dim(trimester_df %>% filter(condition == type, value <= 2 & value >1) %>% select(condition))[1]/(dim(trimester_raw_data)[1]), #Between 1 to 2 weeks
                  dim(trimester_df %>% filter(condition == type, value > 2) %>% select(condition))[1]/(dim(trimester_raw_data)[1])) #More than 2 weeks
  return(stats)
}