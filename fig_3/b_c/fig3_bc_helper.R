calc_eff_size_test = function(counts_long) {
  stats_df = data.frame()
  for (gene in unique(counts_long$variable)) {
    gene_values = counts_long %>% filter(variable == gene) %>% select(value, Class)
    effect_size = cohen.d(value ~ Class, hedges.correction = TRUE, data = gene_values)
    
    fold_change = 2.5 #threshold for fisher's exact
    
    #Get number of samples for sPTB and full-term that are above or below threshold
    p_p_greater = dim(counts_long %>% filter(variable == gene, Class == 'pre', value > fold_change))[1] #num pre-term > FC
    p_p_less = dim(counts_long %>% filter(variable == gene, Class == 'pre', value <= fold_change))[1] #num pre-term <= FC
    t_p_greater = dim(counts_long %>% filter(variable == gene, Class == 'term', value > fold_change))[1] #num term > FC
    t_p_less = dim(counts_long %>% filter(variable == gene, Class == 'term', value <= fold_change))[1] #num term <= FC
    
    #First argument is 2D contingency table, 2nd argument is specified test
    fisher_test = fisher.test(rbind(c(p_p_greater, t_p_greater), c(p_p_less, t_p_less)), alternative = 'greater')
    
    #Add test for each gene to total matrix
    temp = data.frame(gene = gene, estimate = effect_size$estimate, lower_ci = effect_size$conf.int[1], upper_ci = effect_size$conf.int[2], mag = effect_size$magnitude,
                      p_val = fisher_test$p.value)
    stats_df = rbind(stats_df, temp)
  }
  
  #Adjust pval using benjamani hochberg
  stats_df$p_val_adj = p.adjust(stats_df$p_val, method = 'BH')
  
  return(stats_df)
}

plot_PT_T = function(mom_counts, genesToPlot, file_header) {
  
  group.colors = c("term" = "#33a02c", "pre" = "#386cb0")
  
  #Plot dCt per gene vs collection time
  require(ggplot2)
  
  mom_counts_long = melt(mom_counts, id= c('sample_id', 'cohort', 'ga', 'delivery', 'Class'))
  mom_counts_long$cohort = factor(mom_counts_long$cohort, c('T_Penn',  'Danish', 'T_UAB', 'PT_Penn','PT_UAB'))
  
  #Plot and save each gene separately
  for (gene in genesToPlot) {
    toPlot = mom_counts_long %>% filter(variable == gene)
    Plot = ggplot(toPlot, aes(x=cohort, y=log2(value))) +
      geom_boxplot(aes(fill = Class)) +
      scale_fill_manual(values = group.colors) +
      geom_point() +
      labs(x = 'Cohort', y = expression('Multiple of the median'))+
      theme_pub() +
      theme (
        axis.text.x = element_text(angle = 90)
      ) +
      scale_y_continuous(limits = c(-4, 8)) 
    
    pdf(paste0(file_header, gene, '_boxplot.pdf'), useDingbats = FALSE, width = 5, height= 5)
    print(Plot)
    dev.off()
  }
}

get_combos = function(stats, counts) {
  
  combo_results = data.frame()
  fold_change = 2.5
  for (i in 2:4) {
    test_combos = combn(stats$gene, i)
    for (j in 1:dim(test_combos)[2]) {
      test_set = counts[, colnames(counts) %in% c(as.character(test_combos[,j]), 'cohort', 'Class')]
      test_set$pred = ifelse(rowSums(test_set %>% select(-starts_with('cohort'), -starts_with('Class')))>=i*fold_change, 1, 0)
      
      temp = data.frame(combo_num =i, col_num = j, tp = sum(test_set %>% filter(Class == 'pre') %>% select(pred)), fp = sum(test_set %>% filter(Class == 'term') %>% select(pred)))
      combo_results = rbind(combo_results, temp)
    }
  }
  return(combo_results)
}
  
get_results = function (gselect, stats, counts) {
  gselect_results = data.frame()
  fold_change = 2.5 #threshold same as fisher's exact
  
  for (z in 1:dim(gselect)[1]) {
    test_combos = combn(stats$gene, gselect$combo_num[z])
    test_set = counts[, colnames(counts) %in% c(as.character(test_combos[,gselect$col_num[z]]), 'Class', 'sample_id')]
    #If all 3 cfRNA transcripts show MOM > 2.5 then classify as sPTB
    test_set$pred = ifelse(rowSums(test_set %>% select(-starts_with('Class'), -starts_with('sample_id')))>= gselect$combo_num[z]*fold_change, 1, 0)
    
    #Generate temp df with all gene names
    for (e in 1:gselect$combo_num[z]) {
      if (e==1 ) {
        temp = data.frame(gene = as.character(test_combos[,gselect$col_num[z]][e]))
      } else  {
        temp = cbind(temp, data.frame(gene = as.character(test_combos[,gselect$col_num[z]][e])))
      }
    }
    
    #Find rows where prediction failed to classify sPTB correctly
    rows = test_set %>% filter(Class == 'pre') %>% select(sample_id, pred) %>% filter(pred == 0) 
    
    #Label sample that was misclassified as full-term
    if (gselect$tp[z] == 6) {
      temp = cbind(temp, data.frame(miss1 = rows$sample_id[1], miss2 = rows$sample_id[2]))
    } 
    gselect_results = rbind(gselect_results, temp)
  }
  
  #Relabel columns
  colnames(gselect_results) = c('gene1', 'gene2', 'gene3', 'miss1', 'miss2')
  
  return(gselect_results)
}

test_results = function (gselect, stats, counts) {
  predict_test = data.frame(counts %>% select(sample_id, cohort, Class))
  fold_change = 2.5
  
  for (v in 1:dim(gselect)[1]) { #For all relevant combinations
    test_combos = combn(stats$gene, gselect$combo_num[v])  #get relevant combination
    test_set = counts[, colnames(counts) %in% c(as.character(test_combos[,gselect$col_num[v]]))]
    #If all 3 cfRNA transcripts show MOM > 2.5 then classify as sPTB
    test_set$pred = ifelse(rowSums(test_set)>= gselect$combo_num[v]*fold_change, 1, 0)
    temp = data.frame(v = test_set$pred)
    predict_test = cbind(predict_test, temp)
  }
  
  colnames(predict_test) = c('sample_id', 'cohort', 'Class', 1:dim(gselect)[1])
  
  #Probability is the fraction of combinations that classify the sample as sPTB
  predict_test$prob = rowSums(predict_test[,4:(dim(gselect)[1]+3)])/(dim(gselect)[1])
  
  return(predict_test) 
}

plot_roc = function(curve1, curve2, title) {
  curve_df = data.frame(sensitivities = curve1$sensitivities*100, specificities = (1-curve1$specificities)*100)
  curve2_df = data.frame(sensitivities = curve2$sensitivities*100, specificities = (1-curve2$specificities)*100)
  
  y_x = data.frame(x = seq(1, 100, 1), y = seq(1,100, 1))
  
  plot = ggplot() +
    geom_line(data = y_x, aes(x=x, y=y), size = 1.3, color = '#bdbdbd', linetype = 'dashed') +
    geom_line(data = curve_df, aes(x=specificities, y = sensitivities), size = 1.3, color = '#6baed6') +
    geom_line(data = curve2_df, aes(x=specificities, y = sensitivities), size = 1.3, color = '#3182bd') +
    labs(x = 'False Positive Rate', y = 'True Positive Rate') +
    theme_pub()
  
  pdf(title)
  print(plot)
  dev.off()
  
}


