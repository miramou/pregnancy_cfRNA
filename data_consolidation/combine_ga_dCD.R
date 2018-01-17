combineGA_dCD = function(ga_datafile, object, isDanish, counts_danish) {
  ga_data = read.csv(ga_datafile, header = TRUE, strip.white = TRUE)
  ga_data$dCD = ga_data$ga_collection - ga_data$ga_delivery #dCD = delta time between sample collection-delivery
  
  indices = match(object$sample_id, ga_data$sample_id)
  object$dCD= ga_data$dCD[indices]
  object$ga = ga_data$ga_collection[indices]
  object$delivery = ga_data$ga_delivery[indices]
  
  if (isDanish) {
    indices_danish = match(counts_danish$sample_id, object$sample_id)
    object$dCD[indices_danish] = counts_danish$dCD
    object$ga[indices_danish] = counts_danish$gestation_age
    object$delivery[indices_danish] = counts_danish$gestation_age - counts_danish$dCD
  }
  return(object)
}