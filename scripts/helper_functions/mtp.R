# https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/#:~:text=Minimum%20training%20presence,-This%20threshold%20finds&text=Essentially%2C%20it%20assumes%20that%20the,area%20of%20the%20binary%20model.

sdm_threshold <- function(sdm, occs, type = "mtp", threshold = 0.1, binary = FALSE){
  th_val = 1-threshold
  occPredVals <- terra::extract(sdm, occs)
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals[,1]))
  } else if(type == "percentile"){
    if(length(occPredVals[,1]) < 10){
      p10 <- floor(length(occPredVals[,1]) * th_val)
    } else {
      p10 <- ceiling(length(occPredVals[,1]) * th_val)
    }
    thresh <- rev(sort(occPredVals[,1]))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- 0
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}
