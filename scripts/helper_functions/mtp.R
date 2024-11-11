# https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/#:~:text=Minimum%20training%20presence,-This%20threshold%20finds&text=Essentially%2C%20it%20assumes%20that%20the,area%20of%20the%20binary%20model.

sdm_threshold <- function(sdm, occs, type = "mtp", binary = FALSE){
  occPredVals <- terra::extract(sdm, occs)
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals))
  } else if(type == "p10"){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- 0
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}