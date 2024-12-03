#'
#' Weighted Brier Score 
#'
#' Ref: https://onlinelibrary.wiley.com/doi/epdf/10.1002/bimj.200610301
#'
#' Weighted overall prediction error, is a measure of the mean squared error between predicted probabilities and actual outcomes, 
#' adjusted for censoring using Inverse Probability of Censoring Weights (IPCW) 
#' Numerical value from 0 to 1 where lower values indicate better prediction performance.
#' 
#' Function censor.prob.KM.individual to calculate the IPCW.
#' 
#' The idea is to give less weight to individuals who are censored earlier, 
#' as they provide less complete information about the event of interest. 
#' Those who are not censored or censored later during follow up are given higher weights.
#' 
#' @param predictions vector of model predictions.
#' @param time time vector of times. Time-to-event.
#' @param status vector of events. 1) If survival -> [0,...,1]. 2) If competing risks -> [0,...,K].
#' @param tau evaluation time of interest or vector with a range of times [0, max(tau)].
#' @param cause event of interest.
#' @param cens.code value of censoring status, commonly censor patients have status of 0.
#' @param cmprsk logical vector for the presence of competing risks. If TRUE, there are competing risks, otherwise binary outcomes.
#'
#' @return Weigthed Brier Score
#' 
#' 
#' @examples
#' 


WeightedBrierScore2 <- function(predictions, 
                               time, 
                               status, 
                               tau, 
                               cause, 
                               cens.code, 
                               cmprsk = FALSE) {
  # Calculate individual censoring weights
  ind.cens.weights <- censor.prob.KM.individual(time = time, 
                                                status = status, 
                                                cens.code = cens.code)
  
  # Identify the binary outcome for the event of interest
  y_true <- ((time <= tau) & (status == cause))
  
  # Initialize vector for residuals
  residuals <- numeric(length(predictions))
  
  if (cmprsk) {
    ## Consider censoring with competing risks
    # Case 1: Event of interest occurs before tau
    index <- (time <= tau & status == cause)
    residuals[index] <- ((1 - predictions[index])^2) / ind.cens.weights[index]
    
    # Case 2: Competing event occurs before tau
    index <- (time <= tau & status != cause & status != cens.code)
    residuals[index] <- (predictions[index]^2) / ind.cens.weights[index]
    
    # Case 3: Censored before tau
    index <- (time <= tau & status == cens.code)
    residuals[index] <- 0
    
    # Case 4: Event or censoring after tau
    index <- (time > tau)
    residuals[index] <- (predictions[index]^2) / ind.cens.weights[index]
    
  } else {
    ## Consider censoring without competing risks
    # Case 1: Event occurs before tau
    index <- (time <= tau & status == cause)
    residuals[index] <- ((1 - predictions[index])^2) / ind.cens.weights[index]
    
    # Case 2: Censored before tau (including competing risks treated as censored)
    index <- (time <= tau & status == cens.code)
    residuals[index] <- 0
    
    # Case 3: Event or censoring after tau
    index <- (time > tau)
    residuals[index] <- (predictions[index]^2) / ind.cens.weights[index]
  }
  
  # Compute the Brier Score
  WeightedBrierScore <- mean(residuals)
  
  return(WeightedBrierScore)
}
