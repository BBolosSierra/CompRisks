#'
#' Title: Scaled Brier score
#'
#' The scaled Brier score is a normalized version of the Brier score that helps to 
#' interpret the predictive performance of a model relative to a baseline model, 
#' often the null model. 
#' Number from 0 to 1, the higher the better. 
#' It can also be interpreted as the improvement of the evaluated model over the null model.
#' For clarity, the Scaled Brier Score is calculated with the Weighted Brier Score.
#' 
#' Reference: https://arxiv.org/pdf/2212.05157.pdf
#' 
#' \deq{BS_{sc}(\tau) = 1- \frac{BS(\tau)}{BS(\tau)^\text{Null}}}
#' 
#' \deq{BS(\tau)^\text{Null} =\frac{1}{N} \sum_{i=1}^{N} W_{i}(\tau) \cdot \left( D_{i}(\tau) - f_{i}^\text{Null}(\tau) \right)^2}
#' 
#' \deq{f_{i}^\text{Null}(\tau)} is the same for all individuals, as it is derived purely from population-level probabilities at time $\tau$
#' 
#' @param prediction_null Predicted cumulative incidence function or survival probability under the null model (no covariates) 
#' @param predictions vector of model predictions.
#' @param time time vector of times. Time-to-event.
#' @param status vector of events. 1) If survival -> [0,...,1]. 2) If competing risks -> [0,...,K].
#' @param tau evaluation time of interest or vector with a range of times [0, max(tau)].
#' @param cause event of interest.
#' @param cens.code value of censoring status, commonly censor patients have status of 0.
#' @param cmprsk logical vector for the presence of competing risks. If TRUE, there are competing risks, otherwise binary outcomes.
#'
#' @return 
#'  \itemize{
#'   \item weighted.brier.score.null:
#'   \item weighted.brier.score:
#'   \item scaled.brier.score: Numeric value, scaled prediction error ranging from 0 to 1, where higher scores indicate better performance of the evaluated model at specific time tau of interest with respect to the null model.
#'   \item percentage.scaled.brier.score:
#'   \item tau: vector of times of evaluation
#' }
#' 
#' 
#' 

ScaledBrierScore <- function(predictions,
                              predictions_null,
                              tau,
                              time, 
                              status,
                              cause, 
                              cens.code,
                              cmprsk){
  # Number of patients
  n <- length(predictions)
  # Length of predictions null
  pn <- length(predictions_null)
  # If only a single value is given, repeat it for all observations 
  if (pn == 1) {predictions_null <- rep(predictions_null, n)}
  
  BS_weighted_null <- WeightedBrierScore(predictions = predictions_null,
                                            tau = tau,
                                            time = time,
                                            status = status,
                                            cause = cause, 
                                            cens.code = cens.code, 
                                            cmprsk = cmprsk)$weighted.brier.score
  
  BS_weighted <- WeightedBrierScore(predictions = predictions,
                                       tau = tau,
                                       time = time,
                                       status = status,
                                       cause = cause, 
                                       cens.code = cens.code, 
                                       cmprsk = cmprsk)$weighted.brier.score
  
  
  BS_scaled = 1 - (BS_weighted/BS_weighted_null)
  
  return(list(weighted.brier.score.null = BS_weighted_null, 
              weighted.brier.score = BS_weighted, 
              scaled.brier.score = BS_scaled,
              percentage.scaled.brier.score = BS_scaled*100,
              tau = tau))
}
