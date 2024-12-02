#'
#' Title: Brier score Binary
#' 
#' Brier score is the overall prediction error, or the average squared difference between estimated risks (predicted probabilities) and actual outcomes.
#' Numeric values ranging from 0 (perfect predictions) to 1 (maximum error), where lower scores indicate better model performance. 
#'
#' \deq{BS = \frac{1}{N} \sum_{i=1}^{N}(\hat{y}_{i}(t)-y_{i}({t}))^2}
#' 
#' Where N is the number of observations,
#' \deq{\hat{y}_{i}(\tau)} is the predicted probability of the event at time for individual i by the model. 
#' \deq{y_{i}(\tau)} is the event outcome is the actual event outcome (1 if the event occurred, 0 if it did not).
#' 
#' @param predictions: vector of model predictions.
#' @param time time vector of times. Time-to-event.
#' @param status vector of events. Binary [0,...,1]
#' @param tau evaluation time of interest.
#' @param cause event of interest.
#' 
#' @return Numeric value, overall prediction error ranging from 0 to 1, where lower scores indicate better performance at specific time tau of interest.
#' 

#BinaryBrierScore <- function(predictions,
#                        time, 
#                        status, 
#                        tau, 
#                        cens.code){ 
#  
#  y_true = ((time <= tau) * (status != cens.code))
#  
#  BS = mean((predictions-y_true)^2)
#  
#  return(BS)
#  
#}

BinaryBrierScore <- function(predictions,
                             time, 
                             status, 
                             tau, 
                             cause){ 
  
  y_true = ((time <= tau) * (status == cause))
  
  BS = mean((predictions-y_true)^2)
  
  return(BS)
  
}



