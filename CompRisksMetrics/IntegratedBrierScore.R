#'
#' Title: Integrated Brier Score
#' 
#' The IBS is a single, time-averaged value summarizing the prediction error across the entire time range.
#' For survival cases and competing risks cases.
#' It is computed by integrating the Brier Score area under the curve using the trapezoidal rule and normalizing by the length of the time interval.
#' Numeric from 0 to 1, lower values indicate better performance
#'
#' \deq{IBS = \int^{\tau_{max}}_{\tau_{min}} BS(\tau) d\tau}
#' 
#' Can be numerically approximated using trapezoidal rule:
#' 
#' \deq{IBS \approx \frac{1}{\tau_{\text{max}} - \tau_{\text{min}}} \sum_{j=1}^{K-1} \frac{\tau_{j+1} - \tau_j}{2} \cdot \left( BS(\tau_j) + BS(\tau_{j+1}) \right)}
#'
#' Function trapezoidal.integration to calculate IBS.
#' 
#' @param prediction_matrix matrix with predictions from any model at a range of times of interest, where rows are individuals and columns each time of evaluation
#' @param taus vector with a range of evaluation times i.e seq(0,5,0.1)
#' @param time time vector of times. Time-to-event.
#' @param status vector of events. 1) If survival -> [0,...,1]. 2) If competing risks -> [0,...,K].
#' @param cause event of interest. 
#' @param cens.code value of censoring status, commonly censor patients have status of 0.
#' @param cmprsk logical vector for the presence of competing risks. If TRUE, there are competing risks, otherwise binary outcomes.
#' 
#' @return 
#'  \itemize{
#'   \item integrated.brier.score: Numeric value, prediction error across the entire time range ranging from 0 to 1, where lower scores indicate better performance.
#'   \item weighted.brier.score: vector of each weighted brier score at each of the times of evaluation (taus)
#'   \item taus: vector of times of evaluation
#' }



IntegratedBrierScore <- function(prediction_matrix,
                                 taus,
                                 time,
                                 status,
                                 cause, 
                                 cens.code, 
                                 cmprsk){
  
  # Number of times evaluated the brier score
  nt <- dim(prediction_matrix)[2]
  # Initialize a empty vector for brier scores
  bs <- rep(0, nt)
  # Loop through array
  # Rows: individuals, Columns: time points
  for (j in seq(nt)) {
    # Get the model prediction and time of evaluation tau
    tau <- taus[j]
    predictions <- prediction_matrix[,j] 
    # Calculate Weighted Brier Score
    BS_we <- WeightedBrierScore2(predictions = predictions,
                       tau = tau,
                       time =  time,
                       status = status,
                       cause = cause, 
                       cens.code = cens.code, 
                       cmprsk = cmprsk)
    # Store value 
    bs[j] <- BS_we
  }
  # Calculate max and min time to normalize with respect to time interval
  t.max = max(taus)
  t.min = min(taus)
  # Apply trapezoidal integration
  ibs <- trapezoidal.integration(bs, taus)/(t.max - t.min)
  
  return(list(integrated.brier.score = ibs,
              weighted.brier.score = bs,
              taus = taus))
}