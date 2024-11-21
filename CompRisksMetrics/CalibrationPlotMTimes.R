
#' Title: Calibration Plots at Multiple Times

#' The function aims to compute calibration plot with direct use of survival predictions, hazard or cumulative incidence function predicted from any model.
#' The function works for right censored time-to-event with and without competing risks.
#' 
#' Uses the function CalibrationPlots at different time points of interest.
#' 
#' Method 1) Pseudovalue calculation with prodlim (alternatively with pseudoci):  
#'  Step 1.1 Pseudo-observation approximate the individual contribution to the overall CIF or S estimate by using leave-one-out approach:
#'  \deqn{Pseudovalue_{i}(tau, k) = n \hat{CIF}_{k}(tau) - (n-1)\hat{CIF}_{k}^{(-i)}(tau)} 
#'  \deqn{Pseudovalue_{i}(tau, k) = n \hat{S}_{k}(tau) - (n-1)\hat{S}_{k}^{(-i)}(tau)}
#'  First, the CIF is estimated for all the individuals at a given time point. 
#'  Then, the estimate is recalculated removing one individual at a time. 
#'  These estimations are non-parametric estimations of S(tau) computed with Kaplan-Meier, or CIF(tau) computed with Aalen-Johansen estimator.
#' 
#'  Step 1.2 Nearest neighbors kernel smoothing with prodlim:
#'  If kernel smoothing is choosen, nearest neighbors is computed.
#'  The observed values (pseudovalues) would be the mean of nearest neighbors.
#'  Bandwidth for computation is calculated if NULL.
#' 
#' Method 2) Quantiles. Discretized intervals
#' Quantiles, often deciles create discretized intervals of the predictions
#' Prodlim is then fitted with the discretized intervals.
#' The observed values would be computed per each quantile.
#' 
#' Method 3) Subdistribution hazard with Fine and Grey with splines
#' 
#' 
#' @param predictions matrix of model predictions. 
#' @param data dataset.
#' @param time time vector of times. Time-to-event.
#' @param status vector of events. 1) If survival -> [0,...,1]. 2) If competing risks -> [0,...,K].
#' @param tau evaluation time of interest or vector with a range of times [0, max(tau)].
#' @param cause event of interest.
#' @param cens.code value of censoring status, commonly censor patients have status of 0.
#' @param prediction.type 
#' @param method method chosen to obtain aproximated observation values
#' @param quantiles number of quantiles for predicted risk intervals. Often deciles.
#' @param bandwidth parameter kernel smoothing with nearest neighbours. Computed if NULL.
#' @param n.knots number of knots for splines. 3-5 from less to more flexibility. Set to 5 (more flexible)
#' @param graph TRUE if predicted versus observed values are wanted to be plotted
#' 
#' @return Multiple lists per time of interest in taus.
#'  \itemize{
#'   \item graph: ggplot of calibration if graph == TRUE. 
#'   \item values: list with input predictions (pred) and the approximation to observed values (obs)
#' }
#' 
#' @examples
#' # Method 1) Pseudovalues and kernel smoothing
#' # Method 2) Quantiles
#' # Method 3) Subdistribution hazard
#' 
#' 
#'

CalibrationPlotMTimes <- function(predictions, # any model CIF or S predictions
                                            data, 
                                            time, 
                                            status, 
                                            taus, 
                                            cause=1, 
                                            cens.code=0,
                                            predictions.type = c("CIF", "survival"), # only ready for cmprsk
                                            method = c("quantiles", "pseudovalues", "subdistribution"),
                                            bandwidth = NULL,
                                            quantiles=10,
                                            n.knots = 5,
                                            graph = TRUE) {
  result <- list()
  # Sequence along tau
  for (t in seq_along(taus)){ 
    # Get the plots and values
    result[[t]] <- CalibrationPlot(predictions[,t], data, time, status, taus[t], 
                                    cause, cens.code, predictions.type, 
                                    method, bandwidth, n.knots, graph)
  } 
  
  return(result)
  
}