#'
#' Title: Calibration Numerical Measures At Multiple Times
#' 
#' This function calculates numerical calibration measures computed from the absolute error (difference) between the approximation to observed values and predicted. 
#' \deqn{|y_{i} - \hat{y}_{i}|}
#' This function also calculates the squared root of the mean squared error: \deq{\text{RMSE} = \sqrt{\frac{1}{n} \sum_{i=1}^{n} \left( y_i - \hat{y}_i \right)^2}}
#' The input are predicted survival probabilities (predictions.type = 'survival') or Cumulative Incidence Function (CIF) (predictions.type = 'CIF') from any model.
#' 
#' In this function the observed approximation values are calculated with function CalibrationPlots with the option of choosing between different approximation methods:
#' The available methods are 'pseudovalues', 'quantiles' or 'subdistribution'. For more details see CalibrationPlots function.
#' 
#' The function outputs the following calibration measures:
#' 
#' 1) Integrated Calibration Index (ICI) is the mean absolute difference between the observed approximation and predicted.
#' It is the same as the mean absolute error (MAE) and represents the average error.
#' A lower ICI value indicates better calibration, meaning the predicted values are closer to the observed on average. 
#' 
#' 2) 50th Percentile of Absolute Error (E50) is the median absolute error between the observed approximation and predicted. 
#' Represents the typical or central error that is more robust to skewed error distributions or with outliers.
#' 
#' 3) 90th Percentile of Absolute Error (E90) is the upper quantile of the error between observed approximation and predicted. 
#' Represents the how large errors, the upper 10%, behave. If the upper 10% of error has a E90 close to 0, can be interpreted as if even large
#' errors are relatively small. 
#'
#' 4) Maximum Absolute Error (Emax) is the bigger error between observed approximation and predicted. 
#' Represents the highest error, and a low value can be interpreted as the model not having extreme miscalibrations.
#' 
#' 5) Squared Root of the Mean Squared Error (RMSE) between observed approximations and the predicted.
#' It will penalize larger errors more heavily than the smaller errors. Can be used as overall prediction error. 
#' \deq{\text{RMSE} = \sqrt{\frac{1}{n} \sum_{i=1}^{n} \left( y_i - \hat{y}_i \right)^2}}
#' 
#' @param predictions matrix of model predictions. 
#' @param prediction.type string with "CIF" if cumulative incidence functions or "survival" if survival predictions
#' @param data dataset.
#' @param time time vector of times. Time-to-event.
#' @param status vector of events. 1) If survival -> [0,...,1]. 2) If competing risks -> [0,...,K].
#' @param tau evaluation time of interest or vector with a range of times [0, max(tau)].
#' @param cause event of interest.
#' @param cens.code value of censoring status, commonly censor patients have status of 0.
#' @param method method chosen to obtain aproximated observation values
#' @param quantiles number of quantiles for predicted risk intervals. Often deciles.
#' @param bandwidth parameter kernel smoothing with nearest neighbours. Computed if NULL.
#' @param n.knots number of knots for splines. 3-5 from less to more flexibility. Set to 5 (more flexible)
#'
#' @return List of numerical measures at times of interest taus.


CalibrationNumericalMeasuresMTimes <- function(predictions, # any model CIF or S predictions
                                               predictions.type = c("CIF", "survival"), 
                                               data, 
                                               time, 
                                               status, 
                                               taus, 
                                               cause=1, 
                                               cens.code=0,
                                               method = c("quantiles", "pseudovalues", "subdistribution"),
                                               quantiles=10,
                                               bandwidth = NULL,
                                               n.knots = 5){

  values <- list()
  error  <- list()
  measures <- list()
  for (t in seq_along(taus)){ 
    # Get the observed approximation
    values[[t]] <- CalibrationPlot(predictions[,t], data, time, status, taus[t], 
                                    cause, cens.code, predictions.type, 
                                    method, bandwidth, n.knots, graph = FALSE)
    
    
    error[[t]] <- values[[t]]$pred - values[[t]]$obs
    
    # Compute the numerical measures
    measures[[t]] <- c(
      "ICI" = mean(abs(error[[t]])), 
      setNames(quantile(abs(error[[t]]), c(0.5, 0.9)), c("E50", "E90")),
      "Emax" = max(abs(error[[t]])),
      "Root squared bias" = sqrt(mean(error[[t]]^2))
    )
  }
  
  return(result = list(taus=taus, values=values, error=error, calibr.measures=measures))
}

