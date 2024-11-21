
#' Title: Calibration Plot

#' The function aims to compute calibration plot with direct use of survival predictions, hazard or cumulative incidence function predicted from any model.
#' The function works for right censored time-to-event with and without competing risks.
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
#' @param predictions vector of model predictions. 
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
#' @return List if graph == TRUE. Otherwise return values list with input predictions and the approximation to observed values
#'  \itemize{
#'   \item graph: ggplot of calibration
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
CalibrationPlot <- function(predictions, # any model CIF or S predictions
                                      data, 
                                      time, 
                                      status, 
                                      tau, 
                                      cause=1, 
                                      cens.code=0,
                                      predictions.type = c("CIF", "survival"), # only ready for cmprsk
                                      method = c("quantiles", "pseudovalues", "subdistribution"),
                                      #kernel.smoothing = TRUE,
                                      bandwidth = NULL,
                                      quantiles=10,
                                      n.knots = 5,
                                      graph = TRUE) {
  
  if (missing(predictions.type) || is.null(predictions.type)) { stop("Error. Enter the type of predictions 'CIF' or 'survival'")}
  if (predictions.type == "survival") {predictions = 1-predictions}
  
  # 1) Pseudovalue method
  if (method == "pseudovalues"){
    # 1.1 Calculate pseudovalues with prodlim
    margForm <- prodlim::Hist(time, status,
                              cens.code=cens.code)~1
    margFit  <- prodlim::prodlim(margForm,data=data)
    
    pseudo   <- prodlim::jackknife(margFit,
                                   cause=cause,
                                   times=tau)
    
    # Alternatively with pseudoci
    #cause_col <-  paste0("cause", cause)
    #pseudo <- pseudo::pseudoci(time, status, tmax=tau)$pseudo[[cause_col]]
    
    # 1.2 Nearest neighbors kernel smoothing with prodlim
    #predictions<- round(p,2)
    predictions <- na.omit(predictions)
    bw <- bandwidth
    if (is.null(bandwidth)) {
      bw <- prodlim::neighborhood(predictions, bandwidth)$bandwidth
    }
    nbh <- prodlim::meanNeighbors(x=predictions,y=pseudo,bandwidth=bw)
    plotFrame <- data.frame(pred=nbh$uniqueX,obs=nbh$averageY)
    # Next step: Alternative for neareast neighbours ?
  
  # 2) Quantiles method  
  } else if (method == "quantiles") {
    # Create the 10 levels of intervals for predicted risk
    pcut <- cut(predictions,quantiles,include.lowest=TRUE)
    # Calculate the mean per each of the 10 levels
    Pred=tapply(predictions,pcut,mean)
    # Levels as dataset
    newdata=data.frame(pcut=levels(pcut))
    # Fit prodlim. Discrete predictor variable: pcut 
    qfit <- prodlim::prodlim(formula = prodlim::Hist(time,
                                                     status,
                                                     cens.code) ~ pcut, data = data)
  
    # Observed risk 
    obs = stats::predict(qfit, newdata = data.frame(pcut=levels(pcut)), 
                  cause=cause,
                  mode="matrix",
                  times=tau,
                  type="cuminc") ### what if survival?
    
    # Dataframe with both observed and predicted
    plotFrame=data.frame(pred=tapply(predictions,pcut,mean),
                         obs=obs)
    
  # 3) Subdistribution hazard method
  } else if (method == "subdistribution"){
    # Complementary log-log
    cll_pred        <- log(-log(1 - predictions))
    # Splines
    basis           <- splines::ns(cll_pred, df = n.knots + 1)
    colnames(basis) <- paste0("basisf_", colnames(basis))
    cols            <- colnames(basis)
    basis           <- cbind.data.frame(basis, time)
    basis$status    <- status  
    # Calibrate 
    calib_fgr <- riskRegression::FGR(formula = reformulate(response = "prodlim::Hist(time, status)", termlabels = cols),
                             cause = cause,
                             data = basis)
    plotFrame=data.frame(pred=predictions,
                         obs=predict(calib_fgr, times = tau, newdata = basis))
  } else {
    stop("Error. Enter a method to calculate the approximation to observed values. The options are 'pseudovalues', 'quantiles' or 'subdistribution'.")
  }
  
  # Plot and Fit a LOESS smooth curve to the pseudo-values vs predicted CIFs
  if (graph == TRUE) {
    graph <- ggplot(plotFrame, aes(x = pred, y = obs)) +
      geom_smooth(method = "loess", color = "blue") +
      geom_point(alpha = 0.5) +  # Add points for each individual
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Ideal calibration line
      labs(
        x = "Predicted",
        y = paste0("Observed (", method, ")"),
        title =  paste0("Calibration plot for ", predictions.type, " with ", method, " at time ", tau),
      ) +
      theme_minimal()
    
    return(list(graph=graph, values=plotFrame))
    
  } else {
    return(values = plotFrame)
  }
}

