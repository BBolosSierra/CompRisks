## Calculate the censoring probabilities
## Using prodlim
## Controlling for competing risks
## No covariates in this function
censor.prob.KM <- function(time, status, cens.code){
  tmp <- data.frame(time=time)
  # Sets the value in status equal to cens.code to 0, the rest to 1
  # Controlling for competing risks in censoring status
  tmp$censor.status <- ifelse(status == cens.code, 0, 1) # Not handled in riskRegression with cmprsks
  time.max = ceiling(max(time))
  # Fit prodlim, reversed non-parametric survival KM
  # Switches the censoring status to estimate censoring distribution
  fit = prodlim::prodlim(formula=Surv(time,censor.status)~1,
                         data=tmp,
                         reverse=TRUE)
  # Predict weights at specific times
  prob = stats::predict(fit,
                        times=seq(0,time.max,0.1), # switch for: unique(time) ?
                        level.chaos=1,
                        mode="matrix",
                        type="surv")
  
  # Predict weights at subject specific times
  #ipcw.subject.times = prodlim::predictSurvIndividual(fit,lag=1)
  
  #out <- list(ipcw.times=ipcw.times, ipcw.subject.times=ipcw.subject.times)
  
  out <- cbind(seq(0,time.max,0.1), prob)
  
  out <- na.omit(out)
  
  return(out)
}




# For integrated Brier score
trapezoidal.integration <- function(x, y) {
  if (length(x) != length(y)) stop("x and y must have the same length")
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}