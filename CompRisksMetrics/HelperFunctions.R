prepare_data <- function(df, time_col, event_col) {
  time_col  <- deparse(substitute(time_col))
  event_col <- deparse(substitute(event_col))
  
  nm <- names(df)
  nm[nm == time_col]  <- "time"
  nm[nm == event_col] <- "event"
  names(df) <- nm
  
  df$time  <- as.numeric(df$time)
  df$event <- as.integer(df$event)
  
  # Sort by time maybe not necesary 
  df[order(df$time), , drop = FALSE]
}

## this is like DeepHit 
sim_cmprks <- function(n=30000, 
                       p_block = 4,
                       gammas = 10, # y1 = y2 = y3 = 10
                       cens_frac = 0.5, 
                       seed = 123) {
  
  set.seed(seed)
  
  gamma1 = rep(gammas, p_block)
  gamma2 = rep(gammas, p_block)
  gamma3 = rep(gammas, p_block)
  # covariates that are 4 dimensional
  x1 <- matrix(rnorm(n * p_block), n, p_block)
  x2 <- matrix(rnorm(n * p_block), n, p_block)
  x3 <- matrix(rnorm(n * p_block), n, p_block)
  
  # linear effects and matrix multiplications
  g1_x1 <- as.numeric(x1 %*% gamma1)
  g2_x2 <- as.numeric(x2 %*% gamma2)
  g3_x3 <- as.numeric(x3 %*% gamma3)
  
  # quadratic effects 
  lambda1 <- exp(g3_x3**2 + g1_x1)
  lambda2 <- exp(g3_x3**2 + g2_x2)
  
  # hitting times
  T1 = rexp(n, rate=lambda1)
  T2 <- rexp(n, rate = lambda2)
  si <- pmin(T1, T2)
  
  # censoring
  sc <- rep(Inf, n)
  if (cens_frac > 0) {
    cens_idx <- sample.int(n, size = floor(cens_frac * n))
    sc[cens_idx] <- runif(length(cens_idx), 0, si[cens_idx])
  }
  
  # Tobs and cause
  Tobs  <- pmin(T1, T2, sc) 
  cause <- integer(n) # just 30000 zeroes
  is_cens <- sc <= si # if so is censor
  uncens <- !is_cens
  c1 <- T1[uncens] < T2[uncens]
  c2 <- T2[uncens] < T1[uncens]
  
  # tie handling
  tie <- !(c1 | c2)
  
  tmp <- integer(sum(uncens))
  tmp[c1] <- 1L # integer assiged
  tmp[c2] <- 2L
  if (any(tie)) tmp[tie] <- sample(c(1L, 2L), sum(tie), replace = TRUE)
  
  # final cause
  cause[uncens] <- tmp
  status <- as.integer(cause != 0) # just death/no death
  
  # out
  colnames(x1) <- sprintf("x1_%d", seq_len(p_block))
  colnames(x2) <- sprintf("x2_%d", seq_len(p_block))
  colnames(x3) <- sprintf("x3_%d", seq_len(p_block))
  
  #browser()
  data.frame(Tobs, status, cause, x1, x2, x3)
  
}

stratification <- function(time, event, folds = 5) {
  # stratification by event
  df$event <- factor(df$event)
  # sample
  rsample::vfold_cv(df, v = folds, strata = "event")
  
}

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
                        #times=seq(0,time.max,0.1), # switch for: tmp$time ?
                        times = sort(unique(time)), # change to this because of exploding times
                        level.chaos=1,
                        mode="matrix",
                        type="surv")
  
  out <- cbind(sort(unique(time)), prob)
  
  out <- na.omit(out)
  
  return(out)
}

## Calculate the censoring probabilities
## Using prodlim
## Controlling for competing risks
## No covariates in this function
censor.prob.KM.individual <- function(time, status, cens.code, predictions){
  tmp <- data.frame(time=time, status=status, predictions = predictions)
  # Sets the value in status equal to cens.code to 0, the rest to 1
  # Controlling for competing risks in censoring status
  tmp$censor.status <- ifelse(status == cens.code, 0, 1) # Not handled in riskRegression with cmprsks
  # create id 
  rownames(tmp)<- 1:length(time)
  time.max = ceiling(max(time))
  # We need to consider that the prodlim needs ordering 
  tmp <- tmp[order(tmp$time),]
  # Fit prodlim, reversed non-parametric survival KM
  # Switches the censoring status to estimate censoring distribution
  fit = prodlim::prodlim(formula=Surv(time,censor.status)~1,
                         data=tmp,
                         reverse=TRUE)
  # Predict weights at subject specific times
  ipcw.subject.times = prodlim::predictSurvIndividual(fit,lag=0)
  
  tmp$ipcw.subject.times <- ipcw.subject.times
  
  return(tmp)
}

# Effective computation of trapezoidal rule for numerical integration 
# of brier scores across time points of interes
trapezoidal.integration <- function(x, y) {
  if (length(x) != length(y)) stop("x and y must have the same length")
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}
