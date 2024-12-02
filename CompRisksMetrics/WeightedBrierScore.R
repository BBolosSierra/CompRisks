#'
#' Title: Weighted Brier score
#' 
#' Ref: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-022-01679-6
#' Ref: https://square.github.io/pysurvival/metrics/brier_score.html
#'
#' Weighted overall prediction error, is a measure of the mean squared error between predicted probabilities and actual outcomes, 
#' adjusted for censoring using Inverse Probability of Censoring Weights (IPCW) 
#' Numerical value from 0 to 1 where lower values indicate better prediction performance.
#' 
#' \deq{BS = \frac{1}{N} \sum_{i=1}^{N}W_{i}(\tau)(\hat{y}_{i}(\tau)-y_{i}(\tau))^2}
#' 
#' Where N is number of observations, and W are the Inverse Probability Censoring Weights (IPCW).
#' \deq{\hat{y}_{i}(\tau)} is the predicted probability of the event at time for individual i by the model. 
#' \deq{y_{i}(\tau)} is the event outcome is the actual event outcome (1 if the event occurred, 0 if it did not).
#' 
#' The weights per individual are calculated in the following way:
#' \deq{W_{i} = \frac{\mathbb{I}(T_{i} \leq \tau, \delta_{i} = 1)}{\hat{G}(T_{i}^{-})} + \frac{\mathbb{I}(T_{i} > \tau)}{\hat{G}(\tau)}}
#' \deq{\hat G(\tau) } is the censoring probabilities calculated withe the Kaplan-Meier survival probability of being uncensored at time tau. 
#' \deq{\hat G(T_{i}^{-})} is the censoring probabilities right just before time T
#' \deq{\delta_{i}} is an indicator function. Value is 1 when any event occurs, 0 otherwise.
#'
#' Function censor.prob.KM to calculate the IPCW.
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
#' 
#' 
#' @return 
#'  \itemize{
#'   \item brier.score: Unweighted brier score
#'   \item weighted.brier.score: Numeric value, weighted overall prediction error ranging from 0 to 1, where lower scores indicate better performance at specific time tau of interest.
#'   \item tau: vector of times of evaluation
#'   \item n: Number of observations
#'   \item n.risk: Number of observation at risk at time tau.
#'   
#' }
#' 
#' @examples
#' # example code
#' 
#' 

# Same as before but consider censoring (not competing risks)
#WeightedBinaryBrierScore <- function(predictions, time, status, tau, cens.code) {
#  # Calculate individual censoring weights 
#  ind.cens.weights = censor.prob.KM.individual(time = time, 
#                                               status = status, 
#                                               cens.code = cens.code)
#  
#  # Identify the binary outcome for the event of interest
#  y_true = ((time <= tau) & (status != cens.code))
#  
#  # Initialize vector for residuals
#  residuals = numeric(length(predictions))
#  
#  ## Consider censoring without competing risks
#  # Case 1: Event occurs before tau
#  index = (time <= tau & status != cens.code)
#  residuals[index] = ((1 - predictions[index])^2) / ind.cens.weights[index]
#  
#  # Case 2: Censored before tau (including competing risks treated as censored)
#  index = (time <= tau & status == cens.code)
#  residuals[index] = 0
#  
#  # Case 3: Event or censoring after tau
#  index = (time > tau)
#  residuals[index] = (predictions[index]^2) / ind.cens.weights[index]
#  
#  # Compute the Brier Score
#  BrierScoreCensored = mean(residuals)
#  
#  return(BrierScoreCensored)
#}

# Same as before but consider censoring (not competing risks)
WeightedBinaryBrierScore <- function(predictions, time, status, tau, cause, cens.code) {
  # Calculate individual censoring weights 
  ind.cens.weights = censor.prob.KM.individual(time = time, 
                                               status = status, 
                                               cens.code = cens.code)
  
  # Identify the binary outcome for the event of interest
  y_true = ((time <= tau) & (status == cause))
  
  # Initialize vector for residuals
  residuals = numeric(length(predictions))
  
  ## Consider censoring without competing risks
  # Case 1: Event occurs before tau
  index = (time <= tau & status == cause)
  residuals[index] = ((1 - predictions[index])^2) / ind.cens.weights[index]
  
  # Case 2: Censored before tau (including competing risks treated as censored)
  index = (time <= tau & status == cens.code)
  residuals[index] = 0
  
  # Case 3: Event or censoring after tau
  index = (time > tau)
  residuals[index] = (predictions[index]^2) / ind.cens.weights[index]
  
  # Compute the Brier Score
  BrierScore = mean(residuals)
  
  return(BrierScore)
}


## Weighted BRIER-SCORE
WeightedBinaryBrierScore2 <- function(predictions,
                                       tau,
                                       time, 
                                       status,
                                       cause,
                                       cens.code){
  # Extract the censoring probabilities
  G <- censor.prob.KM(time=time, 
                         status=status, 
                         cens.code=cens.code)
  # G[,1] Time points where censoring probabilities are estimated
  # G[,2] Estimated censoring probabilities
  
  # Create a empty vector for weights
  n     = length(predictions)
  W     = rep(0,n)
  
  # Binary indicator for individuals that are still at risk
  Y_tilde = (time > tau)
  # Binary indicator for individuals that have experienced the event at time tau
  Y_true  = (time <= tau) * (status == cause)
  
  # Calculate G1 and G2
  # Compute IPCW
  for( i in 1:n ){
    
    # Get the index of censoring time greater or equal than the individual observed times
    indx1 = which(G[,1] >= time[i])
    # Get the index of time points greater or equal than evaluation time
    indx2 = which(G[,1] >= tau)
    
    # Calculate G1: censoring probabilities just before observed times T 
    if( length(indx1) > 0 ){
      G1 = G[indx1[1],2] # Takes the earliest estimated censoring probabilities 
    } else {
      G1 = 1
    }
    # Calculate G2: censoring probabilities at tau or greater
    if( length(indx2) > 0 ){
      G2 = G[indx2[1],2]
    } else {
      G2 = 1
    }
    
    W[i] = Y_tilde[i]/G2 + Y_true[i]/G1        

  }    
  
  S    = (1/n)*sum(Y_tilde)

  norm        = 1/n
  BS_we       = sum(W*(Y_true - predictions)^2)*norm
  
  return(list(weighted.brier.score = BS_we, 
              tau = tau,
              n=n, 
              n.risk=S*n))
}


#### Competing risks

WeightedCRBrierScore <- function(predictions, time, status, tau, cause, cens.code) {
  # Calculate individual censoring weights
  ind.cens.weights = censor.prob.KM.individual(time=time, 
                                               status=status, 
                                               cens.code=cens.code)
  
  # Identify the cause of interest
  y_true = ((time <= tau) & (status == cause))
  
  # Initialize vecotr for residuals
  residuals = numeric(length(predictions))
  
  ## As in riskRegression, consider competing risks and censoring
  # Case 1: Event of interest occurs before tau
  index = (time <= tau & status == cause)
  residuals[index] = ((1 - predictions[index])^2) / ind.cens.weights[index]
  
  # Case 2: Competing event occurs before tau
  index = (time <= tau & status != cause & status != 0)
  residuals[index] = (predictions[index]^2) / ind.cens.weights[index]
  
  # Case 3: Censored before tau
  index = (time <= tau & status == 0)
  residuals[index] = 0
  
  # Case 4: Event or censoring after tau
  index = (time > tau)
  residuals[index] = (predictions[index]^2) / ind.cens.weights[index]
  
  # Compute the Brier Score
  BrierScoreCompRisks = mean(residuals)
  
  return(BrierScoreCompRisks)
}

WeightedCRBrierScore2 <- function(predictions, 
                                 tau, 
                                 time, 
                                 status, 
                                 cause, 
                                 cens.code) {
  # Extract the censoring probabilities
  G <- censor.prob.KM(time = time, status = status, cens.code = cens.code)
  
  # Create an empty vector for weights
  n <- length(predictions)
  W <- rep(0, n)
  
  # Binary indicator for individuals that are still at risk at tau
  Y_tilde <- (time > tau)
  
  # Binary indicator for individuals that have experienced the event of interest by tau
  Y_true <- (time <= tau & status == cause)
  
  # Create a residuals vector for the weighted Brier score
  residuals <- rep(0, n)
  
  # Calculate G1 and G2
  for (i in 1:n) {
    # Get the index of censoring time greater or equal than the individual's observed time
    indx1 <- which(G[, 1] >= time[i])
    # Get the index of time points greater or equal than evaluation time
    indx2 <- which(G[, 1] >= tau)
    
    # Calculate G1: censoring probabilities just before observed times T
    G1 <- if (length(indx1) > 0) G[indx1[1], 2] else 1
    
    # Calculate G2: censoring probabilities at tau or greater
    G2 <- if (length(indx2) > 0) G[indx2[1], 2] else 1
    
    # Handle weights based on the event type
    if (status[i] == cause && time[i] <= tau) {
      # Case 1: Event of interest occurs before tau
      W[i] <- 1 / G1
      residuals[i] <- (1 - predictions[i])^2 / G1
    } else if (status[i] != cause && status[i] != cens.code && time[i] <= tau) {
      # Case 2: Competing risks before tau
      W[i] <- 1 / G1
      residuals[i] <- (predictions[i]^2) / G1
    } else if (status[i] == cens.code && time[i] <= tau) {
      # Case 3: Censored before tau
      W[i] <- 0
      residuals[i] <- 0
    } else if (time[i] > tau) {
      # Case 4: Event or censoring after tau
      W[i] <- 1 / G2
      residuals[i] <- (predictions[i]^2) / G2
    }
  }
  
  # Normalize residuals to compute the weighted Brier score
  BS_we <- mean(residuals)
  
  return(list(weighted.brier.score = BS_we,
              tau = tau,
              n = n,
              n.risk = sum(Y_tilde)))
}
