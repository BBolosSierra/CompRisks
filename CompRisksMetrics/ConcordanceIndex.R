#'
#' Title: Concordance Index for competing risks
#' 
#' Discrimination error. Time-dependent Concordance Index
#' This measures the proportion of correctly ordered risk pairs for the event k, based 
#' on the predicted risk of the event up to time tau
#' 
#' \deq{C_k(\tau) = \frac{\sum_{i=1}^N \sum_{j=1}^N (A_{ij} + B_{ij}) \cdot Q_{ij} \cdot N_i^k(\tau)}{\sum_{i=1}^N \sum_{j=1}^N (A_{ij} + B_{ij}) \cdot N_i^k(\tau)}}
#' 
#' A == risk ordering of patients, small time means patient 'i' at higher risk than patient 'j' experiencing event of interest 
#' A[i,j] = 0 for tied event times.
#' 
#' B == risk ordering of patients, large time for patient 'i' means lower risk than patient 'j' if not experienced the event of interest.
#' Ties are included in B 
#' 
#' Q == the risk ordering of the subjects, i.e., is subject i assigned a higher risk by the model than the subject j, for event Ek until time t.
#' Q[i,j] = 0 for tied predictions.
#' 
#' N_t == number of subjects with survival time < time point and experience event of interest
#' Tied event times are included
#' 
#' Ref: https://arxiv.org/pdf/1810.11207v1
#' 
#' @param predictions: vector of model predictions.
#' @param time time vector of times. Time-to-event.
#' @param status vector of events. 1) If survival -> [0,...,1]. 2) If competing risks -> [0,...,K].
#' @param tau evaluation time of interest.
#' @param cause event of interest.
#' @param cens.code value of censoring status, commonly censor patients have status of 0.
#' @param method 'survival' if the predictions are survival probabilitites or 'cifs' if they are cumulative incidence functions
#' 
#' @return numerical value, concordance index value as a measure of discrimination error.
#' 
#' 

CIndexCRisks <- function(predictions, 
                                      time,
                                      cens.code, 
                                      status, 
                                      cause, 
                                      tau,
                                      method=c("survival","cifs")){
  
  method = match.arg(method)
  
  censor.status <- ifelse(status == cens.code, 0, 1)
  
  if( method=="survival" ){ predictions=1-predictions; }
  
  n = length(predictions)
  A = matrix(0, nrow=n, ncol=n)
  B = matrix(0, nrow=n, ncol=n)
  Q = matrix(0, nrow=n, ncol=n) 
  N_t = matrix(0, nrow=n, ncol=n)
  Num_mat = matrix(0, nrow=n, ncol=n)
  Den_mat = matrix(0, nrow=n, ncol=n)
  
  Num=0
  Den=0
  for (i in  1:n){
    #print(n)
    A[i,which(time[i] < time)] = 1 
    B[i, intersect(intersect(which((time[i] >= time)),
                             which(status!=cause)), which(censor.status==1))] = 1
    Q[i,which(predictions[i]>predictions)]=1
  }
  
  for (i in 1:n){
    if(time[i]<=tau && status[i]==cause && censor.status[i]==1){
      N_t[i,] = 1
    }
  }
  
  Num_mat = (A+B)*Q*N_t
  Den_mat = (A+B)*N_t
  
  Num = sum(Num_mat)
  Den = sum(Den_mat)
  
  return(Num/Den)
  
}