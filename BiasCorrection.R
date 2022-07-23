# 
# R functions for correcting prevalence estimation for biased 
# sampling with testing errors.
#
# README:
#
# Function is used to correct the prevalence estimation with both
# sampling bias and testing errors, assuming that error rates
# are stratified by syptoms, which are symptomatic individuals
# and asymptomatic individuals. Different groups have their own
# error rates. Besides, the function also assesses the result of
# correction by giving active information for each example.

prev.correct.strata <- function(tildeN, N_T, alpha_0, alpha_1, 
                                beta_0, beta_1){
#
# Arguments:
#  
#  tildeN: the matrix of observed disease status by category of 
#          sysptoms in population
#     N_T: the matrix of observed disease status by category of 
#          sysptoms for each sample
# alpha_0: false positive rate for asymptomatic individuals
# alpha_1: false positive rate for symptomatic individuals
#  beta_0: false negative rate for asymptomatic individuals
#  beta_1: false negative rate for symptomatic individuals
#
# Result components:
#
#       p: real prevalence in population
# tilde_p: naive estimator of prevalence
#   hat_p: corrected estimator of prevalence
#      AI: active information for assessing the correction
  
  tildep_0_0 <- N_T[2,2]/sum(N_T)
  tildep_0_1 <- N_T[2,1]/sum(N_T)
  tildep_1_0 <- N_T[1,2]/sum(N_T)
  tildep_1_1 <- N_T[1,1]/sum(N_T)
  tildep_1 <- tildep_1_0 + tildep_1_1
  tildep_0 <- tildep_0_0 + tildep_0_1
  
  
  barp_T_0_0 <- (tildep_0_0 - beta_0 * tildep_0)/(1 - alpha_0 - beta_0)
  barp_T_1_0 <- (tildep_1_0 - beta_1 * tildep_1)/(1 - alpha_1 - beta_1)
  barp_T_0_1 <- (alpha_0 * tildep_0 - tildep_0_1)/(alpha_0 + beta_0 - 1)
  barp_T_1_1 <- (alpha_1 * tildep_1 - tildep_1_1)/(alpha_1 + beta_1 - 1)
  
  barp_T_1_star <- barp_T_1_0 + barp_T_1_1
  
  if (sum(tildeN[1,]) == sum(N_T[1,])) {
    hatp_1 <- barp_T_1_star/2 * (sum(N_T)/sum(tildeN))
  }else{
    hatp_1 <- barp_T_1_star/2 * (sum(N_T)/sum(tildeN) + 1)
  }
  
  hatp_1_1 <- hatp_1 * N_T[1,1] / sum(N_T[1,])
  
  barp_T_0_star <- barp_T_0_0 + barp_T_0_1
  hatp_0_1 <- barp_T_0_1/barp_T_0_star * (1 - hatp_1)
  hatp <- hatp_0_1 + hatp_1_1
  
  tildep <- sum(N_T[,1])/sum(N_T)
  AI <- log(hatp/p)
  cat("p =", p, "\n")
  cat("tilde_p = ", tildep, "\n")
  cat("hat_p = ", hatp, "\n")
  cat("AI = ", AI, "\n")
}

#
# Paper Sampling Protocols:
#
alpha_0 = 0.001
alpha_1 = 0.005
beta_0 = 0.1
beta_1 = 0.05
tildeN <- matrix(c(1754, 6639, 108, 90731), 2)
# 75% symptomatic and 25% asymptomatic, all symptomatics get sampled.
N_T_1 <- matrix(c(1754, 100, 108, 521), 2) 
# 75% symptomatic and 25% asymptomatic
N_T_2 <- matrix(c(141, 8, 9, 42), 2) 
# 50% symptomatic and 50% asymptomatic
N_T_3 <- matrix(c(94, 16, 6, 84), 2)
# totally random
N_T_4 <- matrix(c(4, 32, 1, 163), 2) 
p <- (round(tildeN[1,1] - tildeN[1,1] * alpha_1 + tildeN[1,2] * beta_1) +
        round(tildeN[2,1] - tildeN[2,1] * alpha_0 + tildeN[2,2] * beta_0)) /
        sum(tildeN)

prev.correct.strata(tildeN, N_T_1, alpha_0, alpha_1, beta_0, beta_1)
prev.correct.strata(tildeN, N_T_2, alpha_0, alpha_1, beta_0, beta_1)
prev.correct.strata(tildeN, N_T_3, alpha_0, alpha_1, beta_0, beta_1)
prev.correct.strata(tildeN, N_T_4, alpha_0, alpha_1, beta_0, beta_1)
