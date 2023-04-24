# 
# R functions for correcting prevalence estimation for biased 
# sampling with testing errors.
# 
# Author: Lili Zhou
#
# Date: 04/24/2023
#
# README:
#
# Function is used to correct the prevalence estimation with both
# sampling bias and testing errors, assuming that error rates
# are stratified by syptoms, which are symptomatic individuals
# and asymptomatic individuals. Different groups have their own
# error rates. 

prev.correct.strata <- function(N, N_T, alpha_0, alpha_1, beta_0, 
                                beta_1, test.err = TRUE, symp.all = FALSE){
  #
  # Arguments:
  #  
  #        N: population size
  #      N_T: the matrix of observed sample with different 
  #           symtoms status and testing results
  #  alpha_0: false positive rate for asymptomatic individuals
  #  alpha_1: false positive rate for symptomatic individuals
  #   beta_0: false negative rate for asymptomatic individuals
  #   beta_1: false negative rate for symptomatic individuals
  # test.err: indicator of if testing errors are considered
  #           in the study
  # symp.all: indicator of if all symptomatic individuals are
  #           sampled from population in the study
  
  #
  # Result components:
  #
  # tilde_p: naive estimator of prevalence
  #   hat_p: corrected estimator of prevalence

  
  if (!test.err) {
    p <- sum(N_T[,1])/sum(N_T)
    
    hatp_T_0_0 <- N_T[2,2]/sum(N_T)
    hatp_T_0_1 <- N_T[2,1]/sum(N_T)
    hatp_T_1_0 <- N_T[1,2]/sum(N_T)
    hatp_T_1_1 <- N_T[1,1]/sum(N_T) 
    hatp_T_1 <- hatp_T_1_0 + hatp_T_1_1
    hatp_T_0 <- hatp_T_0_0 + hatp_T_0_1
    
    if (symp.all) { 
      p_1 <- hatp_T_1 * (sum(N_T)/N)
      hatp_1_1 <- p_1 * (hatp_T_1_1/hatp_T_1)
      hatp_0_1 <- (hatp_T_0_1/hatp_T_0) * (1 - p_1)
      hatp <- hatp_0_1 + hatp_1_1
    } else {
      hatp_1 <- hatp_T_1/2 * (sum(N_T)/N + 1)
      hatp_1_1 <- hatp_1 * (hatp_T_1_1/hatp_T_1)
      hatp_0_1 <- (hatp_T_0_1/hatp_T_0) * (1 - hatp_1)
      hatp <- hatp_0_1 + hatp_1_1
    }
    
    tildep <- hatp_T_0_1 + hatp_T_1_1
    cat("tilde_p = ", tildep, "\n")
    cat("hat_p = ", hatp, "\n")

    
  } else {
    
    tildep <- sum(N_T[,1])/sum(N_T)
    
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
    
    barp_T_1 <- barp_T_1_0 + barp_T_1_1
    barp_T_0 <- barp_T_0_0 + barp_T_0_1
    
    if (symp.all) {
      hatp_1 <- barp_T_1 * (sum(N_T)/N)
    }else{
      hatp_1 <- barp_T_1/2 * (sum(N_T)/N + 1)
    }
    
    hatp_1_1 <- hatp_1 * (barp_T_1_1/barp_T_1)
    hatp_0_1 <- barp_T_0_1/barp_T_0 * (1 - hatp_1)
    hatp <- hatp_0_1 + hatp_1_1
    
    cat("tilde_p = ", tildep, "\n")
    cat("hat_p = ", hatp, "\n")

  }
}

