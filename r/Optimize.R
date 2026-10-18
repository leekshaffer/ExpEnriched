library(tidyverse)
source("R/KeyFunctions.R")

## Approach 1: Set Minimum Overall Power, then Minimize Expected (Or Quantile) Budget

OptimalBudget <- function(Delta,
                          pi_0,pi_1,
                          rho_1,pi_e,rho_2,
                          C_1,C_2,C_e,
                          alpha=0.05,
                          NumCheck=20,
                          quantile=NULL) {
### Parameters:
  # Delta = minimum power (0-1)
  # pi_0 = outcome probability in unexposed individuals
  # pi_1 = outcome probability in exposed individuals
  # rho_1 = retention probability (1-loss to follow-up probability) 
  #         between enrollment and exposure determination
  # pi_e = probability of exposure
  # rho_1 = retention probability (1-loss to follow-up probability) 
  #         between exposure determination and study end
  # C_1 = cost from enrollment to exposure assessment
  # C_2 = cost from exposure to assessment (for unexposed and exposed)
  # C_e = additional cost for exposed individuals
  # alpha = significance level (max. Type I Error Probability)
  # NumCheck = parameter passed to CheckSet.P; if 0, returns approximate value (ignoring variability in N1)
  #            o.w., returns value accounting for additional sources of variability
  #            higher values take slightly longer but may find slightly better result
  # quantile = if NULL, optimize expected budget. 
  #            o.w., optimize that quantile of the budget distribution (0-1)
  
  ## Find initial optima
  V.init <- DeltaToV(pi_0,pi_1,Delta,alpha)
  N.init <- N.star.V(pi_e, pi_0, pi_1, C_1, C_2, C_e, rho_1, rho_2, V.init)
  
  Combos <- tibble(N=rep(floor(N.init),8)+c(rep(0,4),rep(1,4)),
                   L0=rep(c(0,1), times=4),
                   L1=rep(c(0,0,1,1), times=2)) %>%
    dplyr::mutate(N0=floor(NVtoN0(pi_0,pi_1,rho_1,pi_e,rho_2,N,V.init))+L0,
                  N1=floor(N*pi_e)+L1,
                  N1.E=N*pi_e) %>%
    dplyr::select(-c(L0,L1))
  
  Combos$C.N1 <- apply(Combos, 1, FUN=function(row) 
    C.N0N1N(C_1,C_2,C_e,row["N0"], row["N1"], row["N"]))
  Combos$SD.N1 <- apply(Combos, 1, FUN=function(row) 
    sqrt(V.N0N1(pi_0,pi_1,row["N0"],row["N1"],rho_2)))
  Combos$Power.N1 <- apply(Combos, 1, FUN=function(row) 
    VToDelta(pi_0,pi_1,row["SD.N1"]^2,alpha))
  Combos$C.val <- apply(Combos, 1, FUN=function(row) 
    C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile))
  Combos$C.E <- apply(Combos, 1, FUN=function(row) 
    C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=NULL))
  Combos$SD.E <- apply(Combos, 1, FUN=function(row) 
    sqrt(V.N0N(pi_0,pi_1,row["N0"],row["N"],rho_1,pi_e,rho_2)))
  Combos$Power.E <- apply(Combos, 1, FUN=function(row) 
    VToDelta(pi_0,pi_1,row["SD.E"]^2,alpha=0.05))

  if (max(Combos$Power.E) >= Delta) {
    ComboBest <- Combos %>% dplyr::filter(Power.E >= Delta, N0 > 0, N1 > 0, N > (N0+N1))
  } else if (max(Combos$Power.N1) >= Delta) {
    ComboBest <- Combos %>% dplyr::filter(Power.N1 >= Delta, N0 > 0, N1 > 0, N > (N0+N1))
  } else {
    ComboBest <- Combos %>% dplyr::filter(Power.E==max(Combos$Power.E), N0 > 0, N1 > 0, N > (N0+N1))
  }

  ComboBest <- ComboBest %>% dplyr::filter(C.val==min(ComboBest$C.val),
                                           N1==round(N1.E))


  ## Find dynamic optima if requested

  if (NumCheck==0) {
    ComboBest$C.75 <- apply(ComboBest, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=0.75))
    ComboBest$C.90 <- apply(ComboBest, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=0.90))
    return(ComboBest %>%
             dplyr::select(N,N0,Power.E,C.val,
                           C.E,SD.E,
                           C.75,C.90,
                           N1.E,ends_with("N1")))
  } else {
    ComboBest <- ComboBest %>%
      dplyr::select(N,N0) %>%
      dplyr::distinct()
    if (nrow(ComboBest) > 1) {
      Starts <- tibble(N=min(ComboBest$N), N0=min(ComboBest$N0))
      NumCheck <- NumCheck+max(max(ComboBest$N)-min(ComboBest$N),
                               max(ComboBest$N0)-min(ComboBest$N0))
    }
    return(CheckSet.P(Starts,
                      Delta,
                      pi_0,pi_1,
                      rho_1,pi_e,rho_2,
                      C_1,C_2,C_e,
                      alpha=0.05,
                      NumCheck=NumCheck,
                      quantile=NULL))
  }
}

OptimalBudget(Delta=0.8,
              pi_0=0.1, pi_1=0.25,
              rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
              C_1=1000, C_2=500, C_e=200,
              alpha=0.05,
              NumCheck=0,
              quantile=NULL)

OptimalBudget(Delta=0.8,
              pi_0=0.1, pi_1=0.25,
              rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
              C_1=1000, C_2=500, C_e=200,
              alpha=0.05,
              NumCheck=20,
              quantile=NULL)


OptimalBudget(Delta=0.8,
              pi_0=0.1, pi_1=0.25,
              rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
              C_1=1000, C_2=500, C_e=200,
              alpha=0.05,
              NumCheck=20,
              quantile=0.75)


## Approach 2: Set Maximum Expected (Or Quantile) Budget, then Minimize Overall Power

OptimalPower <- function(B,
                         pi_0,pi_1,
                         rho_1,pi_e,rho_2,
                         C_1,C_2,C_e,
                         alpha=0.05,
                         NumCheck=5,
                         quantile=NULL) {
  ### Parameters:
  # B = maximum budget (in the same units as C_1, C_2, C_e)
  # pi_0 = outcome probability in unexposed individuals
  # pi_1 = outcome probability in exposed individuals
  # rho_1 = retention probability (1-loss to follow-up probability) 
  #         between enrollment and exposure determination
  # pi_e = probability of exposure
  # rho_1 = retention probability (1-loss to follow-up probability) 
  #         between exposure determination and study end
  # C_1 = cost from enrollment to exposure assessment
  # C_2 = cost from exposure to assessment (for unexposed and exposed)
  # C_e = additional cost for exposed individuals
  # alpha = significance level (max. Type I Error Probability)
  # NumCheck = parameter passed to CheckSet.B; if 0, returns approximate value (ignoring variability in N1)
  #            o.w., returns value accounting for additional sources of variability
  #            generally does not need to be very high
  # quantile = a quantile of the budget to report [note: not used for optimizing here, just reported]
  #            defaulte (NULL) gives expectation

  ## Find initial optima
  
  N.init.B <- N.star.B(pi_e, pi_0, pi_1, C_1, C_2, C_e, rho_1, B)

  Combos <- tibble(N=rep(floor(N.init.B),8)+c(rep(0,4),rep(1,4)),
                   L0=rep(c(0,1), times=4),
                   L1=rep(c(0,0,1,1), times=2)) %>%
    dplyr::mutate(N0=floor(NBtoN0(rho_1,pi_e,C_1,C_2,C_e,N.init.B,B))+L0,
                  N1=floor(N*pi_e)+L1) %>%
    dplyr::select(-c(L0,L1))
  
  Combos$C.val <- apply(Combos, 1, FUN=function(row) 
    C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile))
  Combos$C.E <- apply(Combos, 1, FUN=function(row) 
    C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=NULL))
  Combos$SD.E <- apply(Combos, 1, FUN=function(row) 
    sqrt(V.N0N(pi_0,pi_1,row["N0"],row["N"],rho_1,pi_e,rho_2)))
  Combos$Power.E <- apply(Combos, 1, FUN=function(row) 
    VToDelta(pi_0,pi_1,row["SD.E"]^2,alpha=0.05))

  ## Find dynamic optima

  if (min(Combos$C.E) <= B) {
    ComboBest <- Combos %>% dplyr::filter(C.E <= B, N0 > 0, N1 > 0, N > (N0+N1))
    ComboBest <- ComboBest %>% dplyr::filter(Power.E==max(ComboBest$Power.E)) %>%
      dplyr::select(-ends_with("N1")) %>%
      distinct()
    
    ComboBest$C.75 <- apply(ComboBest, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=0.75))
    ComboBest$C.90 <- apply(ComboBest, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=0.90))
    
    ComboBest <- ComboBest %>%
      dplyr::mutate(N1.E=N*pi_e,
                    N1=round(N1.E))
    
    ComboBest$C.N1 <- apply(ComboBest, 1, FUN=function(row) 
      C.N0N1N(C_1,C_2,C_e,row["N0"], row["N1"], row["N"]))
    ComboBest$SD.N1 <- apply(ComboBest, 1, FUN=function(row) 
      sqrt(V.N0N1(pi_0,pi_1,row["N0"],row["N1"],rho_2)))
    ComboBest$Power.N1 <- apply(ComboBest, 1, FUN=function(row) 
      VToDelta(pi_0,pi_1,row["SD.N1"]^2,alpha))
    return(ComboBest %>%
             dplyr::select(N,N0,Power.E,C.val,
                           C.E,SD.E,
                           C.75,C.90,
                           N1.E,ends_with("N1")))
  } else if (NumCheck==0) {
    stop("No combination found that meets required budget. Try increasing NumCheck.")
  } else {
    Res <- CheckSet.B(Starts=Combos[1,],
                      B,
                      pi_0,pi_1,
                      rho_1,pi_e,rho_2,
                      C_1,C_2,C_e,
                      alpha=0.05,
                      NumCheck=5)
    Res$C.val <- apply(Res, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile))
    return(Res %>%
             dplyr::select(N,N0,Power.E,C.val,
                           C.E,SD.E,
                           C.75,C.90,
                           N1.E,ends_with("N1")))
  }
}

OptimalPower(B=1900000,
             pi_0=0.1,pi_1=0.25,
             rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
             C_1=1000, C_2=500, C_e=200,
             alpha=0.05,
             NumCheck=5,
             quantile=NULL)
