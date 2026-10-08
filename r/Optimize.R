library(tidyverse)
source("R/KeyFunctions.R")

## Approach 1: Set Minimum Overall Power, then Minimize Expected (Or Quantile) Budget

OptimalBudget <- function(Delta, pi_exp,
                          p0, p1,
                          C1, C2, C3,
                          alpha=0.05,
                          NumCheck=20,
                          quantile=NULL) {
### Parameters:
  # Delta = minimum power (0-1)
  # pi_exp = expected exposure rate
  # p0, p1 = outcome rate in unexposed (p0) and exposed (p1)
  # C1 = cost from enrollment to exposure assessment
  # C2 = additional cost for exposed individuals
  # C3 = cost from exposure to assessment (for unexposed and exposed)
  # NumCheck = parameter passed to CheckSet.P; if 0, returns value for fixed N1
  #            o.w., returns value accounting for N1 variability
  #            higher values take slightly longer but may find slightly better result
  # quantile = if NULL, optimize expected budget. 
  #            o.w., optimize that quantile of the budget distribution (0-1)
  
  ## Find initial optima
  V.init <- Var.FromPower(p0,p1,Delta)
  N.init <- N.star.V(C1,C2,C3,p0,p1,pi_exp,V.init)
  
  Combos <- tibble(N=rep(floor(N.init),8)+c(rep(0,4),rep(1,4)),
                   L0=rep(c(0,1), times=4),
                   L1=rep(c(0,0,1,1), times=2)) %>%
    dplyr::mutate(N0=floor(p0*(1-p0)*N*pi_exp/(N*pi_exp*V.init-p1*(1-p1)))+L0,
                  N1=floor(N*pi_exp)+L1,
                  N1.E=N*pi_exp) %>%
    dplyr::select(-c(L0,L1))
  Combos$B.N1 <- apply(Combos, 1, FUN=function(row) B.N1(C1,C2,C3,row["N0"],row["N1"],row["N"]))
  Combos$SD.N1 <- apply(Combos, 1, FUN=function(row) sqrt(Var.N1(p0,p1,row["N0"],row["N1"])))
  Combos$Power.N1 <- apply(Combos, 1, FUN=function(row) Power.N1(p0,p1,row["N0"],row["N1"],alpha))
  Combos$B.E <- apply(Combos, 1, FUN=function(row) B.E(C1,C2,C3,row["N0"],pi_exp,row["N"]))
  Combos$SD.E <- apply(Combos, 1, FUN=function(row) sqrt(Var.E(p0,p1,row["N0"],pi_exp,row["N"])))
  Combos$Power.E <- apply(Combos, 1, FUN=function(row) Power.E(p0,p1,row["N0"],pi_exp,row["N"],alpha))

  if (max(Combos$Power.E) >= Delta) {
    ComboBest <- Combos %>% dplyr::filter(Power.E >= Delta, N0 > 0, N1 > 0, N > (N0+N1))
  } else if (max(Combos$Power.N1) >= Delta) {
    ComboBest <- Combos %>% dplyr::filter(Power.N1 >= Delta, N0 > 0, N1 > 0, N > (N0+N1))
  } else {
    ComboBest <- Combos %>% dplyr::filter(Power.E==max(Combos$Power.E), N0 > 0, N1 > 0, N > (N0+N1))
  }

  if (is.null(quantile)) {
    ComboBest <- ComboBest %>% dplyr::filter(B.E==min(ComboBest$B.E))
  } else {
    ComboBest$B.quantile <- apply(ComboBest, 1, FUN=function(row) B.quant(C1,C2,C3,row["N0"],pi_exp,row["N"],quantile))
    ComboBest <- ComboBest %>% dplyr::filter(B.quantile==min(ComboBest$B.quantile))
  }


  ## Find dynamic optima if requested

  if (NumCheck==0) {
    ComboBest$B.75 <- apply(ComboBest, 1, FUN=function(row) B.quant(C1,C2,C3,row["N0"],pi_exp,row["N"],0.75))
    ComboBest$B.90 <- apply(ComboBest, 1, FUN=function(row) B.quant(C1,C2,C3,row["N0"],pi_exp,row["N"],0.90))
    return(ComboBest %>%
             dplyr::select(N,N0,ends_with("quantile"),ends_with("E"),B.75,B.90,ends_with("N1")))
  } else {
    return(CheckSet.P(ComboBest,
                      Delta,pi_exp,
                      p0,p1,
                      C1,C2,C3,
                      alpha,
                      NumCheck,
                      quantile) %>% 
             dplyr::select(N,N0,ends_with("quantile"),ends_with("E"),B.75,B.90,ends_with("N1")))
  }
}

OptimalBudget(Delta=0.8,
              pi_exp=0.05,
              p0=0.1, p1=0.25,
              C1=500, C2=200, C3=1000,
              alpha=0.05,
              NumCheck=0,
              quantile=NULL)

OptimalBudget(Delta=0.8,
              pi_exp=0.05,
              p0=0.1, p1=0.25,
              C1=500, C2=200, C3=1000,
              alpha=0.05,
              NumCheck=20,
              quantile=NULL)


OptimalBudget(Delta=0.8,
              pi_exp=0.05,
              p0=0.1, p1=0.25,
              C1=500, C2=200, C3=1000,
              alpha=0.05,
              NumCheck=20,
              quantile=0.75)


## Approach 2: Set Maximum Expected (Or Quantile) Budget, then Minimize Overall Power

OptimalPower <- function(B, pi_exp,
                         p0, p1,
                         C1, C2, C3,
                         alpha=0.05,
                         NumCheck=5) {
  ### Parameters:
  # B = maximum budget (in units of C1,C2,C3)
  # pi_exp = expected exposure rate
  # p0, p1 = outcome rate in unexposed (p0) and exposed (p1)
  # C1 = cost from enrollment to exposure assessment
  # C2 = additional cost for exposed individuals
  # C3 = cost from exposure to assessment (for unexposed and exposed)
  # alpha = significance level (0-1)
  # NumCheck = parameter passed to CheckSet.B; the value does not need to be very high in general

  ## Find initial optima
  
  N.init.B <- N.star.B(C1,C2,C3,p0,p1,pi_exp,B)

  Combos <- tibble(N=rep(floor(N.init.B),8)+c(rep(0,4),rep(1,4)),
                   L0=rep(c(0,1), times=4),
                   L1=rep(c(0,0,1,1), times=2)) %>%
    dplyr::mutate(N0=floor((B-N*(C1+pi_exp*C2+pi_exp*C3))/C3)+L0,
                  N1=floor(N*pi_exp)+L1) %>%
    dplyr::select(-c(L0,L1))
  Combos$B.E <- apply(Combos, 1, FUN=function(row) B.E(C1,C2,C3,row["N0"],pi_exp,row["N"]))
  Combos$Power.E <- apply(Combos, 1, FUN=function(row) Power.E(p0,p1,row["N0"],pi_exp,row["N"],alpha))
  Combos$SD.E <- apply(Combos, 1, FUN=function(row) sqrt(Var.E(p0,p1,row["N0"],pi_exp,row["N"])))

  ## Find dynamic optima

  if (min(Combos$B.E) <= B) {
    ComboBest <- Combos %>% dplyr::filter(B.E <= B, N0 > 0, N1 > 0, N > (N0+N1))
    ComboBest <- ComboBest %>% dplyr::filter(Power.E==max(ComboBest$Power.E)) %>%
      dplyr::select(-ends_with("N1")) %>%
      distinct()
    ComboBest$N1.E <- ComboBest$N*pi_exp
    ComboBest$B.75 <- apply(ComboBest, 1, FUN=function(row) B.quant(C1,C2,C3,row["N0"],pi_exp,row["N"],quantile=0.75))
    ComboBest$B.90 <- apply(ComboBest, 1, FUN=function(row) B.quant(C1,C2,C3,row["N0"],pi_exp,row["N"],quantile=0.90))
    ComboBest$N1 <- round(ComboBest$N1.E)
    ComboBest$B.N1 <- apply(ComboBest, 1, FUN=function(row) B.N1(C1,C2,C3,row["N0"],row["N1"],row["N"]))
    ComboBest$SD.N1 <- apply(ComboBest, 1, FUN=function(row) sqrt(Var.N1(p0,p1,row["N0"],row["N1"])))
    ComboBest$Power.N1 <- apply(ComboBest, 1, FUN=function(row) Power.N1(p0,p1,row["N0"],row["N1"],alpha))
    return(ComboBest)
  } else if (NumCheck==0) {
    stop("No combination found that meets required budget. Try increasing NumCheck.")
  } else {
    return(CheckSet.B(Starts=Combos[1,],
                      B,pi_exp,
                      p0,p1,
                      C1,C2,C3,
                      alpha,
                      NumCheck) %>%
             dplyr::select(N,N0,ends_with(".E"),B.75,B.90,ends_with("N1")))
  }
}

OptimalPower(B=1080000, pi_exp=0.05,
             p0=0.1, p1=0.25,
             C1=500, C2=200, C3=1000,
             alpha=0.05,
             NumCheck=5)
