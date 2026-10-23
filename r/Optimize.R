library(tidyverse)
source("r/KeyFunctions.R")

## Approach 1: Set Minimum Overall Power, then Minimize Expected Budget

OptimalBudget <- function(Delta,
                          theta, s_0, s_1,
                          pi_e, rho_1, rho_2,
                          C_1, C_2, C_e,
                          alpha=0.05,
                          restrict=5,
                          NumCheck=10,
                          nuCheck=2,
                          round=3) {
### Parameters:
  # Delta = minimum power (0-1)
  # theta = desired absolute effect size to power (abs(mu_1-mu_0))
  # s_0 = outcome standard deviation in unexposed individuals
  # s_1 = outcome standard deviation in exposed individuals
  # Note: to power to Cohen's d, set s_0=1, s_1=1, and theta=Cohen's d value
  # pi_e = probability of exposure
  # rho_1 = retention probability (1-loss to follow-up probability) 
  #         between enrollment and exposure determination
  # rho_2 = retention probability (1-loss to follow-up probability) 
  #         between exposure determination and study end
  # C_1 = cost from enrollment to exposure assessment
  # C_2 = cost from exposure to assessment (for unexposed and exposed)
  # C_e = additional cost for exposed individuals
  # alpha = significance level (max. Type I Error Probability)
  # restrict = method to speed up power calculation; higher value is more precise but slower (NULL is exact)
  # NumCheck = parameter passed to CheckSet.D to increment N;
  #            higher values take slightly longer but may find slightly better result
  # nuCheck = parameter passed to CheckSet.D to increment nu;
  #            higher values take slightly longer but may find slightly better result
  # if NumCheck and nuCheck are both set to 0, then the approximate optimum will be returned
  # round = if NULL, uses exact nu value. Otherwise rounds to this many decimal places
  
  ## Find initial optima
  Start <- N.star.V(s_0, s_1, pi_e, rho_1, rho_2, C_1, C_2, C_e,
                                V=DeltaToV(theta, Delta, alpha), round)
  if (NumCheck==0 & nuCheck==0) {
    Opt <- tibble(N=Start["N"],
                  nu=Start["nu"]) %>%
      dplyr::mutate(Power=Delta.N(theta, s_0, s_1, nu, pi_e, rho_1, rho_2, N, alpha, restrict),
                    SD=sqrt(V.N(s_0, s_1, nu, pi_e, rho_1, rho_2, N, restrict)),
                    Cost.E=C.N(C_1, C_2, C_e, nu, pi_e, rho_1, N, quantile=NULL),
                    Ns.E=N*nu,
                    M0.E=Ns.E*rho_1*(1-pi_e)*rho_2,
                    M1.E=N*rho_1*pi_e*rho_2)
  } else {
    Opt <- CheckSet.D(Start["N"],
                      Delta,
                      theta, s_0, s_1,
                      pi_e, rho_1, rho_2,
                      C_1, C_2, C_e,
                      alpha,
                      restrict,
                      NumCheck,
                      nuCheck,
                      round)
  }
  return(Opt)
}

## START HERE

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
      Delta.N0N1(pi_0,pi_1,row["N0"],row["N1"],rho_2,alpha))
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

## Interim Optimization:

### Fixed Budget

OptimalIntBudget <- function(Delta,
                             pi_0, pi_1,
                             rho_2, C_2,
                             N_0, N_1,
                             alpha=0.05,
                             NumCheck=5) {
  Power.Max <- Delta.N0N1(pi_0, pi_1, N_0, N_1, rho_2, alpha)
  if (Power.Max < Delta) {
    warning("The desired power cannot be achieved with these N values. The function will return the maximum power.")
    Checks <- tibble(N0=N_0, N1=N_1, Power.E=Power.Max) %>%
      dplyr::mutate(C.Remaining=(N_0+N_1)*C_2,
                    SD.E=sqrt(V.N0N1(pi_0, pi_1, N_0, N_1, rho_2)))
  } else {
    V_Est <- DeltaToV(pi_0,pi_1,Delta,alpha=0.05)
    if (N_1 <= N_0) {
      M1_Est <- N_1*rho_2
      M0_Est <- (pi_0*(1-pi_0))/(V_Est-pi_1*(1-pi_1)/M1_Est)
      N0_Start <- M0_Est/rho_2
      Checks <- tibble(N0=(floor(N0_Start)):min((floor(N0_Start)+NumCheck),N_0),
                       N1=N_1)
    } else {
      M0_Est <- N_0*rho_2
      M1_Est <- (pi_1*(1-pi_1))/(V_Est-pi_0*(1-pi_0)/M0_Est)
      N1_Start <- N1_Est/rho_2
      Checks <- tibble(N0=N_0,
                       N1=N_1(floor(N1_Start)):min((floor(N1_Start)+NumCheck),N_1))
    }
  
    Checks$Power.E <- apply(Checks, 1, function(row) 
      Delta.N0N1(pi_0, pi_1, row["N0"], row["N1"], rho_2, alpha))
    
    if (max(Checks$Power.E) < Delta) {
      message("No combination found that achieves desired power. Increasing NumCheck.")
      return(OptimalIntBudget(Delta, pi_0, pi_1, rho_2, C_2,
                         N_0, N_1, alpha,
                         2*(NumCheck+1)))
    } else {
      Checks <- Checks %>%
        dplyr::filter(Power.E >= Delta)
      Checks <- Checks %>%
        dplyr::filter(N0==min(Checks$N0),
                      N1==min(Checks$N1))
    }
    Checks$C.Remaining <- (Checks$N0+Checks$N1)*C_2
    Checks$SD.E <- apply(Checks, 1, function(row)
      sqrt(V.N0N1(pi_0,pi_1,row["N0"], row["N1"],rho_2)))
  }
  
  return(Checks)
}

