library(tidyverse)
source("r/Simulate.R")

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
  # Note: if NumCheck and nuCheck are both set to 0, then the approximate optimum will be returned
  # round = if NULL, uses exact nu value. Otherwise rounds to this many decimal places
  
  ## Find initial optima
  Start <- N.star.V(s_0, s_1, pi_e, rho_1, rho_2, C_1, C_2, C_e,
                                V=DeltaToV(theta, Delta, alpha), round)
  if (NumCheck==0 & nuCheck==0) {
    Opt <- tibble(N=Start["N"],
                  nu=Start["nu"]) %>%
      dplyr::mutate(Power=Delta.N(theta, s_0, s_1, nu, pi_e, rho_1, rho_2, N, alpha, restrict),
                    SE=sqrt(V.N(s_0, s_1, nu, pi_e, rho_1, rho_2, N, restrict)),
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

## Approach 2: Set Maximum Expected (Or Quantile) Budget, then Minimize Overall Power

OptimalPower <- function(B,
                         theta, s_0, s_1,
                         pi_e, rho_1, rho_2,
                         C_1, C_2, C_e,
                         alpha=0.05,
                         restrict=5,
                         NumCheck=6,
                         round=3) {
  ### Parameters:
  # B = maximum budget (in the same units as C_1, C_2, C_e)
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
  # Note: if NumCheck and nuCheck are both set to 0, then the approximate optimum will be returned
  # round = if NULL, uses exact nu value. Otherwise rounds to this many decimal places

  ## Find initial optima
  Start <- N.star.B(s_0, s_1, pi_e, rho_1, rho_2, C_1, C_2, C_e,
                    B, round)
  
  if (NumCheck==0 & nuCheck==0) {
    Opt <- tibble(N=Start["N"],
                  nu=Start["nu"]) %>%
      dplyr::mutate(Power=Delta.N(theta, s_0, s_1, nu, pi_e, rho_1, rho_2, N, alpha, restrict),
                    SE=sqrt(V.N(s_0, s_1, nu, pi_e, rho_1, rho_2, N, restrict)),
                    Cost.E=C.N(C_1, C_2, C_e, nu, pi_e, rho_1, N, quantile=NULL),
                    Ns.E=N*nu,
                    M0.E=Ns.E*rho_1*(1-pi_e)*rho_2,
                    M1.E=N*rho_1*pi_e*rho_2)
  } else {
    Opt <- CheckSet.B(Start["N"],
                      B,
                      theta, s_0, s_1,
                      pi_e, rho_1, rho_2,
                      C_1, C_2, C_e,
                      alpha,
                      restrict,
                      NumCheck,
                      round)
  }
  return(Opt)
}

## Combined Computation and Simulation Design Analysis Functions:

AnalyzeDesigns <- function(Designs,
                           theta, s_0, s_1,
                           pi_e, rho_1, rho_2,  
                           C_1, C_2, C_e,
                           alpha=0.05, restrict=5,
                           NoSims=0) {
  Res <- NULL
  for (i in 1:nrow(Designs)) {
    N <- Designs$N[i]
    nu <- Designs$nu[i]
    Res.i <- tibble(Power=Delta.N(theta, s_0, s_1,
                                  nu, pi_e, rho_1, rho_2, N, alpha, restrict),
                    SE=sqrt(V.N(s_0, s_1, nu, pi_e, rho_1, rho_2, N,
                                restrict)),
                    Cost.E=C.N(C_1, C_2, C_e, nu, pi_e, rho_1, N),
                    Ns.E=N*nu,
                    M0.E=Ns.E*rho_1*(1-pi_e)*rho_2,
                    M1.E=N*rho_1*pi_e*rho_2)
    if (NoSims > 0) {
      names(Res.i) <- paste0("Comp.",names(Res.i))
      Sims <- SimMany(NoSims,
                      N, nu, theta,
                      s_0, s_1,
                      pi_e, rho_1, rho_2,
                      C_1, C_2, C_e,
                      alpha)
      SimOut <- c(apply(Sims %>% 
                        dplyr::select(Est, Bias, AbsBias, SE, 
                                      CI.Cvg, Sig, Cost),
                      2, mean),
                  Est.SD=sd(Sims$Est))
      names(SimOut) <- paste0("Sim.", names(SimOut))
      Res.i <- Res.i %>% bind_cols(as_tibble_row(SimOut) %>%
                                     dplyr::rename(Sim.Power=Sim.Sig))
    }
    Res <- Res %>% bind_rows(Res.i)
  }
  return(Designs %>% bind_cols(Res))
}

AnalyzeDesignsBin <- function(Designs,
                           pi_0, pi_1,
                           pi_e, rho_1, rho_2,  
                           C_1, C_2, C_e,
                           alpha=0.05, restrict=5,
                           NoSims=0) {
  Res <- NULL
  
  theta <- abs(log((pi_1/(1-pi_1))/(pi_0/(1-pi_0))))
  s_0 <- sqrt(1/(pi_0*(1-pi_0))) ## Lachin, p. 27
  s_1 <- sqrt(1/(pi_1*(1-pi_1)))
  
  for (i in 1:nrow(Designs)) {
    N <- Designs$N[i]
    nu <- Designs$nu[i]
    Res.i <- tibble(Power=Delta.N(theta, s_0, s_1,
                                  nu, pi_e, rho_1, rho_2, N, alpha, restrict),
                    SE=sqrt(V.N(s_0, s_1, nu, pi_e, rho_1, rho_2, N,
                                restrict)),
                    Cost.E=C.N(C_1, C_2, C_e, nu, pi_e, rho_1, N),
                    Ns.E=N*nu,
                    M0.E=Ns.E*rho_1*(1-pi_e)*rho_2,
                    M1.E=N*rho_1*pi_e*rho_2)
    if (NoSims > 0) {
      names(Res.i) <- paste0("Comp.",names(Res.i))
      Sims <- SimManyBin(NoSims,
                         N, nu, 
                         pi_0, pi_1,
                         pi_e, rho_1, rho_2,
                         C_1, C_2, C_e,
                         alpha)
      SimOut <- c(apply(Sims %>% 
                          dplyr::select(ends_with(c("Est", "Bias", "AbsBias", "SE", 
                                        "CI.Cvg", "Sig", "Cost"))),
                        2, mean, na.rm=TRUE),
                  RD.Est.SD=sd(Sims$RD.Est, na.rm=TRUE),
                  OR.Est.SD=sd(Sims$OR.Est, na.rm=TRUE))
      names(SimOut) <- paste0("Sim.", names(SimOut))
      Res.i <- Res.i %>% bind_cols(as_tibble_row(SimOut) %>%
                                     dplyr::rename(Sim.RD.Power=Sim.RD.Sig,
                                                   Sim.OR.Power=Sim.OR.Sig))
    }
    Res <- Res %>% bind_rows(Res.i)
  }
  return(Designs %>% bind_cols(Res))
}

