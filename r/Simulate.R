library(tidyverse)
source("r/KeyFunctions.R")
source("r/Optimize.R")

## As of M0 and M1:
SimInt <- function(
    M0, M1,
    theta,
    s_0, s_1,
    pi_e, rho_1, rho_2,
    C_1, C_2, C_e,
    alpha=0.05) {
  Y0 <- rnorm(n=M0, mean=0, sd=s_0)
  Y1 <- rnorm(n=M1, mean=theta, sd=s_1)
  
  dat <- tibble(A=rep(c(0,1), times=c(M0,M1)),
                Y=c(Y0,Y1))
  lm <- lm(Y~A, data=dat)
  ci <- confint(lm, parm="A", level=1-alpha)
  return(c(Y0bar=mean(Y0), Y1bar=mean(Y1),
           Y0sd=sd(Y0), Y1sd=sd(Y1),
           Est=unname(coef(lm)["A"]),
           SE=summary(lm)$coefficients["A","Std. Error"],
           p=summary(lm)$coefficients["A","Pr(>|t|)"],
           CIL=ci[1],
           CIU=ci[2]))
}

## As of Beginning of Study:

SimOne <- function(
    N, nu, theta,
    s_0, s_1,
    pi_e, rho_1, rho_2,
    C_1, C_2, C_e,
    alpha=0.05,
    FailStart=0) {
  
  Failures <- FailStart
  
  Ns <- rbinom(1, size=N, prob=nu)
  NsE <- rbinom(1, size=Ns, prob=rho_1) ## The number of Ns still enrolled at exposure determination
  
  if (NsE <= 0) {
    warning("Trial simulation resulted in impermissible Ns size at exposure")
    Failures <- Failures+1
    return(SimOne(N, nu,
                  s_0, s_1,
                  pi_e, rho_1, rho_2,
                  C_1, C_2, C_e,
                  alpha,
                  FailStart=Failures))
  }
  
  N0 <- rbinom(1, size=NsE, prob=1-pi_e)
  N1 <- rbinom(1, size=N-Ns, prob=rho_1*pi_e)+(NsE-N0)
  
  if (N0 <= 0 | N1 <= 0) {
    warning("Trial simulation resulted in impermissible N0 or N1 size")
    Failures <- Failures+1
    return(SimOne(N, nu,
                  s_0, s_1,
                  pi_e, rho_1, rho_2,
                  C_1, C_2, C_e,
                  alpha,
                  FailStart=Failures))
  }
  
  M0 <- rbinom(1, size=N0, prob=rho_2)
  M1 <- rbinom(1, size=N1, prob=rho_2)
  
  if (M0 <= 0 | M1 <= 0) {
    warning("Trial simulation resulted in impermissible M0 or M1 size")
    Failures <- Failures+1
    return(SimOne(N, nu,
                  s_0, s_1,
                  pi_e, rho_1, rho_2,
                  C_1, C_2, C_e,
                  alpha,
                  FailStart=Failures))
  }
  
  return(c(N0=N0, N1=N1, M0=M0, M1=M1,
           Cost=C_1*N + C_2*(N0+N1) + C_e*N1,
           SimInt(M0, M1,
                  theta,
                  s_0, s_1,
                  pi_e, rho_1, rho_2,
                  C_1, C_2, C_e,
                  alpha)))
}

SimMany <- function(NoSims,
                    N, nu, theta,
                    s_0, s_1,
                    pi_e, rho_1, rho_2,
                    C_1, C_2, C_e,
                    alpha=0.05) {
  as_tibble(t(replicate(n=NoSims, expr=SimOne(N, nu, theta,
                                              s_0, s_1,
                                              pi_e, rho_1, rho_2,
                                              C_1, C_2, C_e,
                                              alpha=0.05,
                                              FailStart=0)))) %>%
    dplyr::mutate(Diff=Est-theta,
                  AbsDiff=abs(Est-theta),
                  CI.Incl=(CIL <= theta & theta <= CIU),
                  Sig=(p < alpha))
}

## For Comparison Prospective Cohort Study:

SimOneComp <- function(N, theta,
                       s_0, s_1,
                       pi_e, rho_1, rho_2,
                       C_1, C_2, C_e,
                       alpha=0.05,
                       FailStart=0) {
  M <- rbinom(n=1, size=N, prob=rho_1)
  if (M <= 0) {
    warning("Trial simulation resulted in impermissible M size")
    Failures <- Failures+1
    return(SimOneComp(N,
                      s_0, s_1,
                      pi_e, rho_1, rho_2,
                      C_1, C_2, C_e,
                      alpha,
                      FailStart=Failures))
    
  }
  
  N1 <- rbinom(n=1, size=M, prob=pi_e)
  N0 <- M-N1
  if (N0 <= 0 | N1 <= 0) {
    warning("Trial simulation resulted in impermissible N0 or N1 size")
    Failures <- Failures+1
    return(SimOneComp(N,
                      s_0, s_1,
                      pi_e, rho_1, rho_2,
                      C_1, C_2, C_e,
                      alpha,
                      FailStart=Failures))
    
  }
  
  M1 <- rbinom(n=1, size=N1, prob=rho_2)
  M0 <- rbinom(n=1, size=N0, prob=rho_2)
  
  if (M0 <= 0 | M1 <= 0) {
    warning("Trial simulation resulted in impermissible M0 or M1 size")
    Failures <- Failures+1
    return(SimOneComp(N,
                      s_0, s_1,
                      pi_e, rho_1, rho_2,
                      C_1, C_2, C_e,
                      alpha,
                      FailStart=Failures))
                  
  }
  
  return(c(N0=N0, N1=N1, M0=M0, M1=M1,
           Cost=C_1*N + C_2*(N0+N1) + C_e*N1,
           SimInt(M0, M1,
                  theta,
                  s_0, s_1,
                  pi_e, rho_1, rho_2,
                  C_1, C_2, C_e,
                  alpha)))
}

SimManyComp <- function(NoSims,
                        N, theta,
                        s_0, s_1,
                        pi_e, rho_1, rho_2,
                        C_1, C_2, C_e,
                        alpha=0.05) {
  as_tibble(t(replicate(n=NoSims, expr=SimOneComp(N, theta,
                                                  s_0, s_1,
                                                  pi_e, rho_1, rho_2,
                                                  C_1, C_2, C_e,
                                                  alpha,
                                                  FailStart=0)))) %>%
    dplyr::mutate(Diff=Est-theta,
                  AbsDiff=abs(Est-theta),
                  CI.Incl=(CIL <= theta & theta <= CIU),
                  Sig=(p < alpha))
}
