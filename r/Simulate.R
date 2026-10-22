library(tidyverse)
source("r/KeyFunctions.R")
source("r/Optimize.R")

## At Interim Analysis:

SimOneInt <- function(
    N0, N1,
    pi_0, pi_1,
    rho_2,
    alpha=0.05,
    FailStart=0) {
  
  Failures <- FailStart
  
  M0 <- rbinom(1, size=N0, prob=rho_2)
  M1 <- rbinom(1, size=N1, prob=rho_2)
  
  if (M0 == 0 | M1 == 0) {
    warning("Trial simulation resulted in impermissible M0 or M1 size")
    Failures <- Failures+1
    return(SimOneInt(N0, N1, pi_0, pi_1, rho_2, alpha, FailStart=Failures))
  }
  
  N_Vals <- c(Failures=Failures,
              N0=N0, N1=N1,
              M0=M0, M1=M1,
              SD.M0M1=sqrt(V.M0M1(pi_0,pi_1,M0,M1)),
              Power.M0M1=VToDelta(pi_0,pi_1,
                                  V.M0M1(pi_0,pi_1,M0,M1),
                                  alpha))
  
  Y0 <- rbinom(n=1, size=M0, prob=pi_0)
  Y1 <- rbinom(n=1, size=M1, prob=pi_1)
  
  dat <- tibble(X=rep(c(0,1), times=c(M0,M1)),
                Y=rep(c(0,1,0,1), times=c(M0-Y0,Y0,M1-Y1,Y1)))
  LogReg <- glm(Y~X, family=binomial, data=dat)
  OR.CI <- suppressMessages(exp(confint(LogReg, level=1-alpha, parm="X")))
  Coefs <- summary(LogReg)$coefficients
  Fisher <- fisher.test(x=matrix(c(Y1,M1-Y1,Y0,M0-Y0), nrow=2, ncol=2),
                        conf.level=1-alpha/2)
  
  p0hat <- Y0/N0
  p1hat <- Y1/N1
  RD.SE <- sqrt(p1hat*(1-p1hat)/M1+p0hat*(1-p0hat)/M0)
  RD.CI <- p1hat-p0hat+RD.SE*c(-qnorm(1-alpha/2),qnorm(1-alpha/2))
  
  RD.true <- pi_1-pi_0
  OR.true <- (pi_1/(1-pi_1))/(pi_0/(1-pi_0))
  
  Res_Vals <- c(Y0=Y0, Y1=Y1,
                p0.Est=p0hat, 
                p1.Est=p1hat,
                RD.Est=p1hat-p0hat, 
                RR.Est=p1hat/p0hat, 
                OR.Est=(p1hat/(1-p1hat))/(p0hat/(1-p0hat)),
                RD.SE=RD.SE,
                RD.p=2*pnorm(abs(p1hat-p0hat)/RD.SE, lower.tail=FALSE),
                RD.CIL=RD.CI[1],
                RD.CIU=RD.CI[2],
                RD.CI.inc=(RD.true > RD.CI[1]) & (RD.true < RD.CI[2]),
                OR.SE=unname(Coefs["X","Std. Error"]),
                OR.p=unname(Coefs["X","Pr(>|z|)"]),
                OR.CIL=unname(OR.CI[1]),
                OR.CIU=unname(OR.CI[2]),
                OR.CI.inc=((OR.true > unname(OR.CI[1])) & (OR.true < unname(OR.CI[2]))),
                Fisher.p=Fisher$p.value,
                Fisher.CIL=Fisher$conf.int[1],
                Fisher.CIU=Fisher$conf.int[2],
                Fisher.CI.inc=((OR.true > Fisher$conf.int[1]) & (OR.true < Fisher$conf.int[2])))
  
  return(c(N_Vals,Res_Vals))
}

SimManyInt <- function(NoSims,
                       N0, N1,
                       pi_0, pi_1,
                       rho_2,
                       alpha=0.05) {
  Res <- as_tibble(t(replicate(n=NoSims, expr=SimOneInt(N0, N1,
                                                        pi_0, pi_1,
                                                        rho_2,
                                                        alpha=0.05))))
  
  Res$SD.N1=sqrt(V.N0N1(pi_0,pi_1,N0,N1,rho_2))
  Res$Power.N1=Delta.N0N1(pi_0,pi_1,N0,N1,rho_2,alpha)
  return(Res %>% dplyr::select(Failures,N0,N1,M0,M1,
                               SD.N1,Power.N1,everything()))
}

## As of Beginning of Study:

SimOne <- function(
    N, N0,
    pi_0, pi_1,
    rho_1, pi_e, rho_2,
    C_1, C_2, C_e,
    alpha=0.05,
    FailStart=0) {
  
  Failures <- FailStart
  
  M <- rbinom(1, size=N, prob=rho_1)
  N1 <- rbinom(1, size=M, prob=pi_e)
  
  if (M <= 0 | N1 <= 0 | (N0+N1) > M) {
    warning("Trial simulation resulted in impermissible N1 size")
    Failures <- Failures+1
    return(SimOne(N, N0,
                  pi_0, pi_1,
                  rho_1, pi_e, rho_2,
                  C_1, C_2, C_e,
                  alpha,
                  FailStart=Failures))
  } else {
    return(c(SimOneInt(N0, N1, pi_0, pi_1, rho_2, alpha,
                     FailStart=Failures),
             C.N1=C.N0N1N(C_1, C_2, C_e, N0, N1, N)))
  }
}

SimMany <- function(NoSims,
                    N, N0,
                    pi_0, pi_1,
                    rho_1, pi_e, rho_2,
                    C_1, C_2, C_e,
                    alpha=0.05) {
  as_tibble(t(replicate(n=NoSims, expr=SimOne(N, N0,
                                              pi_0, pi_1,
                                              rho_1, pi_e, rho_2,
                                              C_1, C_2, C_e,
                                              alpha=0.05))))
}

## As of Beginning of Study with Interim Adjustment:
### Warning: this is a bit slow to run (~3 mins per 100 simulations)

SimOneAdj <- function(
    N, N0_Max=NULL, DeltaAdj,
    pi_0, pi_1,
    rho_1, pi_e, rho_2,
    C_1, C_2, C_e,
    alpha=0.05,
    FailStart=0) {
  
  Failures <- FailStart
  
  M <- rbinom(1, size=N, prob=rho_1)
  N1 <- rbinom(1, size=M, prob=pi_e)
  if (is.null(N0_Max)) {
    N0 <- M-N1
  } else {
    N0 <- min(N0_Max, M-N1)
  }
  
  
  if (M <= 0 | N1 <= 0 | N0 <= 0) {
    warning("Trial simulation resulted in impermissible M, N1, or N0 size")
    Failures <- Failures+1
    return(SimOneAdj(N, N0, DeltaAdj,
                  pi_0, pi_1,
                  rho_1, pi_e, rho_2,
                  C_1, C_2, C_e,
                  alpha,
                  FailStart=Failures))
  } else {
    Adj <- OptimalIntBudget(DeltaAdj, pi_0, pi_1, rho_2, C_2, N0, N1, alpha, NumCheck=2)
    return(c(SimOneInt(Adj$N0, Adj$N1, pi_0, pi_1, rho_2, alpha,
                       FailStart=Failures),
             C.N1=N*C_1+N1*C_e+Adj$C.Remaining))
  }
}

SimManyAdj <- function(NoSims,
                    N, N0_Max=NULL, Delta.Adj,
                    pi_0, pi_1,
                    rho_1, pi_e, rho_2,
                    C_1, C_2, C_e,
                    alpha=0.05) {
  as_tibble(t(replicate(n=NoSims, expr=SimOneAdj(N, N0_Max, Delta.Adj,
                                              pi_0, pi_1,
                                              rho_1, pi_e, rho_2,
                                              C_1, C_2, C_e,
                                              alpha=0.05))))
}

## For Comparison Prospective Cohort Study:

SimOneComp <- function(N,
                       pi_0, pi_1,
                       rho_1, pi_e, rho_2,
                       C_1, C_2, C_e,
                       alpha=0.05,
                       FailStart=0) {
  M <- rbinom(n=1, size=N, prob=rho_1)
  N1 <- rbinom(n=1, size=M, prob=pi_e)
  N0 <- M-N1
  C.Full <- C_1*N + C_2*M + C_e*N1
  if (N1 <= 0 | N0 <= 0) {
    Failures = FailStart + 1
    SimOneComp(N, pi_0, pi_1, rho_1, pi_e, rho_2, alpha,
               FailStart=Failures)
  } else {
    Res <- SimOneInt(N0, N1,
                     pi_0, pi_1,
                     rho_2,
                     alpha,
                     FailStart=Failures)
    return(c(Res, C.Full=C.Full))
  }
}

SimManyComp <- function(NoSims,
                        N,
                        pi_0, pi_1,
                        rho_1, pi_e, rho_2,
                        C_1, C_2, C_e,
                        alpha=0.05) {
  as_tibble(t(replicate(n=NoSims, expr=SimOneComp(N,
                                                  pi_0, pi_1,
                                                  rho_1, pi_e, rho_2,
                                                  C_1, C_2, C_e,
                                                  alpha=0.05,
                                                  FailStart=0))))
}
