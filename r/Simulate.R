library(tidyverse)
source("R/KeyFunctions.R")
options(scipen=8)

SimOne <- function(
    N, N0,
    pi_exp,
    p0, p1,
    C1, C2, C3,
    alpha=0.05) {
  Failures <- 0
  N1 <- rbinom(1, size=N, prob=pi_exp)
  
  if (N1 <= 0 | (N0 + N1) > N) {
    warning("Trial simulation resulted in impermissible N1 size")
    Failures <- Failures+1
    return(SimOne(N, N0, pi_exp, p0, p1, C1, C2, C3, alpha))
  }
  
  N_Vals <- c(Failures=Failures,
              N0=N0, N1=N1,
               SD.N1=sqrt(Var.N1(p0,p1,N0,N1)),
               Power.N1=Power.N1(p0,p1,N0,N1,alpha),
               B.N1=B.N1(C1,C2,C3,N0,N1,N))
  
  Y0 <- rbinom(n=1, size=N0, prob=p0)
  Y1 <- rbinom(n=1, size=N1, prob=p1)
  
  dat <- tibble(X=rep(c(0,1), times=c(N0,N1)),
                Y=rep(c(0,1,0,1), times=c(N0-Y0,Y0,N1-Y1,Y1)))
  LogReg <- glm(Y~X, family=binomial, data=dat)
  OR.CI <- suppressMessages(exp(confint(LogReg, level=1-alpha, parm="X")))
  Coefs <- summary(LogReg)$coefficients
  Fisher <- fisher.test(x=matrix(c(Y1,N1-Y1,Y0,N0-Y0), nrow=2, ncol=2),
                        conf.level=1-alpha/2)
  
  p0hat <- Y0/N0
  p1hat <- Y1/N1
  RD.SE <- sqrt(p1hat*(1-p1hat)/N1+p0hat*(1-p0hat)/N0)
  RD.CI <- p1hat-p0hat+RD.SE*c(-qnorm(1-alpha/2),qnorm(1-alpha/2))
  
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
                RD.CI.inc=((p1-p0) > RD.CI[1]) & ((p1-p0) < RD.CI[2]),
                OR.SE=unname(Coefs["X","Std. Error"]),
                OR.p=unname(Coefs["X","Pr(>|z|)"]),
                OR.CIL=unname(OR.CI[1]),
                OR.CIU=unname(OR.CI[2]),
                OR.CI.inc=(((p1/(1-p1))/(p0/(1-p0)) > unname(OR.CI[1])) & ((p1/(1-p1))/(p0/(1-p0)) < unname(OR.CI[2]))),
                Fisher.p=Fisher$p.value,
                Fisher.CIL=Fisher$conf.int[1],
                Fisher.CIU=Fisher$conf.int[2],
                Fisher.CI.inc=(((p1/(1-p1))/(p0/(1-p0)) > Fisher$conf.int[1]) & ((p1/(1-p1))/(p0/(1-p0)) < Fisher$conf.int[2])))
  
  return(c(N_Vals,Res_Vals))
}

SimMany <- function(NoSims,
                    N, N0,
                    pi_exp,
                    p0, p1,
                    C1, C2, C3,
                    alpha=0.05) {
  as_tibble(t(replicate(n=NoSims, expr=SimOne(N, N0, pi_exp, p0, p1, C1, C2, C3, alpha))))
}


set.seed(20261008)
Sim80 <- SimMany(NoSims=2000,
                 N=1594, N0=184,
                 pi_exp=0.05, p0=0.1, p1=0.25,
                 C1=500, C2=200, C3=1000,
                 alpha=0.05)
apply(Sim80 %>% dplyr::select(ends_with("N1"),ends_with(".Est"),ends_with(".SE"),ends_with(".CI.inc")), 2, mean)
apply(Sim80 %>% dplyr::select(ends_with(".p")), 2, function(x) mean(x < 0.05))


