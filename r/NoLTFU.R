## Ignoring Loss To Follow-Up Throughout

library(tidyverse)

## Core Functions
Var.N1 <- function(p0,p1,N0,N1) {p0*(1-p0)/N0+p1*(1-p1)/N1}
Var.E <- function(p0,p1,N0,pi_exp,N) {p0*(1-p0)/N0+p1*(1-p1)*sum(dbinom(x=1:N, size=N, prob=pi_exp)*1/(1:N))}
Power.N1 <- function(p0,p1,N0,N1,alpha=0.05) {pnorm(qnorm(1-alpha/2)-abs(p1-p0)/(sqrt(p1*(1-p1)/N1+p0*(1-p0)/N0)), 0, 1, lower.tail=FALSE)}
Power.E <- function(p0,p1,N0,pi_exp,N,alpha=0.05) {sum(dbinom(x=1:N, size=N, prob=pi_exp)*sapply(1:N, function(x) Power.N1(p0,p1,N0,x,alpha)))}
B.N1 <- function(C1,C2,C3,N0,N1,N) {N*C1+N1*C2+(N1+N0)*C3}
B.E <- function(C1,C2,C3,N0,pi_exp,N) {N*C1+N*pi_exp*C2+(N*pi_exp+N0)*C3}
B.quant <- function(C1,C2,C3,N0,pi_exp,N,quantile) {N*C1+qbinom(p=quantile, size=N, prob=pi_exp)*(C2+C3)+N0*C3}
Var.FromPower <- function(p0,p1,power,alpha=0.05) {(p1-p0)^2/(qnorm(1-alpha/2)-qnorm(1-power))^2}
Power.FromVar <- function(p0,p1,V,alpha=0.05) {1-pnorm(qnorm(1-alpha/2)-abs(p1-p0)/sqrt(V))}

## Non-Random Optimum-Finding Functions
N.star.V <- function(C1,C2,C3,p0,p1,pi_exp,V) {p1*(1-p1)/(pi_exp*V) + sqrt(C3*p0*(1-p0)*p1*(1-p1)/(pi_exp*V^2*(C1+pi_exp*C2+pi_exp*C3)))}
N.star.B <- function(C1,C2,C3,p0,p1,pi_exp,B) {(B*p1*(1-p1)-sqrt(B^2*C3/(C1+pi_exp*C2+pi_exp*C3)*pi_exp*p0*(1-p0)*p1*(1-p1)))/((C1+pi_exp*C2+pi_exp*C3)*p1*(1-p1)-C3*pi_exp*p0*(1-p0))}

## Dynamic Optimum-Finding Functions
CheckSet.P <- function(Starts,
                       Delta,pi_exp,
                       p0,p1,
                       C1,C2,C3,
                       alpha=0.05,
                       NumCheck=20,
                       quantile=NULL) {
  Checks <- cross_join(tibble(N=(Starts$N:(Starts$N+NumCheck))),
                       tibble(N0=(Starts$N0:(Starts$N0+NumCheck))))
  Checks$Power.E <- apply(Checks, 1, FUN=function(row) Power.E(p0,p1,row["N0"],pi_exp,row["N"],alpha))
  if (max(Checks$Power.E) >= Delta) {
    Checks <- Checks %>% dplyr::filter(Power.E >= Delta)
    Checks$B.E <- apply(Checks, 1, FUN=function(row) B.E(C1,C2,C3,row["N0"],pi_exp,row["N"]))
    if (is.null(quantile)) { # Minimizes mean budget
      Checks <- Checks %>% dplyr::filter(B.E==min(Checks$B.E))
    } else {
      Checks$B.quantile <- apply(Checks, 1, FUN=function(row) B.quant(C1,C2,C3,row["N0"],pi_exp,row["N"],quantile))
      Checks <- Checks %>% dplyr::filter(B.quantile==min(Checks$B.quantile))
    }
    Checks$SD.E <- apply(Checks, 1, FUN=function(row) sqrt(Var.E(p0,p1,row["N0"],pi_exp,row["N"])))
    Checks$N1.E <- Checks$N*pi_exp
    Checks$B.75 <- apply(Checks, 1, FUN=function(row) B.quant(C1,C2,C3,row["N0"],pi_exp,row["N"],0.75))
    Checks$B.90 <- apply(Checks, 1, FUN=function(row) B.quant(C1,C2,C3,row["N0"],pi_exp,row["N"],0.90))
    Checks <- Checks %>%
      dplyr::mutate(N1=round(N1.E))
    Checks$B.N1 <- apply(Checks, 1, FUN=function(row) B.N1(C1,C2,C3,row["N0"],row["N1"],row["N"]))
    Checks$SD.N1 <- apply(Checks, 1, FUN=function(row) sqrt(Var.N1(p0,p1,row["N0"],row["N1"])))
    Checks$Power.N1 <- apply(Checks, 1, FUN=function(row) Power.N1(p0,p1,row["N0"],row["N1"],alpha))
    return(Checks)
  } else {
    stop("No combination found that meets required power. Try increasing NumCheck.")
  }
}

CheckSet.B <- function(Starts,
                       B,pi_exp,
                       p0,p1,
                       C1,C2,C3,
                       alpha=0.05,
                       NumCheck=5) {
  Checks <- cross_join(tibble(N=(Starts$N-NumCheck):(Starts$N)),
                       tibble(N0=(Starts$N0-NumCheck):(Starts$N0)))
  Checks$B.E <- apply(Checks, 1, FUN=function(row) B.E(C1,C2,C3,row["N0"],pi_exp,row["N"]))
  if (min(Checks$B.E) <= B) {
    Checks <- Checks %>% dplyr::filter(B.E <= B)
    Checks$Power.E <- apply(Checks, 1, FUN=function(row) Power.E(p0,p1,row["N0"],pi_exp,row["N"],alpha))
    Checks <- Checks %>% dplyr::filter(Power.E==max(Checks$Power.E))
    Checks$SD.E <- apply(Checks, 1, FUN=function(row) sqrt(Var.E(p0,p1,row["N0"],pi_exp,row["N"])))
    Checks$N1.E <- Checks$N*pi_exp
    Checks$B.75 <- apply(Checks, 1, FUN=function(row) B.quant(C1,C2,C3,row["N0"],pi_exp,row["N"],quantile=0.75))
    Checks$B.90 <- apply(Checks, 1, FUN=function(row) B.quant(C1,C2,C3,row["N0"],pi_exp,row["N"],quantile=0.90))
    ComboBest$N1 <- round(ComboBest$N1.E)
    ComboBest$B.N1 <- apply(ComboBest, 1, FUN=function(row) B.N1(C1,C2,C3,row["N0"],row["N1"],row["N"]))
    ComboBest$SD.N1 <- apply(ComboBest, 1, FUN=function(row) sqrt(Var.N1(p0,p1,row["N0"],row["N1"])))
    ComboBest$Power.N1 <- apply(ComboBest, 1, FUN=function(row) Power.N1(p0,p1,row["N0"],row["N1"],alpha))
    return(Checks)
  } else {
    stop("No combination found that meets required budget Try increasing NumCheck.")
  }
}


### Simulate:

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


