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