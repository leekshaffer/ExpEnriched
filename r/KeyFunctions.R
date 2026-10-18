library(tidyverse)

## Core Functions
DeltaToV <- function(pi_0,pi_1,Delta,alpha=0.05) {(pi_1-pi_0)^2/(qnorm(1-alpha/2)-qnorm(1-Delta))^2}
VToDelta <- function(pi_0,pi_1,V,alpha=0.05) {1-pnorm(qnorm(1-alpha/2)-abs(pi_1-pi_0)/sqrt(V))}
V.M0M1 <- function(pi_0,pi_1,M_0,M_1) {pi_0*(1-pi_0)/M_0+pi_1*(1-pi_1)/M_1}
V.N0N1 <- function(pi_0,pi_1,N_0,N_1,rho_2) {
  N_0 <- floor(N_0)
  N_1 <- floor(N_1)
  pi_0*(1-pi_0)*sum(sapply(1:N_0, function(x) 1/x*dbinom(x, N_0, rho_2))) +
    pi_1*(1-pi_1)*sum(sapply(1:N_1, function(x) 1/x*dbinom(x, N_1, rho_2)))
}
V.N0N <- function(pi_0,pi_1,N_0,N,rho_1,pi_e,rho_2) {
  N_0 <- floor(N_0)
  N <- floor(N)
  pi_0*(1-pi_0)*sum(sapply(1:N_0, function(x) 1/x*dbinom(x, N_0, rho_2))) +
    pi_1*(1-pi_1)*sum(sapply(1:N, function(y) 1/y*dbinom(y, N, rho_1*pi_e*rho_2)))
}
C.NN0N1 <- function(C_1, C_2, C_e, N, N_0, N_1) {
  N*C_1+(N_0+N_1)*C_2+N_1*C_e
}
C.NN0 <- function(C_1, C_2, C_e, N, N_0, rho_1, pi_e) {
  N*C_1+N_0*C_2+N*rho_1*pi_e*(C_2+C_e)
}

## Approximate Optima Functions
N.star.V2 <- function(pi_e, pi_0, pi_1, C_1, C_2, C_e, rho_1, rho_2, V) {
  (pi_1*(1-pi_1)+sqrt(C_2*rho_1*pi_e*pi_0*(1-pi_0)*pi_1*(1-pi_1)/(C_1+rho_1*pi_e*C_2+rho_1*pi_e*C_e)))/(V*rho_1*pi_e*rho_2) 
}

N.star.B2 <- function(pi_e, pi_0, pi_1, C_1, C_2, C_e, rho_1, B) {
  B*(pi_1*(1-pi_1)-sqrt(pi_0*(1-pi_0)*pi_1*(1-pi_1)*C_2*rho_1*pi_e/(C_1+rho_1*pi_e*C_2+rho_1*pi_e*C_e)))/((C_1+rho_1*pi_e*C_2+rho_1*pi_e*C_e)*pi_1*(1-pi_1)-C_2*rho_1*pi_e*pi_0*(1-pi_0))
}

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