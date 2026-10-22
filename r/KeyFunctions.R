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
C.N0N1N <- function(C_1, C_2, C_e, N_0, N_1, N) {
  N*C_1+(N_0+N_1)*C_2+N_1*C_e
}
C.N0N <- function(C_1, C_2, C_e, N_0, N, rho_1, pi_e, quantile=NULL) {
  if (is.null(quantile)) {
    N*C_1+N_0*C_2+N*rho_1*pi_e*(C_2+C_e)
  } else {
    N*C_1+N_0*C_2+qbinom(p=quantile, size=N, prob=rho_1*pi_e)*(C_2+C_e)
  }
}
Delta.N0N1 <- function(pi_0, pi_1, N_0, N_1, rho_2, alpha=0.05) {
  mat <- cross_join(tibble(M0=1:N_0), tibble(M1=1:N_1)) %>%
    dplyr::mutate(P0=dbinom(M0, size=N_0, prob=rho_2),
                  P1=dbinom(M1, size=N_1, prob=rho_2))
  mat$V = apply(mat, 1, function(row)
    V.M0M1(pi_0,pi_1,row["M0"],row["M1"]))
  mat$Delta <- apply(mat, 1, function(row)
    VToDelta(pi_0, pi_1, row["V"], alpha))
  return(sum(mat$P0*mat$P1*mat$Delta))
}

## Approximate Optima Functions
NVtoN0 <- function(pi_0,pi_1,rho_1,pi_e,rho_2,N,V) {
  (N*rho_1*pi_e*pi_0*(1-pi_0))/(V*N*rho_1*pi_e*rho_2-pi_1*(1-pi_1))
}

NBtoN0 <- function(rho_1,pi_e,C_1,C_2,C_e,N,B) {
  (B-N*(C_1+rho_1*pi_e*C_2+rho_1*pi_e*C_e))/C_2
}

N.star.V <- function(pi_e, pi_0, pi_1, C_1, C_2, C_e, rho_1, rho_2, V) {
  (pi_1*(1-pi_1)+sqrt(C_2*rho_1*pi_e*pi_0*(1-pi_0)*pi_1*(1-pi_1)/(C_1+rho_1*pi_e*C_2+rho_1*pi_e*C_e)))/(V*rho_1*pi_e*rho_2) 
}

N.star.Delta <- function(pi_e, pi_0, pi_1, C_1, C_2, C_e, rho_1, rho_2, Delta, alpha=0.05) {
  N.star.V(pi_e, pi_0, pi_1, C_1, C_2, C_e, rho_1, rho_2,
           DeltaToV(pi_0, pi_1, Delta, alpha))
}

N.star.B <- function(pi_e, pi_0, pi_1, C_1, C_2, C_e, rho_1, B) {
  B*(pi_1*(1-pi_1)-sqrt(pi_0*(1-pi_0)*pi_1*(1-pi_1)*C_2*rho_1*pi_e/(C_1+rho_1*pi_e*C_2+rho_1*pi_e*C_e)))/((C_1+rho_1*pi_e*C_2+rho_1*pi_e*C_e)*pi_1*(1-pi_1)-C_2*rho_1*pi_e*pi_0*(1-pi_0))
}

## Dynamic Optimum-Finding Functions
CheckSet.P <- function(Starts,
                       Delta,
                       pi_0,pi_1,
                       rho_1,pi_e,rho_2,
                       C_1,C_2,C_e,
                       alpha=0.05,
                       NumCheck=20,
                       quantile=NULL) {
  Checks <- cross_join(tibble(N=(Starts$N:(Starts$N+NumCheck))),
                       tibble(N0=(Starts$N0:(Starts$N0+NumCheck))))
  Checks$Power.E <- apply(Checks, 1, FUN=function(row) VToDelta(pi_0, pi_1,
                                                                V.N0N(pi_0,pi_1,row["N0"],row["N"],
                                                                      rho_1,pi_e,rho_2),
                                                                alpha))
  if (max(Checks$Power.E) >= Delta) {
    Checks <- Checks %>% dplyr::filter(Power.E >= Delta)
    Checks$C.val <- apply(Checks, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile))
    Checks <- Checks %>% dplyr::filter(C.val==min(Checks$C.val))

    Checks$C.E <- apply(Checks, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=NULL))
    Checks$SD.E <- apply(Checks, 1, 
                         FUN=function(row) sqrt(V.N0N(pi_0,pi_1,row["N0"],row["N"],rho_1,pi_e,rho_2)))
    Checks$C.75 <- apply(Checks, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=0.75))
    Checks$C.90 <- apply(Checks, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=0.90))
    
    Checks <- Checks %>%
      dplyr::mutate(N1.E=N*pi_e,
                    N1=round(N1.E))
    Checks$C.N1 <- apply(Checks, 1, FUN=function(row) 
      C.N0N1N(C_1, C_2, C_e, row["N0"], row["N1"], row["N"]))
    Checks$SD.N1 <- apply(Checks, 1, FUN=function(row) 
      sqrt(V.N0N1(pi_0,pi_1,row["N0"],row["N1"],rho_2)))
    Checks$Power.N1 <- apply(Checks, 1, FUN=function(row) 
      Delta.N0N1(pi_0, pi_1, row["N0"],row["N1"],rho_2, alpha))
    return(Checks)
  } else {
    stop("No combination found that meets required power. Try increasing NumCheck.")
  }
}

CheckSet.B <- function(Starts,
                       B,
                       pi_0,pi_1,
                       rho_1,pi_e,rho_2,
                       C_1,C_2,C_e,
                       alpha=0.05,
                       NumCheck=5) {
  Checks <- cross_join(tibble(N=(Starts$N-NumCheck):(Starts$N)),
                       tibble(N0=(Starts$N0-NumCheck):(Starts$N0)))
  Checks$Power.E <- apply(Checks, 1, FUN=function(row) VToDelta(pi_0, pi_1,
                                                                V.N0N(pi_0,pi_1,row["N0"],row["N"],
                                                                      rho_1,pi_e,rho_2),
                                                                alpha))
  Checks$C.val <- apply(Checks, 1, FUN=function(row) 
    C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=NULL))
  
  if (min(Checks$C.val) <= B) {
    Checks <- Checks %>% dplyr::filter(C.val <= B)
    Checks <- Checks %>% dplyr::filter(Power.E==max(Checks$Power.E))
    
    Checks$SD.E <- apply(Checks, 1, 
                         FUN=function(row) sqrt(V.N0N(pi_0,pi_1,row["N0"],row["N"],rho_1,pi_e,rho_2)))
    Checks$C.75 <- apply(Checks, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=0.75))
    Checks$C.90 <- apply(Checks, 1, FUN=function(row) 
      C.N0N(C_1, C_2, C_e, row["N0"], row["N"], rho_1, pi_e, quantile=0.90))
    
    Checks <- Checks %>%
      dplyr::mutate(N1.E=N*pi_e,
                    N1=round(N1.E))
    Checks$C.N1 <- apply(Checks, 1, FUN=function(row) 
      C.N0N1N(C_1, C_2, C_e, row["N0"], row["N1"], row["N"]))
    Checks$SD.N1 <- apply(Checks, 1, FUN=function(row) 
      sqrt(V.N0N1(pi_0,pi_1,row["N0"],row["N1"],rho_2)))
    Checks$Power.N1 <- apply(Checks, 1, FUN=function(row) 
      Delta.N0N1(pi_0, pi_1, row["N0"],row["N1"],rho_2, alpha))
    return(Checks)
  } else {
    stop("No combination found that meets required budget Try increasing NumCheck.")
  }
}
