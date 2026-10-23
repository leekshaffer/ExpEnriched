library(tidyverse)
options(warn=1)

## Core Functions
DeltaToV <- function(theta,Delta,alpha=0.05) {(theta)^2/(qnorm(1-alpha/2)-qnorm(1-Delta))^2}
VToDelta <- function(theta,V,alpha=0.05) {1-pnorm(qnorm(1-alpha/2)-abs(theta)/sqrt(V))}

V.M0M1 <- function(s_0,s_1,M_0,M_1) {s_0^2/M_0+s_1^2/M_1}
V.N0N1 <- function(s_0,s_1,rho_2,N_0,N_1) {
  N_0 <- floor(N_0)
  N_1 <- floor(N_1)
  s_0^2*sum(sapply(1:N_0, function(i) dbinom(i, size=N_0, prob=rho_2)/i)) +
    s_1^2*sum(sapply(1:N_1, function(j) dbinom(j, size=N_1, prob=rho_2)/j))
}
V.N <- function(s_0,s_1,nu,pi_e,rho_1,rho_2,N,
                restrict=NULL) {
  if (is.null(restrict)) {
    sum(sapply(1:(N-1),
               function(k) dbinom(k, size=N, prob=pi_e)*
                 (s_0^2*sum(sapply(1:(N-k), function(i) dbinom(i, size=N-k, prob=nu*rho_1*rho_2)/i))+
                    s_1^2*sum(sapply(1:k, function(j) dbinom(j, size=k, prob=rho_1*rho_2)/j)))))
  } else {
    sum(sapply(max(1,floor(N*pi_e-restrict*sqrt(N*pi_e*(1-pi_e)))):min(N-1,ceiling(N*pi_e+restrict*sqrt(N*pi_e*(1-pi_e)))),
               function(k) dbinom(k, size=N, prob=pi_e)*
                 (s_0^2*sum(sapply(max(1,floor((N-k)*nu*rho_1*rho_2-restrict*sqrt((N-k)*nu*rho_1*rho_2))):
                                     min(N-k,ceiling((N-k)*nu*rho_1*rho_2+restrict*sqrt((N-k)*nu*rho_1*rho_2))), 
                                   function(i) dbinom(i, size=N-k, prob=nu*rho_1*rho_2)/i))+
                    s_1^2*sum(sapply(1:k, function(j) dbinom(j, size=k, prob=rho_1*rho_2)/j)))))
  }
}

C.N0N1N <- function(C_1, C_2, C_e, N_0, N_1, N) {
  N*C_1+(N_0+N_1)*C_2+N_1*C_e
}

C.N <- function(C_1, C_2, C_e, nu, pi_e, rho_1, N, quantile=NULL) {
  if (is.null(quantile)) {
    N*C_1+N*nu*rho_1*(1-pi_e)*C_2+N*rho_1*pi_e*(C_2+C_e)
  } else {
    stop("Quantile not implemented yet for this sampling method.")
    # N*C_1+N_0*C_2+qbinom(p=quantile, size=N, prob=rho_1*pi_e)*(C_2+C_e)
  }
}

Delta.N0N1 <- function(theta,s_0,s_1,rho_2,N_0,N_1,
                       alpha=0.05,
                       restrict=NULL) {
  N_0 <- floor(N_0)
  N_1 <- floor(N_1)
  
  if (is.null(restrict)) {
    v0 <- 1:N_0
    v1 <- 1:N_1
  } else {
    v0 <- max(1,floor(N_0*rho_2-restrict*sqrt(N_0*rho_2*(1-rho_2)))):min(N_0, ceiling(N_0*rho_2+restrict*sqrt(N_0*rho_2*(1-rho_2))))
    v1 <- max(1,floor(N_1*rho_2-restrict*sqrt(N_1*rho_2*(1-rho_2)))):min(N_1, ceiling(N_1*rho_2+restrict*sqrt(N_1*rho_2*(1-rho_2))))
  }
  mat <- cross_join(tibble(M0=v0), tibble(M1=v1)) %>%
    dplyr::mutate(P0=dbinom(M0, size=N_0, prob=rho_2),
                  P1=dbinom(M1, size=N_1, prob=rho_2))
  sum(apply(mat, 1, function(row) row["P0"]*row["P1"]*
              VToDelta(theta, V.M0M1(s_0,s_1,row["M0"],row["M1"]), alpha)))
}

Delta.N <- function(theta, s_0, s_1, nu, pi_e, rho_1, rho_2,
                    N, alpha=0.05,
                    restrict=NULL) {
  N <- floor(N)
  
  if (is.null(restrict)) {
    v <- 1:(N-1)
  } else {
    v <- max(1,floor(N*pi_e-restrict*sqrt(N*pi_e*(1-pi_e)))):min(N-1,ceiling(N*pi_e+restrict*sqrt(N*pi_e*(1-pi_e))))
  }
  
  Dk <- function(k) {
    if (is.null(restrict)) {
      v0 <- 1:(N-k)
      v1 <- 1:k
    } else {
      v0 <- max(1,floor((N-k)*nu*rho_1*rho_2-restrict*sqrt((N-k)*nu*rho_1*rho_2*(1-nu*rho_1*rho_2)))):
        min(N-k,ceiling((N-k)*nu*rho_1*rho_2+restrict*sqrt((N-k)*nu*rho_1*rho_2*(1-nu*rho_1*rho_2))))
      v1 <- max(1,floor(k*rho_1*rho_2-restrict*sqrt(k*rho_1*rho_2*(1-rho_1*rho_2)))):
        min(k,ceiling(k*rho_1*rho_2+restrict*sqrt(k*rho_1*rho_2*(1-rho_1*rho_2))))
    }
    mat <- cross_join(tibble(M0=v0), tibble(M1=v1)) %>%
      dplyr::mutate(P0=dbinom(M0, size=N-k, prob=nu*rho_1*rho_2),
                    P1=dbinom(M1, size=k, prob=rho_1*rho_2))
    sum(apply(mat, 1, function(row) row["P0"]*row["P1"]*
                VToDelta(theta,V.M0M1(s_0, s_1, row["M0"], row["M1"]))))
  }
  
  sum(sapply(v, function(k) dbinom(k, size=N, prob=pi_e)*Dk(k)))
}


## Approximate Optima Functions

nu.B <- function(pi_e, rho_1, C_1, C_2, C_e, N, B,
                 round=NULL) {
  D <- C_1+rho_1*pi_e*(C_2+C_e)
  nu <- (B-N*D)/(N*rho_1*(1-pi_e)*C_2)
  if (is.null(round)) {
    nu
  } else {
    floor(nu*10^round)/(10^round)
  }
}
N.star.B <- function(s_0, s_1, pi_e, rho_1, rho_2, C_1, C_2, C_e,
                     B, round=NULL) {
  D <- C_1+rho_1*pi_e*(C_2+C_e)
  if (D*s_1^2 <= rho_1*pi_e*s_0^2*C_2) {
    warning("Cannot guarantee this is a minimum")
  }
  Nstar <- floor(B*(s_1^2-sqrt(rho_1*pi_e*s_0^2*s_1^2*C_2/D))/(D*s_1^2-rho_1*pi_e*s_0^2*C_2))
  return(c(N=Nstar, nu=nu.B(pi_e, rho_1, C_1, C_2, C_e, Nstar, B, round)))
}

nu.V <- function(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e,
                 round=NULL) {
  D <- C_1+rho_1*pi_e*(C_2+C_e)
  nu <- sqrt(s_0^2*pi_e*D)/((1-pi_e)*sqrt(C_2*s_1^2*rho_1))
  if (is.null(round)) {
    return(nu)
  } else {
    return(ceiling(nu*10^round)/(10^round))
  }
}
N.star.V <- function(s_0, s_1, pi_e, rho_1, rho_2, C_1, C_2, C_e,
                     V, round=NULL) {
  D <- C_1+rho_1*pi_e*(C_2+C_e)
  Nstar <- ceiling((s_1^2+sqrt(C_2*s_0^2*s_1^2*rho_1*pi_e/D))/(V*rho_1*pi_e*rho_2))
  return(c(N=Nstar, nu=nu.V(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e, round)))
}

## Functions for Comparison Prospective Design

Delta.Comp <- function(theta, s_0, s_1, pi_e, rho_1, rho_2, N, 
                       alpha=0.05, restrict=5) {
  Delta.N(theta, s_0, s_1, nu=1, pi_e, rho_1, rho_2,
          N, alpha, restrict)
}

N.star.V.Comp <- function(s_0, s_1, pi_e, rho_1, rho_2, V) {
  (s_0^2*pi_e+s_1^2*(1-pi_e))/(pi_e*(1-pi_e)*rho_1*rho_2*V)
}

N.star.Delta.Comp <- function(theta, s_0, s_1, pi_e, rho_1, rho_2, 
                              Delta, alpha=0.05) {
  N.star.V.Comp(s_0, s_1, pi_e, rho_1, rho_2,
                V=DeltaToV(theta, Delta, alpha))
}

N.star.B.Comp <- function(pi_e, rho_1, 
                          C_1, C_2, C_e,
                          B) {
  return(B/(C_1+C_2*rho_1+C_e*rho_1*pi_e))
}

## Dynamic Optimum-Finding Functions
CheckSet.B <- function(StartN,
                       B,
                       theta, s_0, s_1,
                       pi_e, rho_1, rho_2,
                       C_1, C_2, C_e,
                       alpha=0.05,
                       restrict=5,
                       NumCheck=5,
                       round=3) {
  Checks <- tibble(N=(StartN-NumCheck):(StartN+NumCheck)) %>%
    dplyr::mutate(nu=nu.B(pi_e, rho_1, C_1, C_2, C_e, N, B, round),
                  Cost.E=C.N(C_1, C_2, C_e, nu, pi_e, rho_1, N, quantile=NULL)) %>%
    dplyr::filter(Cost.E <= B)
  Checks$Power <- apply(Checks, 1, function(row) {
    Delta.N(theta, s_0, s_1, row["nu"], pi_e, rho_1, rho_2,
            row["N"], alpha, restrict)
  })
  Checks <- Checks %>% dplyr::filter(Power==max(Checks$Power))
  if (nrow(Checks) > 1) {
    Checks <- Checks %>% dplyr::filter(Cost.E==min(Cost.E))
  }
  if (Checks$N==StartN-NumCheck | Checks$N==StartN+NumCheck) {
    warning("The optimum was at the end of the range. Consider increasing NumCheck or re-running with the output as the new StartN.")
  }
  Checks <- Checks %>% dplyr::mutate(
    SD=sqrt(V.N(s_0, s_1, nu, pi_e, rho_1, rho_2, N, restrict)),
    Ns.E=N*nu,
    M0.E=Ns.E*rho_1*(1-pi_e)*rho_2,
    M1.E=N*rho_1*pi_e*rho_2
  )
  return(Checks %>% dplyr::select(!ends_with(".E"),ends_with(".E")))
}

CheckSet.D <- function(StartN,
                       Delta,
                       theta, s_0, s_1,
                       pi_e, rho_1, rho_2,
                       C_1, C_2, C_e,
                       alpha=0.05,
                       restrict=5,
                       NumCheck=10,
                       nuCheck=2,
                       round=3) {
  Num0 <- NULL
  nu0 <- NULL
  StartD <- Delta.N(theta, s_0, s_1, 
                    nu.V(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e, round),
                    pi_e, rho_1, rho_2, StartN, alpha, restrict)
  
  if (StartD >= Delta) {
    CheckD <- Delta.N(theta, s_0, s_1,
                      nu.V(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e, round)-nuCheck*10^(-1*round),
                      pi_e, rho_1, rho_2, StartN-NumCheck, alpha, restrict)
    
    while (CheckD >= Delta) {
      Num0 <- NumCheck
      nu0 <- nuCheck
      NumCheck <- 2*NumCheck
      nuCheck <- 2*nuCheck
      message(paste0("Increasing Checks Needed: Trying NumCheck=",NumCheck," and nuCheck=",nuCheck))
      CheckD <- Delta.N(theta, s_0, s_1,
                        nu.V(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e, round)-nuCheck*10^(-1*round),
                        pi_e, rho_1, rho_2, StartN-NumCheck, alpha, restrict)
    }
    
    Checks <- cross_join(tibble(N=StartN-(0:NumCheck)),
                         tibble(nu=nu.V(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e, round)-
                                  10^(-1*round)*(0:nuCheck)))
    if (!is.null(Num0)) {
      Checks <- Checks %>% 
        dplyr::filter(N < StartN-Num0 | 
                        nu < nu.V(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e, round)-10^(-1*round)*nu0)
    }
  } else {
    CheckD <- Delta.N(theta, s_0, s_1,
                      nu.V(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e, round)+nuCheck*10^(-1*round),
                      pi_e, rho_1, rho_2, StartN+NumCheck, alpha, restrict)
    
    while (CheckD < Delta) {
      Num0 <- NumCheck
      nu0 <- nuCheck
      NumCheck <- 2*NumCheck
      nuCheck <- 2*nuCheck
      message(paste0("Increasing Checks Needed: Trying NumCheck=",NumCheck," and nuCheck=",nuCheck))
      CheckD <- Delta.N(theta, s_0, s_1,
                        nu.V(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e, round)+nuCheck*10^(-1*round),
                        pi_e, rho_1, rho_2, StartN+NumCheck, alpha, restrict)
    }
    
    Checks <- cross_join(tibble(N=StartN+(0:NumCheck)),
                         tibble(nu=nu.V(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e, round)+
                                  10^(-1*round)*(0:nuCheck)))
    if (!is.null(Num0)) {
      Checks <- Checks %>% 
        dplyr::filter(N > StartN+Num0 | 
                        nu > nu.V(s_0, s_1, pi_e, rho_1, C_1, C_2, C_e, round)+10^(-1*round)*nu0)
    }
  }
  
  Checks$Power <- apply(Checks, 1, function(row) {
    Delta.N(theta, s_0, s_1, row["nu"], pi_e, rho_1, rho_2,
            row["N"], alpha, restrict)
  })
  Checks <- Checks %>% dplyr::filter(Power >= Delta)
  Checks$Cost.E <- apply(Checks, 1, function(row) {
    C.N(C_1, C_2, C_e, row["nu"], pi_e, rho_1, row["N"], quantile=NULL)
  })
  Checks <- Checks %>% dplyr::filter(Cost.E==min(Checks$Cost.E))
  
  if (nrow(Checks) > 1) {
    Checks <- Checks %>% dplyr::filter(Power==max(Power))
  }
  if (Checks$N==StartN-NumCheck | Checks$N==StartN+NumCheck) {
    warning("The optimum was at the end of the range. Consider increasing NumCheck or re-running with the output as the new StartN.")
  }
  Checks <- Checks %>% dplyr::mutate(
    SD=sqrt(V.N(s_0, s_1, nu, pi_e, rho_1, rho_2, N, restrict)),
    Ns.E=N*nu,
    M0.E=Ns.E*rho_1*(1-pi_e)*rho_2,
    M1.E=N*rho_1*pi_e*rho_2
  )
  return(Checks %>% dplyr::select(!ends_with(".E"),ends_with(".E")))
}
