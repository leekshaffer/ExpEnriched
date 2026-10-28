library(tidyverse)
source("r/BaseFunctions.R")

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

SimIntBin <- function(
    M0, M1,
    pi_0, pi_1,
    pi_e, rho_1, rho_2,
    C_1, C_2, C_e,
    alpha=0.05) {
  Y0 <- rbinom(n=M0, size=1, prob=pi_0)
  Y1 <- rbinom(n=M1, size=1, prob=pi_1)
  
  dat <- tibble(A=rep(c(0,1), times=c(M0,M1)),
                Y=c(Y0,Y1))
  
  lm <- tryCatch({
    suppressWarnings(glm(Y~A, family=binomial(link="identity"), data=dat))
  }, error=function (e) {
    message("Error in One Simulation Identity GLM. Returning NA for that simulation.")
    NULL
  })
  
  if (is.null(lm)) {
    RD.Est <- NA_real_
    RD.SE <- NA_real_
    RD.p <- NA_real_
    RD.CIL <- NA_real_
    RD.CIU <- NA_real_
  } else {
    RD.Est <- unname(coef(lm)["A"])
    RD.SE <- summary(lm)$coefficients["A","Std. Error"]
    RD.p <- summary(lm)$coefficients["A","Pr(>|z|)"]
    ci <- tryCatch({
      suppressWarnings(suppressMessages(confint(lm, parm="A", level=1-alpha)))
    }, error=function(e) {
      message("Error in One Simulation Identity GLM CI. Returning NA for that simulation.")
      NULL
    })
    if (is.null(ci)) {
      RD.CIL <- NA_real_
      RD.CIU <- NA_real_
    } else {
      RD.CIL <- unname(ci[1])
      RD.CIU <- unname(ci[2])
    }
  }
  
  glm <- tryCatch({
    suppressWarnings(glm(Y~A, family=binomial(link="logit"), data=dat))
  }, error=function(e) {
    message("Error in One Simulation Logit GLM. Returning NA for that simulation.")
    NULL
    }
  )
  
  if (is.null(glm)) {
    OR.Est <- NA_real_
    LOR.SE <- NA_real_
    OR.p <- NA_real_
    OR.CIL <- NA_real_
    OR.CIU <- NA_real_
  } else {
    OR.Est <- exp(unname(coef(glm)["A"]))
    LOR.SE <- summary(glm)$coefficients["A","Std. Error"]
    OR.p <- summary(glm)$coefficients["A","Pr(>|z|)"]
    gci <- tryCatch({
      suppressWarnings(suppressMessages(confint(glm, parm="A", level=1-alpha)))
    }, error=function(e) {
      message("Error in One Simulation Logit GLM CI. Returning NA for that simulation.")
      NULL
      })
    
    if (is.null(gci)) {
      OR.CIL <- NA_real_
      OR.CIU <- NA_real_
    } else {
      OR.CIL <- exp(unname(gci[1]))
      OR.CIU <- exp(unname(gci[2]))
    }
  }
  
  return(c(Y0bar=mean(Y0), Y1bar=mean(Y1),
           Y0sd=sd(Y0), Y1sd=sd(Y1),
           RD.Est=RD.Est,
           RD.SE=RD.SE,
           RD.p=RD.p,
           RD.CIL=RD.CIL,
           RD.CIU=RD.CIU,
           OR.Est=OR.Est,
           LOR.SE=LOR.SE,
           OR.p=OR.p,
           OR.CIL=OR.CIL,
           OR.CIU=OR.CIU))
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

SimOneBin <- function(
    N, nu, 
    pi_0, pi_1,
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
    return(SimOneBin(N, nu,
                     pi_0, pi_1,
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
    return(SimOneBin(N, nu,
                     pi_0, pi_1,
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
    return(SimOneBin(N, nu,
                     pi_0, pi_1,
                  pi_e, rho_1, rho_2,
                  C_1, C_2, C_e,
                  alpha,
                  FailStart=Failures))
  }
  
  return(c(N0=N0, N1=N1, M0=M0, M1=M1,
           Cost=C_1*N + C_2*(N0+N1) + C_e*N1,
           SimIntBin(M0, M1,
                     pi_0, pi_1,
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
    dplyr::mutate(Bias=Est-theta,
                  AbsBias=abs(Est-theta),
                  CI.Cvg=(CIL <= theta & theta <= CIU),
                  Sig=(p < alpha))
}

SimManyBin <- function(NoSims,
                    N, nu, 
                    pi_0, pi_1,
                    pi_e, rho_1, rho_2,
                    C_1, C_2, C_e,
                    alpha=0.05) {
  RD <- pi_1 - pi_0
  OR <- (pi_1/(1-pi_1))/(pi_0/(1-pi_0))
  as_tibble(t(replicate(n=NoSims, expr=SimOneBin(N, nu, 
                                              pi_0, pi_1,
                                              pi_e, rho_1, rho_2,
                                              C_1, C_2, C_e,
                                              alpha=0.05,
                                              FailStart=0)))) %>%
    dplyr::mutate(RD.Bias=RD.Est-RD,
                  RD.AbsBias=abs(RD.Est-RD),
                  RD.CI.Cvg=(RD.CIL <= RD & RD <= RD.CIU),
                  RD.Sig=(RD.p < alpha),
                  OR.Bias=OR.Est-OR,
                  OR.AbsBias=abs(OR.Est-OR),
                  OR.CI.Cvg=(OR.CIL <= OR & OR <= OR.CIU),
                  OR.Sig=(OR.p < alpha))
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
    dplyr::mutate(Bias=Est-theta,
                  AbsBias=abs(Est-theta),
                  CI.Cvg=(CIL <= theta & theta <= CIU),
                  Sig=(p < alpha))
}
