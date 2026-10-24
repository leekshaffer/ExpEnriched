source("R/Simulate.R")
source("R/Optimize.R")
options(scipen=8)
options(warn=1)
set.seed(20261023)

## Exposure Prevalence:
pi_e <- 0.05

alpha <- 0.05

## Rough Cost Numbers:
C_1 <- 40
C_2 <- 75*3+30*2
C_e <- 175

## Rough LTFU Numbers:
rho_1 <- 0.95
rho_2 <- 0.96

## Aim 1, Continuous Outcome:
s_0 <- 1
s_1 <- 1
theta <- 0.4

Design.1 <- OptimalBudget(Delta=0.9,
              theta, s_0, s_1,
              pi_e, rho_1, rho_2,
              C_1, C_2, C_e,
              alpha=0.05,
              restrict=5,
              NumCheck=10,
              nuCheck=2,
              round=3)
Design.1

Designs.Aim1 <- tibble(Design=1:3,
                       N=c(Design.1$N, 2000, 2100),
                       nu=c(Design.1$nu, .15, .12))

Designs.Aim1 <- AnalyzeDesigns(Designs.Aim1,
                               theta, s_0, s_1,
                               pi_e, rho_1, rho_2,  
                               C_1, C_2, C_e,
                               alpha=0.05, restrict=5,
                               NoSims=2000)

Designs.Aim1 %>% dplyr::select(ends_with(c("Design","N","nu","Cvg","Power")))  

## Aim 2, Binary Outcome:
pi_0 <- 0.2
pi_1 <- 0.35
s_0.bin <- sqrt(pi_0*(1-pi_0))
s_1.bin <- sqrt(pi_1*(1-pi_1))
theta.bin <- abs(pi_1-pi_0)
Delta.bin <- 0.8

Design.4 <- OptimalBudget(Delta.bin,
                          theta.bin, s_0.bin, s_1.bin,
                          pi_e, rho_1, rho_2,
                          C_1, C_2, C_e,
                          alpha=0.05,
                          restrict=5,
                          NumCheck=5,
                          nuCheck=2,
                          round=3)

Designs.Aim2 <- Designs.Aim1 %>% 
  dplyr::select(Design,N,nu) %>%
  add_row(Design=4, N=Design.4$N, nu=Design.4$nu)

Designs.Aim2 <- AnalyzeDesignsBin(Designs.Aim2,
                                  pi_0, pi_1,
                                  pi_e, rho_1, rho_2,  
                                  C_1, C_2, C_e,
                                  alpha=0.05, restrict=5,
                                  NoSims=2000)

Designs.Aim2 %>% dplyr::select(ends_with(c("Design","N","nu","Cvg","Power")))                               

save(list=c("Designs.Aim1", "Designs.Aim2"),
     file="out/Design_Results.Rda")