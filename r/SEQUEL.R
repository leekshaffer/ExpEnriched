source("R/Simulate.R")
options(scipen=8)
options(warn=1)

## Exposure Prevalence:
pi_e <- 0.05

alpha <- 0.05

## Rough Cost Numbers:
C_1 <- 40
C_2 <- 75*3+30*2
# C_2.alt <- 125
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

## Fixing N=2000 at nu=15%:
Design.1.alt <- tibble(N=2000, nu=0.15) %>%
  dplyr::mutate(Power=Delta.N(theta, s_0, s_1,
                              nu, pi_e, rho_1, rho_2, N, alpha=0.05, restrict=5),
                SD=sqrt(V.N(s_0,s_1,nu,pi_e,rho_1,rho_2,N,
                             restrict=5)),
                Cost.E=C.N(C_1, C_2, C_e, nu, pi_e, rho_1, N),
                Ns.E=N*nu,
                M0.E=Ns.E*rho_1*(1-pi_e)*rho_2,
                M1.E=N*rho_1*pi_e*rho_2)
Design.1.alt

set.seed(20261023)

Sim.1 <- SimMany(1000, N=Design.1$N, nu=Design.1$nu,
                 theta, s_0, s_1, pi_e, rho_1, rho_2, C_1, C_2, C_e)
Sim.1.alt <- SimMany(1000, N=Design.1.alt$N, nu=Design.1.alt$nu,
                 theta, s_0, s_1, pi_e, rho_1, rho_2, C_1, C_2=C_2, C_e)

SRes.1 <- c(apply(Sim.1 %>% 
                  dplyr::select(Est, Diff, AbsDiff, SE, CI.Incl, Sig, Cost),
                2, mean),
            SD_Est=sd(Sim.1$Est))
SRes.1.alt <- c(apply(Sim.1.alt %>% 
                      dplyr::select(Est, Diff, AbsDiff, SE, CI.Incl, Sig, Cost),
                    2, mean),
                SD_Est=sd(Sim.1.alt$Est))

Designs.1 <- tibble(Design=c("1","1.alt")) %>% bind_cols(
  Design.1 %>% bind_rows(Design.1.alt))
SRes.1 <- tibble(Design=c("1","1.alt")) %>% bind_cols(
  SRes.1 %>% bind_rows(SRes.1.alt))

save(list=c("Designs.1", "SRes.1",
            "Sim.1","Sim.1.alt"),
     file="out/Aim1_Designs.Rda")


## Aim 2, Binary Outcome:
pi_0 <- 0.2
pi_1 <- 0.35
s_0 <- sqrt(pi_0*(1-pi_0))
s_1 <- sqrt(pi_1*(1-pi_1))
theta <- abs(pi_1-pi_0)

Design.2 <- OptimalBudget(Delta=0.8,
                          theta, s_0, s_1,
                          pi_e, rho_1, rho_2,
                          C_1, C_2, C_e,
                          alpha=0.05,
                          restrict=5,
                          NumCheck=0,
                          nuCheck=0,
                          round=3)

Designs.2 <- Designs.1 %>%
  bind_rows(Design.2 %>% mutate(Design="2")) %>%
  dplyr::select(Design,N,nu)

Designs.2$Power <- apply(as.matrix(Designs.2 %>% dplyr::select(N,nu)), 
                         1, function(row) {
  Delta.N(theta, s_0, s_1, row["nu"], pi_e, rho_1, rho_2, row["N"],
          alpha, restrict=5)
})
Designs.2$SD <- apply(as.matrix(Designs.2 %>% dplyr::select(N,nu)), 
                      1, function(row) {
                        sqrt(V.N(s_0, s_1, row["nu"], pi_e, rho_1, rho_2, row["N"], restrict=5))
                      })
Designs.2$Cost.E <- apply(as.matrix(Designs.2 %>% dplyr::select(N,nu)), 
                          1, function(row) {
                            C.N(C_1, C_2, C_e, row["nu"], pi_e, rho_1, row["N"], quantile=NULL)
                          })
Designs.2 <- Designs.2 %>%
  dplyr::mutate(Ns.E=N*nu,
                M0.E=Ns.E*rho_1*(1-pi_e)*rho_2,
                M1.E=N*rho_1*pi_e*rho_2)

set.seed(10232026)

Sim.2.1 <- SimManyBin(NoSims=1000,
                      N=unname(Designs.2 %>% dplyr::filter(Design=="1") %>% pull(N)), 
                      nu=unname(Designs.2 %>% dplyr::filter(Design=="1") %>% pull(nu)),
                      pi_0, pi_1, pi_e, rho_1, rho_2, C_1, C_2, C_e, alpha)
SRes.2.1 <- c(apply(Sim.2.1 %>% 
                    dplyr::select(ends_with(c("Est", "Diff", "AbsDiff", "SE", "CI.Incl", "Sig", "Cost"))),
                  2, mean),
            RD.SD_Est=sd(Sim.2.1$RD.Est),
            OR.SD_Est=sd(Sim.2.1$OR.Est))

Sim.2.1.alt <- SimManyBin(NoSims=1000,
                      N=unname(Designs.2 %>% dplyr::filter(Design=="1.alt") %>% pull(N)), 
                      nu=unname(Designs.2 %>% dplyr::filter(Design=="1.alt") %>% pull(nu)),
                      pi_0, pi_1, pi_e, rho_1, rho_2, C_1, C_2, C_e, alpha)
SRes.2.1.alt <- c(apply(Sim.2.1.alt %>% 
                      dplyr::select(ends_with(c("Est", "Diff", "AbsDiff", "SE", "CI.Incl", "Sig", "Cost"))),
                    2, mean),
              RD.SD_Est=sd(Sim.2.1.alt$RD.Est),
              OR.SD_Est=sd(Sim.2.1.alt$OR.Est))

Sim.2.2 <- SimManyBin(NoSims=1000,
                      N=unname(Designs.2 %>% dplyr::filter(Design=="2") %>% pull(N)), 
                      nu=unname(Designs.2 %>% dplyr::filter(Design=="2") %>% pull(nu)),
                      pi_0, pi_1, pi_e, rho_1, rho_2, C_1, C_2, C_e, alpha)
SRes.2.2 <- c(apply(Sim.2.2 %>% 
                      dplyr::select(ends_with(c("Est", "Diff", "AbsDiff", "SE", "CI.Incl", "Sig", "Cost"))),
                    2, mean),
              RD.SD_Est=sd(Sim.2.2$RD.Est),
              OR.SD_Est=sd(Sim.2.2$OR.Est))

SRes.2 <- tibble(Design=c("1","1.alt","2")) %>% bind_cols(
  SRes.2.1 %>% bind_rows(SRes.2.1.alt, SRes.2.2))

save(list=c("Designs.2", "SRes.2",
            "Sim.2.1","Sim.2.1.alt","Sim.2.2"),
     file="out/Aim2_Designs.Rda")
