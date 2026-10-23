source("R/Simulate.R")
options(scipen=8)

## Brought Over:
sim <- SimMany(2000, 1925, .218, .4, 1, 1, .05, 1, .96, 120, 90, 180)
apply(sim %>% dplyr::select(Est, Diff, AbsDiff, SE, CI.Incl, Sig, Cost), 2, mean)
sd(sim$Est)

simComp <- SimManyComp(2000, 1925, .4, 1, 1, .05, 1, .96, 120, 90, 180)
apply(simComp %>% dplyr::select(Est, Diff, AbsDiff, SE, CI.Incl, Sig, Cost), 2, mean)
sd(simComp$Est)



## Numbers to use:
OptimalBudget(Delta=0.8,
              pi_0=0.2, pi_1=0.35,
              rho_1=0.92, pi_e=0.05, rho_2=0.964,
              C_1=90, C_2=90, C_e=160,
              alpha=0.05,
              NumCheck=20,
              quantile=NULL)



## Optimizing:

OptimalBudget(Delta=0.8,
              pi_0=0.1, pi_1=0.25,
              rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
              C_1=1000, C_2=500, C_e=200,
              alpha=0.05,
              NumCheck=0,
              quantile=NULL)

OptimalBudget(Delta=0.8,
              pi_0=0.2, pi_1=0.35,
              rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
              C_1=1000, C_2=500, C_e=200,
              alpha=0.05,
              NumCheck=20,
              quantile=NULL)

OptimalBudget(Delta=0.8,
              pi_0=0.1, pi_1=0.25,
              rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
              C_1=1000, C_2=500, C_e=200,
              alpha=0.05,
              NumCheck=20,
              quantile=0.75)

OptimalPower(B=1900000,
             pi_0=0.1,pi_1=0.25,
             rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
             C_1=1000, C_2=500, C_e=200,
             alpha=0.05,
             NumCheck=5,
             quantile=NULL)

OptimalIntBudget(Delta=0.8,
                 pi_0=0.1,pi_1=0.25,
                 rho_2=0.96^2, C_2=500,
                 N_0=369, N_1=83,
                 alpha=0.05,
                 NumCheck=5)

OptimalIntBudget(Delta=0.8,
                 pi_0=0.1,pi_1=0.25,
                 rho_2=0.96^2, C_2=500,
                 N_0=369, N_1=79,
                 alpha=0.05,
                 NumCheck=5)

### Simulations:
set.seed(20261022)

Sim80Int <- SimManyInt(NoSims=2000,
                       N0=237, N1=83,
                       pi_0=0.1, pi_1=0.25,
                       rho_2=0.96^2,
                       alpha=0.05)
sum(Sim80Int$Failures)
apply(Sim80Int %>% dplyr::select(ends_with("N1"),ends_with(".Est"),ends_with(".SE"),ends_with(".CI.inc")), 2, mean)
apply(Sim80Int %>% dplyr::select(ends_with(".p")), 2, function(x) mean(x < 0.05))


Sim80 <- SimMany(NoSims=2000,
                 N=1651, N0=369,
                 pi_0=0.1, pi_1=0.25,
                 rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
                 C_1=1000, C_2=500, C_e=200,
                 alpha=0.05)
sum(Sim80$Failures)
apply(Sim80 %>% dplyr::select(ends_with("N1"),ends_with(".Est"),ends_with(".SE"),ends_with(".CI.inc")), 2, mean)
apply(Sim80 %>% dplyr::select(ends_with(".p")), 2, function(x) mean(x < 0.05))
mean(Sim80$C.N1)
save(Sim80, file="out/Sim80.Rda")


Var2M <- OptimalPower(2000000, .1, .25, .96, .05, .96^2, 1000, 500, 200)
Sim2M <- SimMany(NoSims=2000,
                 N=1750, N0=382,
                 pi_0=0.1, pi_1=0.25,
                 rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
                 C_1=1000, C_2=500, C_e=200,
                 alpha=0.05)
sum(Sim2M$Failures)
apply(Sim2M %>% dplyr::select(ends_with("N1"),ends_with(".Est"),ends_with(".SE"),ends_with(".CI.inc")), 2, mean)
apply(Sim2M %>% dplyr::select(ends_with(".p")), 2, function(x) mean(x < 0.05))
mean(Sim2M$C.N1)
save(Sim2M, file="out/Sim80.Rda")


Sim80Adj <- SimManyAdj(NoSims=1000,
                       N=1651, N0_Max=500, Delta.Adj=.8,
                       pi_0=0.1, pi_1=0.25,
                       rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
                       C_1=1000, C_2=500, C_e=200,
                       alpha=0.05)
sum(Sim80Adj$Failures)
apply(Sim80Adj %>% dplyr::select(ends_with("N1"),ends_with(".Est"),ends_with(".SE"),ends_with(".CI.inc")), 2, mean)
apply(Sim80Adj %>% dplyr::select(ends_with(".p")), 2, function(x) mean(x < 0.05))
save(Sim80Adj, file="out/Sim80Adj.Rda")

## Comparisons:

## At 80%:
Sim80Comp <- SimManyComp(NoSims=2000,
                         N=ceiling(N.star.Delta.Comp(0.1, 0.25, 0.96, 0.05, 0.96^2, 0.8)),
                         pi_0=0.1, pi_1=0.25,
                         rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
                         C_1=1000, C_2=500, C_e=200,
                         alpha=0.05)
sum(Sim80Comp$Failures)
apply(Sim80Comp %>% dplyr::select(ends_with("N1"),ends_with(".Est"),ends_with(".SE"),ends_with(".CI.inc")), 2, mean)
apply(Sim80Comp %>% dplyr::select(ends_with(".p")), 2, function(x) mean(x < 0.05))
mean(Sim80Comp$C.Full)

save(Sim80Comp, file="out/Sim80Comp.Rda")

## At Budget = 2,000,000 (just above avg budget from exposure-enriched 80% power)
Sim2MComp <- SimManyComp(NoSims=2000,
                         N=floor(BToN.Comp(2000000, 0.96, 0.05, 
                                           1000, 500, 200)),
                         pi_0=0.1, pi_1=0.25,
                         rho_1=0.96, pi_e=0.05, rho_2=0.96^2,
                         C_1=1000, C_2=500, C_e=200,
                         alpha=0.05)
sum(Sim2MComp$Failures)
apply(Sim2MComp %>% dplyr::select(ends_with("N1"),ends_with(".Est"),ends_with(".SE"),ends_with(".CI.inc")), 2, mean)
apply(Sim2MComp %>% dplyr::select(ends_with(".p")), 2, function(x) mean(x < 0.05))
mean(Sim2MComp$C.Full)

save(Sim2MComp, file="out/Sim2MComp.Rda")
