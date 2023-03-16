# KELPER2
# Vital rate regressions
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse)
library(glue); library(brms); library(ipmr)
library(doParallel); library(foreach)
source("code/00_fn.R")

reg.dir <- "out/regr/"
out.dir <- "out/0_baseline/"
cores <- 10
s_DD <- 0.007

dir.create(out.dir, recursive=T)



# simulation settings -----------------------------------------------------

nSim <- 1000
post_draw <- sample.int(4000, nSim)
sim.i <- list(
  size_rng=log(c(0.1, 50)),
  nYrs=50,
  harvestYr=10,
  s_DD=s_DD,
  s_DD_pow=0.5,
  ndraws=1
)
mesh.ls=list(
  L=sim.i$size_rng[1],
  U=sim.i$size_rng[2],
  n_mesh_p=30,
  n_yr=1
)

saveRDS(sim.i, glue("{out.dir}/sim_i.rds"))

cl <- makeCluster(cores)
registerDoParallel(cl)
foreach(i=1:nSim,
        .export=c("reg.dir", "out.dir", "sim.i"),
        .packages=c("tidyverse", "glue", "brms", "ipmr", "doParallel"), 
        .combine="c") %dopar% {
          
          source("code/00_fn.R")
          out.g <- readRDS(glue("{reg.dir}opt_g.rds"))
          out.s <- readRDS(glue("{reg.dir}opt_s.rds"))
          out.r <- readRDS(glue("{reg.dir}opt_r.rds"))
          out.r_z <- readRDS(glue("{reg.dir}opt_r_z.rds"))
          out.hf_wt <- readRDS(glue("{reg.dir}opt_allom_hw.rds"))
          
          ipmPar.ls <- list(
            out_g=out.g,
            out_s=out.s,
            s_DD=sim.i$s_DD,
            s_DD_pow=sim.i$s_DD_pow,
            out_r=out.r,
            out_r_z=out.r_z,
            draw=post_draw[i],
            adultThresh=10
          )
          
          lam.df <- tibble(year=1:sim.i$nYrs, 
                           OA=NA,
                           TURF=NA)
          pop.ls <- vector("list", length=sim.i$nYrs)
          
          N.init <- c(1, rep(0, mesh.ls$n_mesh_p-1))
          ipm.OA <- build_IPM(c(ipmPar.ls, mgmt="OA"), mesh.ls)
          ipm.TURF <- build_IPM(c(ipmPar.ls, mgmt="TURF"), mesh.ls)
          N.OA <- N.init
          N.TURF <- N.init
          
          for(j in 1:sim.i$nYrs) {
            if(j == sim.i$harvestYr) {
              N.OA <- N.OA * 0.05
              N.TURF <- N.TURF * 0.05
            }
            OA.j <- ipm.OA %>% 
              define_pop_state(n_z=N.OA) %>% 
              make_ipm(iterate=TRUE, iterations=1)
            TURF.j <- ipm.TURF %>% 
              define_pop_state(n_z=N.TURF) %>% 
              make_ipm(iterate=TRUE, iterations=1)
            
            mesh_info <- int_mesh(OA.j)
            
            lam.df$OA[j] <- c(lambda(OA.j))
            lam.df$TURF[j] <- c(lambda(TURF.j))
            
            pop.ls[[j]] <- bind_rows(
              tibble(year=j,
                     size_ln=mesh_info$z_1[1:mesh.ls$n_mesh_p],
                     N=c(OA.j$pop_state$n_z[,2]),
                     Management="OA"),
              tibble(year=j,
                     size_ln=mesh_info$z_1[1:mesh.ls$n_mesh_p],
                     N=c(TURF.j$pop_state$n_z[,2]),
                     Management="TURF")
            )
            N.OA <- c(OA.j$pop_state$n_z[,2])
            N.TURF <- c(TURF.j$pop_state$n_z[,2])
          }
          
          lam.df %>%
            pivot_longer(2:3, names_to="Management", values_to="lambda") %>%
            mutate(sim=i,
                   s_DD=sim.i$s_DD) %>%
            saveRDS(glue("{out.dir}/lambda_{i}.rds"))
          pop.df <- pop.ls %>% 
            do.call('rbind', .) %>%
            mutate(sim=i,
                   s_DD=sim.i$s_DD,
                   size=exp(size_ln),
                   stage=if_else(size>10, "Adult", "Juvenile"),
                   sizeClass=case_when(size < 10 ~ "0-10",
                                       between(size, 10, 20) ~ "10-20",
                                       between(size, 20, 30) ~ "20-30",
                                       between(size, 30, 40) ~ "30-40",
                                       between(size, 40, 50) ~ "40-50",
                                       size > 50 ~ "50+"))
          pop.df <- pop.df %>%
            mutate(weightClass_ln=c(posterior_epred(out.hf_wt, 
                                                    newdata=pop.df %>% rename(logHoldfast=size_ln), 
                                                    re_formula=NA, draw_ids=ipmPar.ls$draw)))
          
          saveRDS(pop.df, glue("{out.dir}/pop_{i}.rds"))
          
          paste("Finished", i)
        }
stopCluster(cl)



lam.df <- map_dfr(dir(out.dir, "lambda", recursive=T, full.names=T), readRDS)
pop.df <- map_dfr(dir(out.dir, "pop", recursive=T, full.names=T), readRDS)


ggplot(lam.df, aes(year, lambda, colour=Management, group=paste(Management, sim))) +
  geom_line() +
  facet_wrap(~s_DD)
lam.df %>% group_by(year, Management, s_DD) %>%
  summarise(lambda_mn=exp(mean(log(lambda))),
            lambda_q10=quantile(lambda, 0.05),
            lambda_q90=quantile(lambda, 0.95)) %>%
  ggplot(aes(year, lambda_mn, ymin=lambda_q10, ymax=lambda_q90, colour=Management, fill=Management)) +
  geom_ribbon(alpha=0.5, colour=NA) +
  geom_line() +
  facet_wrap(~s_DD)

pop.df %>%
  filter(year >= 10) %>%
  group_by(year, Management, stage, sim) %>%
  summarise(N=sum(N)) %>%
  # mutate(N=log(N)) %>%
  group_by(year, Management, stage) %>%
  summarise(N_mn=median(N),
            N_lo1=HDInterval::hdi(N, 0.8)[1],
            N_hi1=HDInterval::hdi(N, 0.8)[2],
            N_lo2=HDInterval::hdi(N, 0.95)[1],
            N_hi2=HDInterval::hdi(N, 0.95)[2]) %>%
  ggplot(aes(year, N_mn, colour=Management, fill=Management)) +
  geom_ribbon(aes(ymin=N_lo1, ymax=N_hi1), alpha=0.25, colour=NA) +
  geom_ribbon(aes(ymin=N_lo2, ymax=N_hi2), alpha=0.25, colour=NA) +
  geom_line() +
  scale_y_continuous(trans="log1p") +
  facet_grid(stage~., scales="free_y")

pop.df %>%
  filter(year >= 10) %>%
  group_by(year, Management, stage, sim) %>%
  summarise(weight=sum(N * exp(weightClass_ln))) %>%
  group_by(year, Management, stage) %>%
  summarise(weight_mn=median(weight),
            weight_lo1=HDInterval::hdi(weight, 0.8)[1],
            weight_hi1=HDInterval::hdi(weight, 0.8)[2],
            weight_lo2=HDInterval::hdi(weight, 0.90)[1],
            weight_hi2=HDInterval::hdi(weight, 0.90)[2]) %>%
  ggplot(aes(year, weight_mn, colour=stage, fill=stage)) +
  geom_ribbon(aes(ymin=weight_lo1, ymax=weight_hi1), alpha=0.25, colour=NA) +
  geom_ribbon(aes(ymin=weight_lo2, ymax=weight_hi2), alpha=0.25, colour=NA) +
  geom_line() +
  # scale_y_continuous(trans="log") +
  facet_grid(Management~.)

