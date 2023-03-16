# KELPER2
# Vital rate regressions
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse)
library(glue); library(brms); library(ipmr)
library(doParallel); library(foreach)
source("code/00_fn.R")

reg.dir <- "out/regr/"
h_i <- read_csv("data/harvest_i.csv")
cores <- 10
s_DD <- 1e-4
s_DD_pow <- 1
nYrs <- 30




# simulation settings -----------------------------------------------------

nSim <- 50
post_draw <- sample.int(4000, nSim)
N.init_rcr <- 5

for(h in 1:nrow(h_i)) {
  
  dir.create(h_i$outDir[h], recursive=T)
  
  sim.i <- list(
    size_rng=log(c(0.1, 50)),
    nYrs=nYrs,
    H_freq=h_i$H_freq[h],
    H_sizeMin=h_i$H_sizeMin[h],
    H_sizeMax=h_i$H_sizeMax[h],
    H_prop=h_i$H_prop[h],
    s_DD=s_DD,
    s_DD_pow=s_DD_pow,
    ndraws=1
  )
  mesh.ls=list(
    L=sim.i$size_rng[1],
    U=sim.i$size_rng[2],
    n_mesh_p=30,
    n_yr=1
  )
  
  saveRDS(sim.i, glue("{h_i$outDir[h]}/sim_i.rds"))
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i=1:nSim,
          .export=c("reg.dir", "h_i", "h", "sim.i"),
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
            H.df <- tibble(year=1:sim.i$nYrs, 
                           OA=0,
                           TURF=0)
            
            N.init <- dnorm(seq(mesh.ls$L+0.125, mesh.ls$U-0.125, length.out=mesh.ls$n_mesh_p),
                            colMeans(posterior_linpred(out.r_z, re.form=NA, newdata=tibble(x=1))), 
                            colMeans(as.matrix(out.r_z, variable="sigma"))
                            ) * N.init_rcr
            
            ipm.OA <- build_IPM(c(ipmPar.ls, mgmt="OA"), mesh.ls)
            ipm.TURF <- build_IPM(c(ipmPar.ls, mgmt="TURF"), mesh.ls)
            N.OA <- N.init
            N.TURF <- N.init
            
            for(j in 1:sim.i$nYrs) {
              if(j > 1 & j %% sim.i$H_freq == 0) {
                size_H <- which(between(size, sim.i$H_sizeMin, sim.i$H_sizeMax))
                H.df$OA[j] <- sum(N.OA[size_H] * wtClass[size_H] * sim.i$H_prop)
                H.df$TURF[j] <- sum(N.TURF[size_H] * wtClass[size_H] * sim.i$H_prop)
                N.OA[size_H] <- N.OA[size_H] * (1-sim.i$H_prop)
                N.TURF[size_H] <- N.TURF[size_H] * (1-sim.i$H_prop)
              }
              OA.j <- ipm.OA %>% 
                define_pop_state(n_z=N.OA) %>% 
                make_ipm(iterate=TRUE, iterations=1)
              TURF.j <- ipm.TURF %>% 
                define_pop_state(n_z=N.TURF) %>% 
                make_ipm(iterate=TRUE, iterations=1)
              
              if(j==1) {
                mesh_info <- int_mesh(OA.j)
                size_ln <- mesh_info$z_1[1:mesh.ls$n_mesh_p]
                size <- exp(size_ln)
                wtClass_ln=c(posterior_epred(out.hf_wt, 
                                             newdata=tibble(logHoldfast=size_ln), 
                                             re_formula=NA, draw_ids=ipmPar.ls$draw))
                wtClass <- exp(wtClass_ln)
              }
              
              lam.df$OA[j] <- c(lambda(OA.j))
              lam.df$TURF[j] <- c(lambda(TURF.j))
              
              pop.ls[[j]] <- bind_rows(
                tibble(year=j,
                       size_ln=size_ln,
                       N=c(OA.j$pop_state$n_z[,2]),
                       Management="OA"),
                tibble(year=j,
                       size_ln=size_ln,
                       N=c(TURF.j$pop_state$n_z[,2]),
                       Management="TURF")
              )
              N.OA <- c(OA.j$pop_state$n_z[,2])
              N.TURF <- c(TURF.j$pop_state$n_z[,2])
            }
            
            lam.df %>%
              pivot_longer(2:3, names_to="Management", values_to="lambda") %>%
              mutate(sim=i, 
                     scenario=h) %>%
              saveRDS(glue("{h_i$outDir[h]}/lambda_{i}.rds"))
            H.df %>%
              pivot_longer(2:3, names_to="Management", values_to="H_grams") %>%
              mutate(sim=i, 
                     scenario=h) %>%
              saveRDS(glue("{h_i$outDir[h]}/H_{i}.rds"))
            pop.df <- pop.ls %>% 
              do.call('rbind', .) %>%
              mutate(sim=i, 
                     scenario=h,
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
            
            saveRDS(pop.df, glue("{h_i$outDir[h]}/pop_{i}.rds"))
            
            paste("Finished", i)
          }
  stopCluster(cl)
}





