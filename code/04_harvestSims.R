# KELPER2
# Vital rate regressions
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse)
library(glue); library(brms); library(ipmr)
library(doParallel); library(foreach)
source("code/00_fn.R")

reg.dir <- "out/regr/"
h_i <- read_csv("data/harvest_i.csv")[1:3,]
N.init_rcr <- 5



# simulation settings -----------------------------------------------------

cores <- 70
s_DD <- -0.000789
s_DD_pow <- 1
rcr_max <- 100
nYrs <- 30
nSim <- 1000



# run simulations ---------------------------------------------------------

ipm.scenarios <- expand_grid(mgmt=c("OA", "TURF"), upwell=c("N", "Y"))
post_draw <- sample.int(4000, nSim)

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
    ndraws=1,
    rcr_max=rcr_max
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
            out.hf_wt <- readRDS(glue("{reg.dir}opt_allom_hw.rds"))
            
            ipmPar.ls <- list(
              out_g=readRDS(glue("{reg.dir}opt_g.rds")),
              out_s=readRDS(glue("{reg.dir}opt_s.rds")),
              out_r=readRDS(glue("{reg.dir}opt_r.rds")),
              out_r_z=readRDS(glue("{reg.dir}opt_r_z.rds")),
              s_DD=sim.i$s_DD,
              s_DD_pow=sim.i$s_DD_pow,
              draw=post_draw[i],
              adultThresh=10,
              rcr_max=sim.i$rcr_max
            )
            
            pop.ls <- map(1:nrow(ipm.scenarios), ~vector("list", length=sim.i$nYrs))
            H.ls <- vector("list", nrow(ipm.scenarios))
            
            N.init <- dnorm(seq(mesh.ls$L+0.125, mesh.ls$U-0.125, length.out=mesh.ls$n_mesh_p),
                            colMeans(posterior_linpred(ipmPar.ls$out_r_z, re.form=NA, newdata=tibble(x=1))), 
                            colMeans(as.matrix(ipmPar.ls$out_r_z, variable="sigma"))
            ) 
            N.init <- N.init/sum(N.init) * N.init_rcr
            
            # ipm scenarios -----------------------------------------------
            for(j in 1:nrow(ipm.scenarios)) {
              pars.j <- c(ipmPar.ls,
                          mgmt=ipm.scenarios$mgmt[j],
                          upwell=ipm.scenarios$upwell[j])
              ipm_init <- build_IPM2(pars.j, mesh.ls)
              N <- N.init
              H <- rep(0, sim.i$nYrs)
              
              for(k in 1:sim.i$nYrs) {
                if(k > 1 & k %% sim.i$H_freq == 0) {
                  size_H <- which(between(size, sim.i$H_sizeMin, sim.i$H_sizeMax))
                  H[k] <- sum(N[size_H] * wtClass[size_H] * sim.i$H_prop)
                  N[size_H] <- N[size_H] * (1-sim.i$H_prop)
                }
                
                # Implement growth and survival
                ipm.j <- ipm_init %>%
                  define_pop_state(n_z=N) %>%
                  make_ipm(iterate=TRUE, iterations=1)
                
                if(k==1) {
                  mesh_info <- int_mesh(ipm.j)
                  size_ln <- mesh_info$z_1[1:mesh.ls$n_mesh_p]
                  size <- exp(size_ln)
                  juv_i <- which(size < pars.j$adultThresh)
                  adult_i <- which(size >= pars.j$adultThresh)
                  wtClass_ln=c(posterior_epred(out.hf_wt, 
                                               newdata=tibble(logHoldfast=size_ln), 
                                               re_formula=NA, draw_ids=ipmPar.ls$draw))
                  wtClass <- exp(wtClass_ln)
                }
                
                # Ipmlement recruitment
                ipm.j$pop_state$n_z[,2] <- calc_new_recruits(ipm.j, pars.j, size_ln, adult_i)
                
                pop.ls[[j]][[k]] <- tibble(year=k,
                                           size_ln=mesh_info$z_1[1:mesh.ls$n_mesh_p],
                                           N=c(ipm.j$pop_state$n_z[,2]),
                                           Management=ipm.scenarios$mgmt[j],
                                           Upwelling=ipm.scenarios$upwell[j])
                N <- c(ipm.j$pop_state$n_z[,2])
              }
              H.ls[[j]] <- tibble(year=1:sim.i$nYrs,
                                  H=H,
                                  Management=ipm.scenarios$mgmt[j],
                                  Upwelling=ipm.scenarios$upwell[j])
            }
            
            pop.df <- map_dfr(pop.ls, 
                              ~do.call('rbind', .x) %>%
                                mutate(sim=i,
                                       scenario=h,
                                       s_DD=sim.i$s_DD, 
                                       s_DD_pow=sim.i$s_DD_pow,
                                       size=exp(size_ln),
                                       stage=if_else(size>10, "Adult", "Juvenile"),
                                       sizeClass=case_when(size < 10 ~ "0-10",
                                                           between(size, 10, 20) ~ "10-20",
                                                           between(size, 20, 30) ~ "20-30",
                                                           between(size, 30, 40) ~ "30-40",
                                                           between(size, 40, 50) ~ "40-50",
                                                           size > 50 ~ "50+")))
            pop.df %>%
              mutate(weightClass_ln=c(posterior_epred(out.hf_wt, 
                                                      newdata=pop.df %>% rename(logHoldfast=size_ln), 
                                                      re_formula=NA, draw_ids=ipmPar.ls$draw))) %>%
              saveRDS(glue("{h_i$outDir[h]}/pop_{i}.rds"))
            do.call('rbind', H.ls) %>%
              mutate(H_cumul=cumsum(H), 
                     scenario=h,
                     sim=i,
                     s_DD=sim.i$s_DD, 
                     s_DD_pow=sim.i$s_DD_pow) %>%
            saveRDS(glue("{h_i$outDir[h]}/harvest_{i}.rds"))
       
            paste("Finished", i)
          }
  stopCluster(cl)
}





