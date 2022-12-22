# KELPER2
# Vital rate regressions
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse)
library(glue); library(brms); library(ipmr)
library(doParallel); library(foreach)
source("code/00_fn.R")

reg.dir <- "temp/"
s_DD.seq <- seq(0, -0.015, length.out=31)

for(s in seq_along(s_DD.seq)) {
  
  out.dir <- glue("out/dd_tune/sdd-{str_pad(s, 2, 'left', '0')}/")
  dir.create(out.dir, recursive=T)
  
  nSim=200
  dd.i <- list(
    size_rng=log(c(0.1, 50)),
    nYrs=40,
    harvestYr=10,
    s_DD=s_DD.seq[s],
    ndraws=1
  )
  mesh.ls=list(
    L=dd.i$size_rng[1] * 1.2,
    U=dd.i$size_rng[2] * 1.2,
    n_mesh_p=30,
    n_yr=1
  )
  
  saveRDS(dd.i, glue("{out.dir}/dd_i.rds"))
  
  cl <- makeCluster(7)
  registerDoParallel(cl)
  foreach(i=1:nSim,
          .export=c("reg.dir", "out.dir", "dd.i"),
          .packages=c("tidyverse", "glue", "brms", "ipmr", "doParallel"), 
          .combine="c") %dopar% {
            
            source("code/00_fn.R")
            out.g <- readRDS(glue("{reg.dir}out_g.rds"))
            out.s <- readRDS(glue("{reg.dir}out_s.rds"))
            out.r <- readRDS(glue("{reg.dir}out_r.rds"))
            out.r_z <- readRDS(glue("{reg.dir}out_r_z.rds"))
            
            ipmPar.ls <- list(
              out_g=out.g,
              out_s=out.s,
              s_DD=dd.i$s_DD,
              out_r=out.r,
              out_r_z=out.r_z,
              draw=sample.int(4000, 1),
              adultThresh=10
            )
            
            lam.df <- tibble(year=1:dd.i$nYrs, 
                             OA=NA,
                             TURF=NA)
            pop.ls <- vector("list", length=dd.i$nYrs)
            
            N.init <- c(1, rep(0, mesh.ls$n_mesh_p-1))
            ipm.OA <- build_IPM(c(ipmPar.ls, mgmt="OA"), mesh.ls)
            ipm.TURF <- build_IPM(c(ipmPar.ls, mgmt="TURF"), mesh.ls)
            N.OA <- N.init
            N.TURF <- N.init
            
            for(j in 1:dd.i$nYrs) {
              if(j == dd.i$harvestYr) {
                N.OA <- N.OA * 0.1
                N.TURF <- N.TURF * 0.1
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
                     s_DD=dd.i$s_DD) %>%
              saveRDS(glue("{out.dir}/lambda_{i}.rds"))
            pop.ls %>% 
              do.call('rbind', .) %>%
              mutate(sim=i,
                     s_DD=dd.i$s_DD,
                     size=exp(size_ln),
                     stage=if_else(size>10, "Adult", "Juvenile"),
                     sizeClass=case_when(size < 10 ~ "0-10",
                                         between(size, 10, 20) ~ "10-20",
                                         between(size, 20, 30) ~ "20-30",
                                         between(size, 30, 40) ~ "30-40",
                                         between(size, 40, 50) ~ "40-50",
                                         size > 50 ~ "50+")) %>%
              
              saveRDS(glue("{out.dir}/pop_{i}.rds"))
            
            paste("Finished", i)
          }
  stopCluster(cl)
  
  
  
}





