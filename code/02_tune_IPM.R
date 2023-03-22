# KELPER2
# Vital rate regressions
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse)
library(glue); library(brms); library(ipmr)
library(doParallel); library(foreach)
source("code/00_fn.R")

reg.dir <- "out/regr/"
nSim <- 500
cores <- 70
N.init_rcr <- 5
rcr_max <- 100
# s_DD * Area^s_DD_pow
nPow <- 5
s_DD.pars <- list(pow=seq(0.5, 1, length.out=nPow),
                  low_lim=-exp(seq(log(0.01), log(0.0025), length.out=nPow)))
s_DD.df <- map2_dfr(s_DD.pars$pow, s_DD.pars$low_lim,
                    ~tibble(s_DD=seq(.y, 0, length.out=20),
                            s_DD_pow=.x) %>%
                      mutate(s_DD_i=row_number())) %>%
  mutate(pow.dir=glue("out/dd_tune_pow_{s_DD_pow}/"),
         out.dir=glue("{pow.dir}sdd-{str_pad(s_DD_i, 2, 'left', '0')}"))
ipm.scenarios <- expand_grid(mgmt=c("OA", "TURF"), upwell=c("N", "Y"))
post_draw <- sample.int(4000, nSim)



# run simulations ---------------------------------------------------------

for(s in 1:nrow(s_DD.df)) {
  
  dir.create(s_DD.df$out.dir[s], recursive=T, showWarnings=F)
  dd.i <- list(size_rng=log(c(0.1, 50)),
               nYrs=40,
               harvestYr=10,
               s_DD=s_DD.df$s_DD[s],
               s_DD_pow=s_DD.df$s_DD_pow[s],
               ndraws=1,
               rcr_max=rcr_max
  )
  mesh.ls=list(L=dd.i$size_rng[1],
               U=dd.i$size_rng[2],
               n_mesh_p=30,
               n_yr=1
  )
  saveRDS(dd.i, glue("{s_DD.df$out.dir[s]}/dd_i.rds"))
  
  
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach(i=1:nSim,
          .export=c("reg.dir", "ipm.scenarios", "s_DD.df", "dd.i", "s", "post_draw"),
          .packages=c("tidyverse", "glue", "brms", "ipmr", "doParallel"), 
          .combine="c",
          .errorhandling="pass") %dopar% {
            
            # setup -------------------------------------------------------

            source("code/00_fn.R")
            
            ipmPar.ls <- list(
              out_g=readRDS(glue("{reg.dir}opt_g.rds")),
              out_s=readRDS(glue("{reg.dir}opt_s.rds")),
              out_r=readRDS(glue("{reg.dir}opt_r.rds")),
              out_r_z=readRDS(glue("{reg.dir}opt_r_z.rds")),
              s_DD=dd.i$s_DD,
              s_DD_pow=dd.i$s_DD_pow,
              draw=post_draw[i],
              adultThresh=10,
              rcr_max=dd.i$rcr_max
            )
            
            pop.ls <- map(1:nrow(ipm.scenarios), ~vector("list", length=dd.i$nYrs))
            
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

              for(k in 1:dd.i$nYrs) {
                if(k == dd.i$harvestYr) {
                  N <- N * 0.1
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
                }
                
                # Ipmlement recruitment
                ipm.j$pop_state$n_z[,2] <- calc_new_recruits(ipm.j, pars.j, size_ln)

                pop.ls[[j]][[k]] <- tibble(year=k,
                                           size_ln=mesh_info$z_1[1:mesh.ls$n_mesh_p],
                                           N=c(ipm.j$pop_state$n_z[,2]),
                                           Management=ipm.scenarios$mgmt[j],
                                           Upwelling=ipm.scenarios$upwell[j])
                N <- c(ipm.j$pop_state$n_z[,2])
              }
            }
            
            map_dfr(pop.ls, 
                    ~do.call('rbind', .x) %>%
                      mutate(sim=i,
                             s_DD_i=s_DD.df$s_DD_i[s],
                             s_DD=dd.i$s_DD, 
                             s_DD_pow=dd.i$s_DD_pow,
                             size=exp(size_ln),
                             stage=if_else(size>10, "Adult", "Juvenile"),
                             sizeClass=case_when(size < 10 ~ "0-10",
                                                 between(size, 10, 20) ~ "10-20",
                                                 between(size, 20, 30) ~ "20-30",
                                                 between(size, 30, 40) ~ "30-40",
                                                 between(size, 40, 50) ~ "40-50",
                                                 size > 50 ~ "50+"))) %>%
              saveRDS(glue("{s_DD.df$out.dir[s]}/pop_{i}.rds"))
            
          }
  stopCluster(cl)
}






# summarise ---------------------------------------------------------------

dir.create("out/tune_processed")
s_DD.pows <- s_DD.df %>% group_by(s_DD_pow) %>% slice_head(n=1)
walk(1:nrow(s_DD.pows), 
    ~dir(s_DD.pows$pow.dir[.x], "pop", recursive=T, full.names=T) %>%
      map_dfr(readRDS) %>%
      saveRDS(glue("out/tune_processed/pop_pow_{s_DD.pows$s_DD_pow[.x]}.rds")))



