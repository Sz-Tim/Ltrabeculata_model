# KELPER2
# Vital rate regressions
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse)
library(glue); library(brms); library(ipmr)
library(doParallel); library(foreach)
source("code/00_fn.R")

reg.dir <- "out/regr/"
tune.dir <- "out/dd_tune_noSqrt/"
s_DD.seq <- seq(0, -0.00025, length.out=16)



# run simulations ---------------------------------------------------------

for(s in 1:length(s_DD.seq)) {
  
  out.dir <- glue("{tune.dir}/sdd-{str_pad(s, 2, 'left', '0')}/")
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
    L=dd.i$size_rng[1],
    U=dd.i$size_rng[2],
    n_mesh_p=30,
    n_yr=1
  )
  
  saveRDS(dd.i, glue("{out.dir}/dd_i.rds"))
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  foreach(i=1:nSim,
          .export=c("reg.dir", "out.dir", "dd.i"),
          .packages=c("tidyverse", "glue", "brms", "ipmr", "doParallel"), 
          .combine="c",
          .errorhandling="pass") %dopar% {
            
            source("code/00_fn.R")
            out.g <- readRDS(glue("{reg.dir}opt_g.rds"))
            out.s <- readRDS(glue("{reg.dir}opt_s.rds"))
            out.r <- readRDS(glue("{reg.dir}opt_r.rds"))
            out.r_z <- readRDS(glue("{reg.dir}opt_r_z.rds"))
            
            ipmPar.ls <- list(
              out_g=out.g,
              out_s=out.s,
              s_DD=dd.i$s_DD,
              out_r=out.r,
              out_r_z=out.r_z,
              draw=i,
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






# visualise ---------------------------------------------------------------

# obs.df <- read_csv(glue("{data.dir}Kelp_transect.csv")) %>%
#   mutate(source="NERC_bl",
#          Density_m2=Density) %>%
#   select(source, Management, Zone, Date, MAE, Site, Transect, Station, Depth, Stage, Density_m2) %>%
#   filter(!is.na(Density_m2)) %>%
#   rename(stage=Stage)
# obs.sum <- obs.df %>% 
#   group_by(Management, stage) %>%
#   summarise(N_md=median(Density_m2),
#             N_mn=mean(Density_m2),
#             N_lo1=quantile(Density_m2, 0.25),
#             N_hi1=quantile(Density_m2, 0.75),
#             N_lo2=quantile(Density_m2, 0.05),
#             N_hi2=quantile(Density_m2, 0.95),) %>%
#   ungroup 
# 
# lam.df <- map_dfr(dir(tune.dir, "lambda", recursive=T, full.names=T), readRDS)
# pop.df <- map_dfr(dir(tune.dir, "pop", recursive=T, full.names=T), readRDS)
# 
# 
# ggplot(lam.df, aes(year, lambda, colour=Management, group=paste(Management, sim))) + 
#   geom_line() + 
#   facet_wrap(~s_DD)
# lam.df %>% group_by(year, Management, s_DD) %>%
#   summarise(lambda_mn=exp(mean(log(lambda))),
#             lambda_q10=quantile(lambda, 0.05),
#             lambda_q90=quantile(lambda, 0.95)) %>%
#   ggplot(aes(year, lambda_mn, ymin=lambda_q10, ymax=lambda_q90, colour=Management, fill=Management)) + 
#   geom_ribbon(alpha=0.5, colour=NA) + 
#   geom_line() + 
#   facet_wrap(~s_DD)
# 
# pop.df %>%
#   filter(year > 30) %>%
#   group_by(year, Management, s_DD, stage, sim) %>%
#   summarise(N=sum(N)) %>%
#   group_by(year, Management, s_DD, stage) %>%
#   summarise(N_mn=median(N),
#             N_lo1=quantile(N, 0.25),
#             N_hi1=quantile(N, 0.75),
#             N_lo2=quantile(N, 0.05),
#             N_hi2=quantile(N, 0.95)) %>%
#   ggplot(aes(year, N_mn, colour=stage, fill=stage)) + 
#   geom_ribbon(aes(ymin=N_lo1, ymax=N_hi1), alpha=0.25, colour=NA) +
#   geom_ribbon(aes(ymin=N_lo2, ymax=N_hi2), alpha=0.25, colour=NA) +
#   geom_line() + 
#   geom_point(data=obs.sum %>% mutate(year=45)) +
#   geom_linerange(data=obs.sum %>% mutate(year=45), aes(ymin=N_lo1, ymax=N_hi1), size=1.25) +
#   geom_errorbar(data=obs.sum %>% mutate(year=45), aes(ymin=N_lo2, ymax=N_hi2), width=2) +
#   facet_grid(Management~s_DD)
# 
# validation.sum <- pop.df %>% filter(year > 30) %>%
#   mutate(simTransect=(sim-1) %/% 20) %>%
#   group_by(year, Management, s_DD, stage, simTransect) %>%
#   summarise(N=sum(N)/20) %>%
#   group_by(Management, s_DD, stage) %>%
#   summarise(N_md=median(N),
#             N_mn=mean(N),
#             N_lo1=quantile(N, 0.25),
#             N_hi1=quantile(N, 0.75),
#             N_lo2=quantile(N, 0.05),
#             N_hi2=quantile(N, 0.95)) %>%
#   ungroup %>%
#   full_join(obs.sum, suffix=c(".sim", ".obs"), by=c("Management", "stage"))
# validation.sum %>%
#   ggplot() + 
#   geom_ribbon(aes(x=s_DD, ymin=N_lo1.obs, ymax=N_hi1.obs), colour=NA, fill="grey", alpha=0.5) +
#   geom_ribbon(aes(x=s_DD, ymin=N_lo2.obs, ymax=N_hi2.obs), colour=NA, fill="grey", alpha=0.5) +
#   geom_hline(aes(yintercept=N_md.obs)) +
#   geom_hline(aes(yintercept=N_mn.obs), linetype=2) +
#   geom_point(aes(s_DD, N_md.sim)) + 
#   geom_point(aes(s_DD, N_mn.sim), shape=1) + 
#   geom_linerange(aes(x=s_DD, ymin=N_lo1.sim, ymax=N_hi1.sim), size=1) +
#   geom_linerange(aes(x=s_DD, ymin=N_lo2.sim, ymax=N_hi2.sim)) +
#   facet_grid(Management~stage) +
#   # scale_y_continuous(trans="log1p") +
#   labs(title="Avg of yr 30-40 (harvest yr 10), simulated 20m2",
#        x="Survival density dependence parameter", y="Density (N/m2)")
# ggsave("figs/tuning_sDD_simVsObs.png", width=6, height=6, units="in", dpi=300)
# 
# validation.yr <- pop.df %>% filter(year > 30) %>%
#   group_by(year, Management, s_DD, stage, sim) %>%
#   summarise(N=sum(N)) %>%
#   group_by(year, Management, s_DD, stage) %>%
#   summarise(N_md=median(N),
#             N_mn=mean(N),
#             N_lo1=quantile(N, 0.25),
#             N_hi1=quantile(N, 0.75),
#             N_lo2=quantile(N, 0.05),
#             N_hi2=quantile(N, 0.95)) %>%
#   ungroup %>%
#   full_join(obs.sum, suffix=c(".sim", ".obs"), by=c("Management", "stage"))

