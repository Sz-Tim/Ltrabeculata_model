



library(tidyverse); library(glue)
theme_set(theme_bw())

data.dir <- "data/chile/"
tune.dir <- "out/dd_tune/"


obs.df <- read_csv(glue("{data.dir}Kelp_transect.csv")) %>%
  mutate(source="NERC_bl",
         Density_m2=Density) %>%
  select(source, Management, Zone, Date, MAE, Site, Transect, Station, Depth, Stage, Density_m2) %>%
  filter(!is.na(Density_m2)) %>%
  rename(stage=Stage)
obs.sum <- obs.df %>% 
  group_by(Management, stage) %>%
  summarise(N_md=median(Density_m2),
            N_mn=mean(Density_m2),
            N_lo1=quantile(Density_m2, 0.25),
            N_hi1=quantile(Density_m2, 0.75),
            N_lo2=quantile(Density_m2, 0.05),
            N_hi2=quantile(Density_m2, 0.95),) %>%
  ungroup 

lam.df <- map_dfr(dir(tune.dir, "lambda", recursive=T, full.names=T), readRDS)
pop.df <- map_dfr(dir(tune.dir, "pop", recursive=T, full.names=T), readRDS)


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
  filter(year > 30) %>%
  group_by(year, Management, s_DD, stage, sim) %>%
  summarise(N=sum(N)) %>%
  group_by(year, Management, s_DD, stage) %>%
  summarise(N_mn=median(N),
            N_lo1=quantile(N, 0.25),
            N_hi1=quantile(N, 0.75),
            N_lo2=quantile(N, 0.05),
            N_hi2=quantile(N, 0.95)) %>%
  ggplot(aes(year, N_mn, colour=stage, fill=stage)) + 
  geom_ribbon(aes(ymin=N_lo1, ymax=N_hi1), alpha=0.25, colour=NA) +
  geom_ribbon(aes(ymin=N_lo2, ymax=N_hi2), alpha=0.25, colour=NA) +
  geom_line() + 
  geom_point(data=obs.sum %>% mutate(year=45)) +
  geom_linerange(data=obs.sum %>% mutate(year=45), aes(ymin=N_lo1, ymax=N_hi1), size=1.25) +
  geom_errorbar(data=obs.sum %>% mutate(year=45), aes(ymin=N_lo2, ymax=N_hi2), width=2) +
  facet_grid(Management~s_DD)

validation.sum <- pop.df %>% filter(year > 30) %>%
  mutate(simTransect=(sim-1) %/% 20) %>%
  group_by(year, Management, s_DD, stage, simTransect) %>%
  summarise(N=sum(N)/20) %>%
  group_by(Management, s_DD, stage) %>%
  summarise(N_md=median(N),
            N_mn=mean(N),
            N_lo1=quantile(N, 0.25),
            N_hi1=quantile(N, 0.75),
            N_lo2=quantile(N, 0.05),
            N_hi2=quantile(N, 0.95)) %>%
  ungroup %>%
  full_join(obs.sum, suffix=c(".sim", ".obs"), by=c("Management", "stage"))
validation.sum %>%
  ggplot() + 
  geom_ribbon(aes(x=s_DD, ymin=N_lo1.obs, ymax=N_hi1.obs), colour=NA, fill="grey", alpha=0.5) +
  geom_ribbon(aes(x=s_DD, ymin=N_lo2.obs, ymax=N_hi2.obs), colour=NA, fill="grey", alpha=0.5) +
  geom_hline(aes(yintercept=N_md.obs)) +
  geom_hline(aes(yintercept=N_mn.obs), linetype=2) +
  geom_point(aes(s_DD, N_md.sim)) + 
  geom_point(aes(s_DD, N_mn.sim), shape=1) + 
  geom_linerange(aes(x=s_DD, ymin=N_lo1.sim, ymax=N_hi1.sim), size=1) +
  geom_linerange(aes(x=s_DD, ymin=N_lo2.sim, ymax=N_hi2.sim)) +
  facet_grid(Management~stage) +
  # scale_y_continuous(trans="log1p") +
  labs(title="Avg of yr 30-40 (harvest yr 10), simulated 20m2",
       x="Survival density dependence parameter", y="Density (N/m2)")
ggsave("figs/tuning_sDD_simVsObs.png", width=6, height=6, units="in", dpi=300)

validation.yr <- pop.df %>% filter(year > 30) %>%
  group_by(year, Management, s_DD, stage, sim) %>%
  summarise(N=sum(N)) %>%
  group_by(year, Management, s_DD, stage) %>%
  summarise(N_md=median(N),
            N_mn=mean(N),
            N_lo1=quantile(N, 0.25),
            N_hi1=quantile(N, 0.75),
            N_lo2=quantile(N, 0.05),
            N_hi2=quantile(N, 0.95)) %>%
  ungroup %>%
  full_join(obs.sum, suffix=c(".sim", ".obs"), by=c("Management", "stage"))

library(gganimate) 
anim <- ggplot(validation.yr) +
  geom_ribbon(aes(x=s_DD, ymin=N_lo1.obs, ymax=N_hi1.obs), colour=NA, fill="grey", alpha=0.5) +
  geom_ribbon(aes(x=s_DD, ymin=N_lo2.obs, ymax=N_hi2.obs), colour=NA, fill="grey", alpha=0.5) +
  geom_hline(aes(yintercept=N_md.obs)) +
  geom_hline(aes(yintercept=N_mn.obs), linetype=2) +
  geom_point(aes(s_DD, N_md.sim)) + 
  geom_point(aes(s_DD, N_mn.sim), shape=1) + 
  geom_linerange(aes(x=s_DD, ymin=N_lo1.sim, ymax=N_hi1.sim), size=1) +
  geom_linerange(aes(x=s_DD, ymin=N_lo2.sim, ymax=N_hi2.sim)) +
  facet_grid(Management~stage) +
  scale_y_continuous(trans="log1p") +
  transition_states(year, 1, 0) +
  labs(title="{closest_state}", 
       x="Survival density dependence parameter", y="Density (N/m2)") 
anim_save("figs/tuning_sDD_simVsObs.gif", anim, width=6, height=6, units="in", res=300)










out.dir <- "temp/sdd-0025/"


lam.df <- map_dfr(dir(out.dir, "lambda_", full.names=T), 
                  ~readRDS(.x) %>%
                    mutate(sim=str_sub(str_remove(.x, "temp/lambda_"), 1, -5)))
lam.df %>%
  mutate(lambda=log(lambda)) %>%
  group_by(Management, year) %>%
  summarise(lam_mn=mean(lambda), 
            lam_lo=quantile(lambda, 0.025),
            lam_hi=quantile(lambda, 0.975)) %>%
  mutate(across(starts_with("lam"), exp)) %>%
  ggplot(aes(year, lam_mn, colour=Management, fill=Management)) + 
  geom_ribbon(aes(ymin=lam_lo, ymax=lam_hi), alpha=0.25, colour=NA) +
  geom_line(linewidth=1) +
  labs(x="Year", y="lambda (mean + 95% CI)")


pop.df <- map_dfr(dir(out.dir, "pop_", full.names=T), readRDS)
# ggplot(pop.df, aes(year, N, colour=size, group=paste0(sim,size))) + 
#   geom_line(alpha=0.5) + scale_colour_viridis_c() + facet_wrap(~Management)
pop.df %>%
  group_by(sim, Management, year, sizeClass) %>%
  summarise(N=sum(N)) %>%
  group_by(Management, year, sizeClass) %>%
  mutate(N=log(N)) %>%
  summarise(N_mn=mean(N), N_lo=quantile(N, 0.025), N_hi=quantile(N, 0.975)) %>%
  mutate(across(starts_with("N"), exp)) %>%
  ggplot(aes(year, N_mn, colour=sizeClass, fill=sizeClass, group=sizeClass)) + 
  geom_ribbon(aes(ymin=N_lo, ymax=N_hi), alpha=0.25, colour=NA) +
  geom_line(linewidth=1) + 
  scale_colour_viridis_d() + scale_fill_viridis_d() +
  facet_wrap(~Management) +
  labs(x="Year", y="N/m2 (mean + 95% CIs)")
pop.df %>%
  group_by(sim, Management, year, stage) %>%
  summarise(N=sum(N)) %>%
  group_by(Management, year, stage) %>%
  mutate(N=log(N)) %>%
  summarise(N_mn=mean(N), N_lo=quantile(N, 0.025), N_hi=quantile(N, 0.975)) %>%
  mutate(across(starts_with("N"), exp)) %>%
  ggplot(aes(year, N_mn, colour=stage, fill=stage, group=stage)) + 
  geom_ribbon(aes(ymin=N_lo, ymax=N_hi), alpha=0.25, colour=NA) +
  geom_line(linewidth=1) + 
  facet_wrap(~Management) +
  labs(x="Year", y="N/m2 (mean + 95% CIs)")
pop.df %>%
  group_by(sim, Management, year, stage) %>%
  summarise(N=sum(N)) %>%
  group_by(Management, year, stage) %>%
  mutate(N=log(N)) %>%
  summarise(N_mn=mean(N), N_lo=quantile(N, 0.025), N_hi=quantile(N, 0.975)) %>%
  mutate(across(starts_with("N"), exp)) %>%
  ggplot(aes(year, N_mn, colour=Management, fill=Management)) + 
  geom_ribbon(aes(ymin=N_lo, ymax=N_hi), alpha=0.25, colour=NA) +
  geom_line(linewidth=1) + 
  facet_wrap(~stage, scales="free_y") +
  labs(x="Year", y="N/m2 (mean + 95% CIs)")
pop.df %>%
  group_by(sim, Management, year, stage) %>%
  summarise(N=sum(N)) %>%
  ggplot(aes(year, N, colour=stage, group=paste(sim, stage))) + 
  geom_line(alpha=0.5) + 
  facet_wrap(~Management)
pop.df %>%
  group_by(sim, Management, year, stage) %>%
  summarise(N=sum(N)) %>%
  ggplot(aes(year, N, colour=stage, group=paste(sim, stage))) + 
  geom_line(alpha=0.5) + 
  facet_wrap(~Management) +
  scale_y_continuous(trans="log1p") 

pop.df %>%
  group_by(year, N, size, Management) %>%
  summarise(N=mean(N)) %>%
  arrange(Management, year, desc(size)) %>%
  ggplot(aes(year, N, fill=size)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_viridis_c() +
  facet_wrap(~Management)
pop.df %>%
  group_by(year, N, size_ln, Management) %>%
  summarise(N=mean(N)) %>%
  arrange(Management, year, desc(size_ln)) %>%
  ggplot(aes(year, N, fill=size_ln)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_viridis_c() +
  facet_wrap(~Management)
pop.df %>%
  group_by(year, N, sizeClass, Management) %>%
  summarise(N=mean(N)) %>%
  arrange(Management, year, desc(sizeClass)) %>%
  ggplot(aes(year, N, fill=sizeClass)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_viridis_d() +
  facet_wrap(~Management)
pop.df %>%
  group_by(year, N, stage, Management) %>%
  summarise(N=mean(N)) %>%
  arrange(Management, year) %>%
  ggplot(aes(year, N, fill=stage)) +
  geom_bar(stat="identity", position="fill") +
  facet_wrap(~Management)

pop.df %>%
  mutate(area=N * pi * (size/2)^2) %>%
  group_by(sim, Management, year) %>%
  summarise(area=sum(area)) %>%
  ggplot(aes(year, area, group=sim)) + geom_line() + facet_wrap(~Management)

