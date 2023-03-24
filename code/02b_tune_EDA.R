



library(tidyverse); library(glue)
theme_set(theme_bw())

data.dir <- "data/chile/"

pop.df <- dir("out/tune_processed", "pop_pow_.*rds", full.names=T) %>%
  map_dfr(readRDS)

obs.df <- read_csv(glue("{data.dir}Kelp_transect.csv")) %>%
  filter(Depth < 15) %>%
  mutate(source="NERC_bl",
         Density_m2=Density,
         Upwelling=if_else(Upwelling=="upwelling", "Y", "N")) %>%
  select(source, Management, Upwelling, Date, MAE, Site, Transect, Station, Depth, Stage, Density_m2) %>%
  filter(!is.na(Density_m2)) %>%
  rename(stage=Stage)
obs.sum <- obs.df %>% 
  group_by(Upwelling, Management, stage) %>%
  # group_by(stage) %>%
  summarise(N_md=median(Density_m2),
            N_mn=mean(Density_m2),
            N_lo1=quantile(Density_m2, 0.25),
            N_hi1=quantile(Density_m2, 0.75),
            N_lo2=quantile(Density_m2, 0.025),
            N_hi2=quantile(Density_m2, 0.975),
            N_lo3=min(Density_m2),
            N_hi3=max(Density_m2)) %>%
  ungroup 


# No differences in NERC_bl adult densities by Upwelling * Management
# 
# library(brms)
# obs.brm <- obs.df %>%
#   filter(stage=="Adult") %>%
#   mutate(lDens=log1p(Density_m2))
# out1 <- brm(bf(lDens ~ Upwelling * Management + (1|Site/Transect)),
#             data=obs.brm, cores=4, control=list(adapt_delta=0.95),
#             prior=prior(normal(0,1), class="b"))
# out2 <- brm(bf(lDens ~ Upwelling * Management + (1|Site/Transect),
#                hu ~ Upwelling*Management + (1|Site/Transect)),
#             family=hurdle_lognormal(),
#             prior=c(prior(normal(0,1), class="b"),
#                     prior(normal(0,1), class="b", dpar="hu")),
#             data=obs.brm, cores=4, control=list(adapt_delta=0.95))
# out3 <- brm(bf(Density_m2 ~ Upwelling * Management + (1|Site/Transect),
#                hu ~ Upwelling*Management + (1|Site/Transect)),
#             family=hurdle_lognormal(),
#             prior=c(prior(normal(0,1), class="b"),
#                     prior(normal(0,1), class="b", dpar="hu")),
#             data=obs.brm, cores=4, control=list(adapt_delta=0.95))



pop.df %>%
  filter(year > 30) %>%
  filter(stage=="Adult") %>%
  group_by(year, Upwelling, Management, s_DD_pow, s_DD_i, stage, sim) %>%
  summarise(N=sum(N)) %>%
  group_by(year, Upwelling, Management, s_DD_pow, s_DD_i, stage) %>%
  summarise(N_mn=median(N),
            N_lo1=quantile(N, 0.25),
            N_hi1=quantile(N, 0.75),
            N_lo2=quantile(N, 0.05),
            N_hi2=quantile(N, 0.95)) %>%
  ggplot(aes(year, N_mn, colour=Management, fill=Management, linetype=Upwelling)) + 
  # geom_ribbon(aes(ymin=N_lo1, ymax=N_hi1), alpha=0.25, colour=NA) +
  # geom_ribbon(aes(ymin=N_lo2, ymax=N_hi2), alpha=0.25, colour=NA) +
  geom_line() + 
  geom_point(data=obs.sum %>% mutate(year=45)) +
  geom_linerange(data=obs.sum %>% mutate(year=45), aes(ymin=N_lo1, ymax=N_hi1), size=1.25) +
  geom_errorbar(data=obs.sum %>% mutate(year=45), aes(ymin=N_lo2, ymax=N_hi2), width=2) +
  facet_grid(s_DD_pow~s_DD_i)

pop.df %>%
  # filter(year > 30) %>%
  group_by(year, Upwelling, Management, s_DD_pow, s_DD_i, stage, sim) %>%
  summarise(N=sum(N)) %>%
  group_by(year, Upwelling, Management, s_DD_pow, s_DD_i, stage) %>%
  summarise(N_mn=median(N),
            N_lo1=quantile(N, 0.25),
            N_hi1=quantile(N, 0.75),
            N_lo2=quantile(N, 0.05),
            N_hi2=quantile(N, 0.95)) %>%
  ggplot(aes(year, N_mn, colour=Management, fill=Management, linetype=Upwelling)) + 
  # geom_ribbon(aes(ymin=N_lo1, ymax=N_hi1), alpha=0.25, colour=NA) +
  geom_line() + 
  geom_point(data=obs.sum %>% mutate(year=45)) + scale_y_continuous(trans="log1p") +
  geom_linerange(data=obs.sum %>% mutate(year=45), aes(ymin=N_lo1, ymax=N_hi1), size=1.25) +
  geom_errorbar(data=obs.sum %>% mutate(year=45), aes(ymin=N_lo2, ymax=N_hi2), width=2) +
  facet_grid(stage~s_DD_i, scales="free_y")

validation.sum <- pop.df %>% filter(year > 30) %>%
  mutate(simTransect=(sim-1) %/% 20) %>%
  group_by(year, Upwelling, Management, s_DD_pow, s_DD, s_DD_i, stage, simTransect) %>%
  summarise(N=sum(N)/20) %>%
  group_by(Upwelling, Management, s_DD_pow, s_DD, s_DD_i, stage) %>%
  summarise(N_md=median(N),
            N_mn=mean(N),
            N_lo1=quantile(N, 0.025),
            N_hi1=quantile(N, 0.975),
            N_lo1=min(N),
            N_hi1=max(N),
            N_lo2=quantile(N, 0.25),
            N_hi2=quantile(N, 0.75)) %>%
  ungroup %>%
  # full_join(obs.sum, suffix=c(".sim", ".obs"), by=c("Management", "stage")) %>%
  full_join(obs.sum, suffix=c(".sim", ".obs"), by=c("Upwelling", "Management", "stage")) %>%
  mutate(Upwelling=if_else(Upwelling=="Y", "upwelling", "non-upwelling"))
validation.sum %>%
  filter(stage=="Adult") %>%
  filter(between(s_DD_i, 14, 18)) %>%
  mutate(s_DD_pow=factor(s_DD_pow)) %>%
  ggplot() + 
  geom_ribbon(aes(x=s_DD_i, ymin=N_lo1.obs, ymax=N_hi1.obs), colour=NA, fill="grey", alpha=0.5) +
  geom_ribbon(aes(x=s_DD_i, ymin=N_lo2.obs, ymax=N_hi2.obs), colour=NA, fill="grey", alpha=0.5) +
  geom_hline(aes(yintercept=N_md.obs)) +
  geom_hline(aes(yintercept=N_mn.obs), linetype=2) +
  geom_line(aes(s_DD_i, N_md.sim, colour=s_DD_pow), position=position_dodge(width=0.5)) + 
  geom_point(aes(s_DD_i, N_md.sim, colour=s_DD_pow), position=position_dodge(width=0.5)) + 
  geom_point(aes(s_DD_i, N_mn.sim, colour=s_DD_pow), shape=1, position=position_dodge(width=0.5)) +
  geom_linerange(aes(x=s_DD_i, ymin=N_lo1.sim, ymax=N_hi1.sim, colour=s_DD_pow), 
                 position=position_dodge(width=0.5)) +
  geom_linerange(aes(x=s_DD_i, ymin=N_lo2.sim, ymax=N_hi2.sim, colour=s_DD_pow), 
                 size=1, position=position_dodge(width=0.5)) +
  scale_colour_viridis_d() +
  facet_grid(Management~Upwelling) +
  scale_y_continuous(trans="log1p") +
  labs(title="Adults: Avg of yr 30-40 (harvest yr 10), simulated 20m2",
       x="Survival density dependence parameter", y="Density (N/m2)")
ggsave("figs/tuning_sDD_simVsObs.png", width=6, height=6, units="in", dpi=300)

validation.sum %>%
  filter(stage=="Adult") %>%
  mutate(s_DD_pow=factor(s_DD_pow)) %>%
  ggplot() + 
  geom_ribbon(aes(x=s_DD, ymin=N_lo1.obs, ymax=N_hi1.obs), colour=NA, fill="grey", alpha=0.5) +
  geom_ribbon(aes(x=s_DD, ymin=N_lo2.obs, ymax=N_hi2.obs), colour=NA, fill="grey", alpha=0.5) +
  geom_hline(aes(yintercept=N_md.obs)) +
  geom_hline(aes(yintercept=N_mn.obs), linetype=2) +
  geom_line(aes(s_DD, N_md.sim, colour=s_DD_pow), position=position_dodge(width=0)) + 
  geom_point(aes(s_DD, N_md.sim, colour=s_DD_pow), position=position_dodge(width=0)) + 
  geom_point(aes(s_DD, N_mn.sim, colour=s_DD_pow), shape=1, position=position_dodge(width=0)) +
  geom_linerange(aes(x=s_DD, ymin=N_lo1.sim, ymax=N_hi1.sim, colour=s_DD_pow), 
                 position=position_dodge(width=0)) +
  geom_linerange(aes(x=s_DD, ymin=N_lo2.sim, ymax=N_hi2.sim, colour=s_DD_pow), 
                 size=1, position=position_dodge(width=0)) +
  scale_colour_viridis_d() +
  facet_grid(Management~Upwelling) +
  scale_y_continuous(trans="log1p") +
  labs(title="Adults: Avg of yr 30-40 (harvest yr 10), simulated 20m2",
       x="Survival density dependence parameter", y="Density (N/m2)")

validation.sum %>%
  filter(stage=="Juvenile") %>%
  mutate(s_DD_pow=factor(s_DD_pow)) %>%
  ggplot() + 
  geom_ribbon(aes(x=s_DD_i, ymin=N_lo1.obs, ymax=N_hi1.obs), colour=NA, fill="grey", alpha=0.5) +
  geom_ribbon(aes(x=s_DD_i, ymin=N_lo2.obs, ymax=N_hi2.obs), colour=NA, fill="grey", alpha=0.5) +
  geom_hline(aes(yintercept=N_md.obs)) +
  geom_hline(aes(yintercept=N_mn.obs), linetype=2) +
  geom_line(aes(s_DD_i, N_md.sim, colour=s_DD_pow), position=position_dodge(width=0.5)) + 
  geom_point(aes(s_DD_i, N_md.sim, colour=s_DD_pow), position=position_dodge(width=0.5)) + 
  geom_point(aes(s_DD_i, N_mn.sim, colour=s_DD_pow), shape=1, position=position_dodge(width=0.5)) +
  geom_linerange(aes(x=s_DD_i, ymin=N_lo1.sim, ymax=N_hi1.sim, colour=s_DD_pow), 
                 position=position_dodge(width=0.5)) +
  geom_linerange(aes(x=s_DD_i, ymin=N_lo2.sim, ymax=N_hi2.sim, colour=s_DD_pow), 
                 size=1, position=position_dodge(width=0.5)) +
  scale_colour_viridis_d() +
  facet_grid(Upwelling~Management) +
  scale_y_continuous(trans="log1p") +
  labs(title="Juveniles: Avg of yr 30-40 (harvest yr 10), simulated 20m2",
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


pop.df %>%
  filter(between(s_DD_pow, 0.8, 1)) %>%
  filter(s_DD_i==8) %>%
  filter(year > 10) %>%
  group_by(year, Upwelling, Management, stage, sim) %>%
  summarise(N=sum(N)) %>%
  group_by(year, Upwelling, Management, stage) %>%
  summarise(N_md=median(N),
            N_mn=mean(N),
            N_lo1=quantile(N, 0.25),
            N_hi1=quantile(N, 0.75),
            N_lo2=quantile(N, 0.025),
            N_hi2=quantile(N, 0.975)) %>%
            # N_lo2=quantile(N, 0.05),
            # N_hi2=quantile(N, 0.95)) %>%
  ungroup %>%
  full_join(obs.sum, suffix=c(".sim", ".obs"), by=c("Management", "Upwelling", "stage")) %>%
  ggplot(aes(year, N_md.sim, colour=Management, fill=Management)) + 
  geom_line() + 
  geom_ribbon(aes(year, ymin=N_lo1.sim, ymax=N_hi1.sim), alpha=0.2, colour=NA) +
  # geom_ribbon(aes(year, ymin=N_lo2.sim, ymax=N_hi2.sim), alpha=0.2, colour=NA) +
  # scale_y_log10() +
  facet_grid(stage~Upwelling, scales="free_y")

pop.df %>%
  filter(s_DD == -0.008) %>%
  filter(year > 10) %>%
  group_by(year, Management, stage, sim) %>%
  summarise(N=pmin(sum(N), 100)) %>%
  ggplot(aes(year, N, group=sim)) + 
  geom_line(alpha=0.1) + 
  facet_grid(stage~Management, scales="free_y")


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

