# KELPER2
# Vital rate regressions
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse)
library(glue); library(brms); library(ipmr)
source("code/00_fn.R")

reg.dir <- "out/regr/"
h_i <- read_csv("data/harvest_i.csv")




# load output -------------------------------------------------------------

out.f <- map(dir("out", "_H"), ~dir(glue("out/{.x}/"), "harvest_.*rds", 
                                    recursive=T, full.names=T)) %>%
  unlist
H.df <- map_dfr(out.f, ~readRDS(.x) %>% mutate(fname=.x)) %>% 
  mutate(scenario=str_sub(fname, 8, 9)) %>%
  left_join(., h_i, by="scenario")

out.f <- map(dir("out", "_H"), ~dir(glue("out/{.x}/"), "pop_.*rds", 
                                    recursive=T, full.names=T)) %>%
  unlist
pop.df <- map_dfr(out.f, ~readRDS(.x) %>% mutate(fname=.x)) %>% 
  mutate(scenario=str_sub(fname, 8, 9)) %>%
  left_join(., h_i, by="scenario")



H.df %>%
  mutate(H_cumul=H_cumul/1e3) %>%
  # mutate(H_cumul=log1p(H_cumul)) %>%
  group_by(scenario, Upwelling, Management, year, description) %>%
  summarise(H_cumul_mn=median(H_cumul),
            H_cumul_lo=quantile(H_cumul, probs=0.25),
            H_cumul_hi=quantile(H_cumul, probs=0.75)) %>%
  ggplot(aes(year, H_cumul_mn, ymin=H_cumul_lo, ymax=H_cumul_hi, colour=Management, fill=Management)) + 
  geom_line() + 
  geom_ribbon(colour=NA, alpha=0.5) +
  facet_grid(Upwelling~description)

H.df %>%
  mutate(H_cumul=H_cumul/1e3) %>%
  # mutate(H_cumul=log1p(H_cumul)) %>%
  group_by(scenario, Upwelling, Management, year, description) %>%
  summarise(H_cumul_mn=median(H_cumul),
            H_cumul_lo=quantile(H_cumul, probs=0.25),
            H_cumul_hi=quantile(H_cumul, probs=0.75)) %>%
  ggplot(aes(year, H_cumul_mn, ymin=H_cumul_lo, ymax=H_cumul_hi, colour=description, fill=description)) + 
  geom_line() + 
  geom_ribbon(colour=NA, alpha=0.5) +
  scale_colour_brewer(type="qual", palette="Paired") +
  scale_fill_brewer(type="qual", palette="Paired") +
  facet_grid(Upwelling~Management)

pop.df %>%
  group_by(scenario, description, sim, Upwelling, Management, stage, year) %>%
  summarise(N=log(sum(N))) %>%
  group_by(scenario, description, Upwelling, Management, stage, year) %>%
  summarise(N=exp(mean(N))) %>%
  ungroup() %>%
  # mutate(N=if_else(N < 0.1, 0.01, N)) %>%
  ggplot(aes(year, N, fill=stage)) +
  geom_bar(stat="identity", position="fill") + 
  # scale_y_log10() +
  facet_grid(Management*Upwelling~description, scales="free_y")

pop.df %>%
  group_by(scenario, description, sim, Upwelling, Management, stage, year) %>%
  summarise(N=log(sum(N))) %>%
  group_by(scenario, description, Upwelling, Management, stage, year) %>%
  summarise(N=exp(mean(N))) %>%
  ungroup() %>%
  # mutate(N=if_else(N < 0.1, 0.01, N)) %>%
  ggplot(aes(year, N, colour=stage, linetype=Upwelling)) +
  geom_line() + 
  # scale_y_log10() +
  facet_grid(Management~description, scales="free_y")

pop.df %>%
  group_by(scenario, description, sim, Upwelling, Management, sizeClass, year) %>%
  summarise(N=log(sum(N))) %>%
  group_by(scenario, description, Upwelling, Management, sizeClass, year) %>%
  summarise(N=exp(mean(N))) %>%
  ungroup() %>%
  # mutate(N=if_else(N < 0.1, 0.01, N)) %>%
  ggplot(aes(year, N, colour=sizeClass, linetype=Upwelling)) +
  geom_line() + 
  scale_fill_viridis_d() +
  # scale_y_log10() +
  facet_grid(Management~description, scales="free_y")


# 
# lam.df <- map_dfr(dir(h_i$outDir[h], "lambda", recursive=T, full.names=T), readRDS)
# pop.df <- map_dfr(dir(h_i$outDir[h], "pop", recursive=T, full.names=T), readRDS)
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
#   filter(year >= 10) %>%
#   group_by(year, Management, stage, sim) %>%
#   summarise(N=sum(N)) %>%
#   # mutate(N=log(N)) %>%
#   group_by(year, Management, stage) %>%
#   summarise(N_mn=median(N),
#             N_lo1=HDInterval::hdi(N, 0.8)[1],
#             N_hi1=HDInterval::hdi(N, 0.8)[2],
#             N_lo2=HDInterval::hdi(N, 0.95)[1],
#             N_hi2=HDInterval::hdi(N, 0.95)[2]) %>%
#   ggplot(aes(year, N_mn, colour=Management, fill=Management)) +
#   geom_ribbon(aes(ymin=N_lo1, ymax=N_hi1), alpha=0.25, colour=NA) +
#   geom_ribbon(aes(ymin=N_lo2, ymax=N_hi2), alpha=0.25, colour=NA) +
#   geom_line() +
#   scale_y_continuous(trans="log1p") +
#   facet_grid(stage~., scales="free_y")
# 
# pop.df %>%
#   filter(year >= 10) %>%
#   group_by(year, Management, stage, sim) %>%
#   summarise(weight=sum(N * exp(weightClass_ln))) %>%
#   group_by(year, Management, stage) %>%
#   summarise(weight_mn=median(weight),
#             weight_lo1=HDInterval::hdi(weight, 0.8)[1],
#             weight_hi1=HDInterval::hdi(weight, 0.8)[2],
#             weight_lo2=HDInterval::hdi(weight, 0.90)[1],
#             weight_hi2=HDInterval::hdi(weight, 0.90)[2]) %>%
#   ggplot(aes(year, weight_mn, colour=stage, fill=stage)) +
#   geom_ribbon(aes(ymin=weight_lo1, ymax=weight_hi1), alpha=0.25, colour=NA) +
#   geom_ribbon(aes(ymin=weight_lo2, ymax=weight_hi2), alpha=0.25, colour=NA) +
#   geom_line() +
#   # scale_y_continuous(trans="log") +
#   facet_grid(Management~.)
# 
