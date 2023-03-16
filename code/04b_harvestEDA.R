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




# load output -------------------------------------------------------------

out.f <- map(dir("out", "_H"), ~dir(glue("out/{.x}/"), "H_.*rds", 
                                    recursive=T, full.names=T)) %>%
  unlist
H.df <- map_dfr(out.f, ~readRDS(.x) %>% mutate(fname=.x)) %>% 
  mutate(scenario=str_sub(fname, 8, 9)) %>%
  left_join(., h_i, by="scenario") %>%
  group_by(sim, scenario, Management) %>%
  mutate(Htot=cumsum(H_grams))

out.f <- map(dir("out", "_H"), ~dir(glue("out/{.x}/"), "pop_.*rds", 
                                    recursive=T, full.names=T)) %>%
  unlist
pop.df <- map_dfr(out.f, ~readRDS(.x) %>% mutate(fname=.x)) %>% 
  mutate(scenario=str_sub(fname, 8, 9)) %>%
  left_join(., h_i, by="scenario")



H.df %>%
  mutate(Htot=Htot/1e3) %>%
  # mutate(Htot=log1p(Htot)) %>%
  group_by(scenario, Management, year, description) %>%
  summarise(Htot_mn=median(Htot),
            Htot_lo=quantile(Htot, probs=0.25),
            Htot_hi=quantile(Htot, probs=0.75)) %>%
  ggplot(aes(year, Htot_mn, ymin=Htot_lo, ymax=Htot_hi, colour=Management, fill=Management)) + 
  geom_line() + 
  geom_ribbon(colour=NA, alpha=0.5) +
  facet_wrap(~description)

H.df %>%
  mutate(Htot=Htot/1e3) %>%
  # mutate(Htot=log1p(Htot)) %>%
  group_by(scenario, Management, year, description) %>%
  summarise(Htot_mn=median(Htot),
            Htot_lo=quantile(Htot, probs=0.25),
            Htot_hi=quantile(Htot, probs=0.75)) %>%
  ggplot(aes(year, Htot_mn, ymin=Htot_lo, ymax=Htot_hi, colour=description, fill=description)) + 
  geom_line() + 
  geom_ribbon(colour=NA, alpha=0.5) +
  facet_wrap(~Management)

pop.df %>%
  group_by(scenario, description, sim, Management, stage, year) %>%
  summarise(N=sum(N), 
            wt=sum(N*exp(weightClass_ln)/1e3))


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
