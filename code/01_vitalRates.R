# KELPER2
# Vital rate regressions
# Tim Szewczyk




# setup -------------------------------------------------------------------

library(tidyverse); library(brms); library(glue); library(lubridate); library(loo)

data.dir <- "data/chile/"


Growth_Surv <- read_csv(glue("{data.dir}Growth&Survival.csv")) %>%
  mutate(Site=factor(Site),
         Management=factor(Management),
         Plant=factor(Plant),
         alive=as.numeric(Holdfast > 0),
         Date_Harvest=Date_harvest_dmy %>%
           str_replace("Dic", "01-12") %>%
           str_replace("Nov", "01-11") %>%
           str_replace("Ene", "01-01") %>%
           str_replace("Feb", "01-02") %>%
           dmy(),
         Date_Obs=Date_monitoring %>%
           str_replace("Abr", "01-04") %>%
           str_replace("Oct", "01-10") %>%
           str_replace("May", "01-05") %>%
           dmy()) %>%
  mutate(days_since_harvest=as.numeric(Date_Obs - Date_Harvest)) %>%
  filter(!is.na(Site))

survey.i <- Growth_Surv %>%
  group_by(Site, Management, Patch, Monitoring) %>%
  summarise(Date_Harvest=min(Date_Harvest),
            Date_Obs=min(Date_Obs)) %>%
  arrange(Management, Site, Patch, Monitoring)
survey.i <- survey.i %>%
  bind_rows(., survey.i %>% slice_head(n=1) %>%
              mutate(Monitoring=0,
                     Date_Obs=Date_Harvest)) %>%
  arrange(Management, Site, Patch, Monitoring) %>%
  mutate(days_since_harvest=as.numeric(Date_Obs - Date_Harvest),
         days_since_prevObs=as.numeric(Date_Obs - lag(Date_Obs)),
         days_since_prevObs=replace_na(days_since_prevObs, 0),
         days_to_nextObs=as.numeric(lead(Date_Obs) - Date_Obs))

Growth_Surv <- Growth_Surv %>%
  select(Site, Management, Patch, Monitoring, 
         Plant, Holdfast, Stipes, Total_Length) %>%
  left_join(., survey.i) %>%
  arrange(Site, Patch, Plant, Monitoring) %>%
  mutate(prevSurv=Monitoring-1) %>%
  arrange(Site, Patch, Plant, Monitoring) %>%
  group_by(Site, Patch, Plant) %>%
  mutate(doesSurvive=lead(Holdfast>0),
         nextHoldfast=lead(Holdfast),
         nextStipes=lead(Stipes),
         nextLength=lead(Total_Length),
         maxAge_days=cumsum(days_since_prevObs),
         maxAge_yrs=maxAge_days/365,
         surv_01=as.numeric(doesSurvive),
         logHoldfast=log(Holdfast)) %>%
  ungroup





# survival ----------------------------------------------------------------

surv.df <- Growth_Surv %>% filter(!is.na(doesSurvive))

form.s <- list(
  m1=bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
        pSurvYr ~ logHoldfast + (1|Site/Patch/Plant),
        nl=T),
  m2=bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
        pSurvYr ~ logHoldfast + Management + (1|Site/Patch/Plant),
        nl=T),
  m3=bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
        pSurvYr ~ logHoldfast*Management + (1|Site/Patch/Plant),
        nl=T),
  m4=bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
        pSurvYr ~ logHoldfast + (1+logHoldfast|Site/Patch/Plant),
        nl=T),
  m5=bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
        pSurvYr ~ logHoldfast + Management + (1+logHoldfast+Management|Site/Patch/Plant),
        nl=T),
  m6=bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
        pSurvYr ~ logHoldfast*Management + (1+logHoldfast*Management|Site/Patch/Plant),
        nl=T)
)

out.s <- map(form.s, 
             ~brm(.x, data=surv.df, cores=4, save_pars=save_pars(all=T), 
                  family="bernoulli", init=0, control=list(adapt_delta=0.9),
                  prior=c(prior(normal(0, 0.5), nlpar="pSurvYr"),
                          prior(normal(-1, 1), nlpar="pSurvYr", coef="Intercept"),
                          prior(normal(0, 0.5), nlpar="pSurvYr", class="sd", lb=0, 
                                group="Site"),
                          prior(normal(0, 0.1), nlpar="pSurvYr", class="sd", lb=0,
                                group="Site:Patch"),
                          prior(normal(0, 0.25), nlpar="pSurvYr", class="sd", lb=0,
                                group="Site:Patch:Plant"))))
out.s <- map(out.s, ~add_criterion(.x, "loo", moment_match=T))
loo.s <- map(out.s, ~loo(.x, cores=4)) %>% loo_compare()

opt.s <- loo.s %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -4) %>%
  arrange(model)

saveRDS(out.s[[opt.s$model[1]]], "out/regr/opt_s.rds")



out.s <- brm(bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
                pSurvYr ~ logHoldfast*Management + 
                  (1+logHoldfast*Management|Site/Patch/Plant),
                nl=T),
             data=surv.df, init=0, family="bernoulli", chains=4, cores=4,
             prior=c(prior(normal(0, 0.5), nlpar="pSurvYr"),
                     prior(normal(-1, 1), nlpar="pSurvYr", coef="Intercept"),
                     prior(normal(0, 0.5), nlpar="pSurvYr", class="sd", lb=0, 
                           group="Site"),
                     prior(normal(0, 0.1), nlpar="pSurvYr", class="sd", lb=0,
                           group="Site:Patch"),
                     prior(normal(0, 0.25), nlpar="pSurvYr", class="sd", lb=0,
                           group="Site:Patch:Plant")), 
             file="out/regr/survival.rds")




# growth ------------------------------------------------------------------

grow.df <- Growth_Surv %>% 
  filter(doesSurvive) %>%
  mutate(annHoldfastGrowth=(nextHoldfast-Holdfast)/days_to_nextObs*365,
         annStipeGrowth=(nextStipes-Stipes)/days_to_nextObs*365,
         annLengthGrowth=(nextLength-Total_Length)/days_to_nextObs*365,
         nextHoldfastAnnual=Holdfast+annHoldfastGrowth,
         logNextHoldfastAnnual=log(nextHoldfastAnnual))

form.g <- list(
  m1=bf(logNextHoldfastAnnual ~ logHoldfast + (1|Site/Patch/Plant)),
  m2=bf(logNextHoldfastAnnual ~ logHoldfast + Management + (1|Site/Patch/Plant)),
  m3=bf(logNextHoldfastAnnual ~ logHoldfast * Management + (1|Site/Patch/Plant)),
  m4=bf(logNextHoldfastAnnual ~ logHoldfast + (1 + logHoldfast|Site/Patch/Plant)),
  m5=bf(logNextHoldfastAnnual ~ logHoldfast + Management + (1 + logHoldfast + Management|Site/Patch/Plant)),
  m6=bf(logNextHoldfastAnnual ~ logHoldfast * Management + (1 + logHoldfast * Management|Site/Patch/Plant))
)

out.g <- map(form.g, 
                 ~brm(.x, data=grow.df, cores=4, save_pars=save_pars(all=T), 
                      control=list(adapt_delta=0.9),
                      prior=c(prior(normal(0, 0.1), class="sd", lb=0, group="Site"),
                              prior(normal(0, 0.1), class="sd", lb=0, group="Site:Patch"),
                              prior(normal(0, 0.1), class="sd", lb=0, group="Site:Patch:Plant"),
                              prior(normal(0, 1), class="b"),
                              prior(normal(2, 1), class="Intercept"))))
out.g <- map(out.g, ~add_criterion(.x, "loo", moment_match=T))
loo.g <- map(out.g, ~loo(.x, cores=4)) %>% loo_compare()

opt.g <- loo.g %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -4) %>%
  arrange(model)

saveRDS(out.g[[opt.g$model[1]]], "out/regr/opt_g.rds")








# recruitment rate --------------------------------------------------------

rcrDens.df <- read_csv(glue("{data.dir}Patches_quad.csv")) %>%
  filter(Species=="Lessonia trabeculata") %>%
  mutate(Site=factor(Site),
         Zone=factor(Zone),
         Management=factor(Management),
         Upwelling=factor(Upwelling),
         Stage=factor(Stage)) %>% 
  filter(Quad.area==1) %>%
  select(-...1) %>%
  pivot_wider(names_from=Stage, values_from=Density) %>%
  mutate(N=Juvenile+Adult,
         Int=1)


# Problem: This uses Adult[t] to predict Juvenile[t]
# But I want 1) new recruits, and 2) Adult[t-1] to predict [t]
form.r <- list(
  m0=bf(Juvenile ~ 0 + Int + (1|Site),
        zi ~ 0 + Int + (1|Site)),
  m1=bf(Juvenile ~ 0 + Int + Adult + (1|Site),
        zi ~ 0 + Int + (1|Site)),
  m2=bf(Juvenile ~ 0 + Int + Adult + (1|Site),
        zi ~ 0 + Int + Adult + (1|Site)),
  m3=bf(Juvenile ~ 0 + Int + Adult + Management + (1|Site),
        zi ~ 0 + Int + Adult + (1|Site)),
  m4=bf(Juvenile ~ 0 + Int + Adult + Management + (1|Site),
        zi ~ 0 + Int + Adult + Management + (1|Site))
)
out.r <- map(form.r, 
               ~brm(.x, data=rcrDens.df, cores=4, save_pars=save_pars(all=T), 
                    family=zero_inflated_poisson(), control=list(adapt_delta=0.9),
                    prior=c(prior(normal(0, 0.1), class="sd", lb=0, group="Site"),
                            prior(normal(-1, 1), class="b"),
                            prior(normal(0.5, 1), class="b", coef="Int"),
                            prior(normal(0, 0.1), class="sd", lb=0, group="Site", dpar="zi"),
                            prior(normal(0.5, 0.5), class="b", dpar="zi"),
                            prior(normal(0, 1), class="b", coef="Int", dpar="zi"))))
out.r <- map(out.r, ~add_criterion(.x, "loo", moment_match=T))
loo.r <- map(out.r, ~loo(.x, cores=4)) %>% loo_compare()

opt.r <- loo.r %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -4) %>%
  arrange(model)

saveRDS(out.r[[opt.r$model[1]]], "out/regr/opt_r.rds")

out.r <- brm(bf(Juvenile ~ Adult + (1|Site),
                zi ~ Adult + (1|Site)), 
             family=zero_inflated_poisson(), data=rcrDens.df, cores=4,
             prior=c(prior(normal(0, 0.1), class="sd", lb=0, group="Site"),
                     prior(normal(0.5, 1), class="Intercept"),
                     prior(normal(-1, 1), class="b", ub=0),
                     prior(normal(0, 0.1), class="sd", lb=0, group="Site", dpar="zi"),
                     prior(normal(0, 1), class="Intercept", dpar="zi"),
                     prior(normal(0.5, 0.5), class="b", dpar="zi")), 
             file="temp/out_r_ZIP.rds")
pp_check(out.r, ndraws=100) + scale_x_continuous(trans="log1p")


rcrSize.df <- Growth_Surv %>% 
  select(Site, Management, Patch, Monitoring, 
         Plant, Holdfast, Stipes, Total_Length) %>%
  left_join(., survey.i) %>%
  arrange(Site, Patch, Plant, Monitoring) %>%
  mutate(prevSurv=Monitoring-1) %>%
  arrange(Site, Patch, Plant, Monitoring) %>%
  filter(Monitoring==1,
         Holdfast < 10) %>% 
  mutate(logHoldfast=log(Holdfast))

form.r_z <- list(
  m0=bf(logHoldfast ~ 1),
  m1=bf(logHoldfast ~ 1 + (1|Site)),
  m2=bf(logHoldfast ~ Management),
  m3=bf(logHoldfast ~ Management + (1|Site)),
  m4=bf(logHoldfast ~ Management + (1+Management|Site))
)

out.r_z <- map(form.r_z, 
             ~brm(.x, data=rcrSize.df, cores=4, save_pars=save_pars(all=T), 
                  control=list(adapt_delta=0.9)))
out.r_z <- map(out.r_z, ~add_criterion(.x, "loo", moment_match=T))
loo.r_z <- map(out.r_z, ~loo(.x, cores=4)) %>% loo_compare()

opt.r_z <- loo.r_z %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -4) %>%
  arrange(model)

saveRDS(out.r_z[[opt.r_z$model[1]]], "out/regr/opt_r_z.rds")

out.r_z <- brm(logHoldfast ~ Management + (1|Site), 
               data=rcrSize.df, cores=4,
               prior=c(prior(normal(0, 0.5), class="sd", lb=0),
                       prior(normal(0, 1), class="b")),
               file="temp/out_r_z.rds")








# Allometry ---------------------------------------------------------------

allom.df <- list(
  Cata_ph=read_csv(glue("{data.dir}Morfo_DW_C.csv")) %>%
    select(-c(1,11,20:30)) %>%
    rename(Holdfast=Diam_Holdfast,
           Height=High_Holdfast,
           LengthTot=LT,
           WeightTot=TW_sum,
           Weight=Weight_Str,
           Length=LT_stipes,
           Dichot=Dichot_stipes,
           Diam=Diam_stipes) %>%
    pivot_wider(names_from="Structure", values_from=c(8:17)) %>%
    select(!where(~all(is.na(.x)))) %>%
    select(-starts_with("Structure_")) %>%
    mutate(source="Cata_ph") %>%
    rename(Stipes=Stipes_Stipes),
  IFOP_allom=read_csv(glue("{data.dir}IFOP_allometry.csv")) %>% 
    select(-16) %>%
    mutate(Date=ymd(paste(Year, Month, Day, sep="-")),
           Weight=1e3*Weight,
           source="IFOP_allom") %>%
    rename(LengthTot=TL,
           WeightTot=Weight)
) %>%
  do.call('bind_rows', .) %>%
  mutate(logHoldfast=log(Holdfast),
         logWeight=log(WeightTot),
         logLength=log(LengthTot),
         Management=case_when(Management!="AMERB" ~ Management,
                              Management=="AMERB" ~ "TURF")) %>%
  mutate(Site=str_to_lower(Site)) %>%
  mutate(Site2=factor(Site, levels=unique(Site)),
         Site2=factor(Site2, labels=janitor::make_clean_names(levels(Site2))))



# From holdfast diameter to biomass

# candidate models of increasing complexity
form.hf_wt <- list(
  m0=bf(logWeight ~ logHoldfast),
  m1=bf(logWeight ~ logHoldfast + (1|Site2)),
  m2=bf(logWeight ~ logHoldfast + Management + (1|Site2)),
  m3=bf(logWeight ~ logHoldfast * Management + (1|Site2)),
  m4=bf(logWeight ~ logHoldfast + (1+logHoldfast|Site2)),
  m5=bf(logWeight ~ logHoldfast + Management + (1 + logHoldfast + Management|Site2)),
  m6=bf(logWeight ~ logHoldfast * Management + (1 + logHoldfast * Management|Site2))
)

out.hf_wt <- map(form.hf_wt, 
                  ~brm(.x, data=allom.df, cores=4, save_pars=save_pars(all=T), 
                       control=list(adapt_delta=0.9)))
out.hf_wt <- map(out.hf_wt, ~add_criterion(.x, "loo", moment_match=T))
loo.hf_wt <- map(out.hf_wt, ~loo(.x, cores=4)) %>% loo_compare()

opt.hf_wt <- loo.hf_wt %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -4) %>%
  arrange(model)

saveRDS(out.hf_wt[[opt.hf_wt$model[1]]], "out/regr/opt_allom_hw.rds")




# From holdfast diameter to total length

form.hf_len <- list(
  m0=bf(logLength ~ logHoldfast),
  m1=bf(logLength ~ logHoldfast + (1|Site2)),
  m2=bf(logLength ~ logHoldfast + I(logHoldfast^2) + (1|Site2)),
  m3=bf(logLength ~ logHoldfast + I(logHoldfast^2) + Management + (1|Site2)),
  m4=bf(logLength ~ logHoldfast + (1+logHoldfast|Site2)),
)

out.hf_len <- map(form.hf_len, 
                  ~brm(.x, data=allom.df, cores=4, save_pars=save_pars(all=T), 
                       control=list(adapt_delta=0.9)))
out.hf_len <- map(out.hf_len, ~add_criterion(.x, "loo", moment_match=T))
loo.hf_len <- map(out.hf_len, loo) %>% loo_compare()

opt.hf_len <- loo.hf_len %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -4) %>%
  arrange(model)

saveRDS(out.hf_wt[[opt.hf_wt$model[1]]], "out/regr/opt_alom_hw.rds")
