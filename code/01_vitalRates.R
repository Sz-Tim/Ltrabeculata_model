# KELPER2
# Vital rate regressions
# Tim Szewczyk




# setup -------------------------------------------------------------------

library(tidyverse); library(brms); library(glue); library(lubridate); library(loo)

cores <- 10
options(mc.cores=cores)
data.dir <- "data/chile/"


Growth_Surv <- read_csv(glue("{data.dir}Growth&Survival.csv")) %>%
  left_join(read_csv("data/Cata_site_i.csv") %>% select(-Site_og)) %>%
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
  group_by(Site, Upwelling, Management, Patch, Monitoring) %>%
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
  select(Site, Upwelling, Management, Patch, Monitoring, 
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
        nl=T),
  m7=bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
        pSurvYr ~ logHoldfast + Upwelling + (1+logHoldfast+Upwelling|Site/Patch/Plant),
        nl=T),
  m8=bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
        pSurvYr ~ logHoldfast*Upwelling + (1+logHoldfast*Upwelling|Site/Patch/Plant),
        nl=T),
  m9=bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
        pSurvYr ~ logHoldfast + Management + Upwelling + (1+logHoldfast+Management+Upwelling|Site/Patch/Plant),
        nl=T),
  m10=bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
        pSurvYr ~ logHoldfast*Management*Upwelling + (1+logHoldfast*Management*Upwelling|Site/Patch/Plant),
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
out.s <- map(out.s, ~add_criterion(.x, "loo", moment_match=T, cores=cores, resp="pSurvYr"))
loo.s <- map(out.s, ~loo(.x, cores=cores, resp="pSurvYr")) %>% loo_compare()
saveRDS(loo.s, "out/regr/loo_s.rds")

opt.s <- loo.s %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -2) %>%
  mutate(model2=paste0("m", str_pad(str_sub(model, 2, -1), 2, "left", "0"))) %>%
  arrange(model2)
opt2.s <- loo.s %>% as_tibble(rownames="model") %>% 
  arrange(desc(as.numeric(elpd_diff)))

# saveRDS(out.s[[opt.s$model[1]]], "out/regr/opt_s.rds")
# saveRDS(out.s$m6, "out/regr/opt_s.rds")
saveRDS(out.s$m10, "out/regr/opt_s.rds")
saveRDS(out.s[[opt2.s$model[1]]], "out/regr/opt2_s.rds")
saveRDS(out.s, "out/regr/all_s.rds")
saveRDS(out.s$m10, "out/regr/full_s.rds")



# out.s <- brm(bf(surv_01 ~ inv_logit(pSurvYr)^maxAge_yrs,
#                 pSurvYr ~ logHoldfast*Management + 
#                   (1+logHoldfast*Management|Site/Patch/Plant),
#                 nl=T),
#              data=surv.df, init=0, family="bernoulli", chains=4, cores=4,
#              prior=c(prior(normal(0, 0.5), nlpar="pSurvYr"),
#                      prior(normal(-1, 1), nlpar="pSurvYr", coef="Intercept"),
#                      prior(normal(0, 0.5), nlpar="pSurvYr", class="sd", lb=0, 
#                            group="Site"),
#                      prior(normal(0, 0.1), nlpar="pSurvYr", class="sd", lb=0,
#                            group="Site:Patch"),
#                      prior(normal(0, 0.25), nlpar="pSurvYr", class="sd", lb=0,
#                            group="Site:Patch:Plant")), 
#              file="out/regr/survival.rds")




# growth ------------------------------------------------------------------

grow.df <- Growth_Surv %>% 
  filter(doesSurvive) %>%
  mutate(annHoldfastGrowth=(nextHoldfast-Holdfast)/days_to_nextObs*365,
         annStipeGrowth=(nextStipes-Stipes)/days_to_nextObs*365,
         annLengthGrowth=(nextLength-Total_Length)/days_to_nextObs*365,
         nextHoldfastAnnual=Holdfast+annHoldfastGrowth,
         logNextHoldfastAnnual=log(nextHoldfastAnnual)) %>%
  filter(logHoldfast > -2) # outlier in OA non-upwelling with huge weight

form.g <- list(
  m1=bf(logNextHoldfastAnnual ~ logHoldfast + (1|Site/Patch/Plant)),
  m2=bf(logNextHoldfastAnnual ~ logHoldfast + Management + (1|Site/Patch/Plant)),
  m3=bf(logNextHoldfastAnnual ~ logHoldfast * Management + (1|Site/Patch/Plant)),
  m4=bf(logNextHoldfastAnnual ~ logHoldfast + (1 + logHoldfast|Site/Patch/Plant)),
  m5=bf(logNextHoldfastAnnual ~ logHoldfast + Management + (1 + logHoldfast + Management|Site/Patch/Plant)),
  m6=bf(logNextHoldfastAnnual ~ logHoldfast * Management + (1 + logHoldfast * Management|Site/Patch/Plant)),
  m7=bf(logNextHoldfastAnnual ~ logHoldfast + Upwelling + (1 + logHoldfast + Upwelling|Site/Patch/Plant)),
  m8=bf(logNextHoldfastAnnual ~ logHoldfast * Upwelling + (1 + logHoldfast * Upwelling|Site/Patch/Plant)),
  m9=bf(logNextHoldfastAnnual ~ logHoldfast + Management + Upwelling + (1 + logHoldfast + Management + Upwelling|Site/Patch/Plant)),
  m10=bf(logNextHoldfastAnnual ~ logHoldfast * Management * Upwelling + (1 + logHoldfast * Management * Upwelling|Site/Patch/Plant))
)

out.g <- map(form.g, 
                 ~brm(.x, data=grow.df, cores=4, save_pars=save_pars(all=T), 
                      control=list(adapt_delta=0.9),
                      prior=c(prior(normal(0, 0.1), class="sd", lb=0, group="Site"),
                              prior(normal(0, 0.1), class="sd", lb=0, group="Site:Patch"),
                              prior(normal(0, 0.1), class="sd", lb=0, group="Site:Patch:Plant"),
                              prior(normal(0, 1), class="b"),
                              prior(normal(2, 1), class="Intercept"))))
out.g <- map(out.g, ~add_criterion(.x, "loo", moment_match=T, cores=cores))
loo.g <- map(out.g, ~loo(.x, cores=cores)) %>% loo_compare()
saveRDS(loo.g, "out/regr/loo_g.rds")

opt.g <- loo.g %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -2) %>%
  mutate(model2=paste0("m", str_pad(str_sub(model, 2, -1), 2, "left", "0"))) %>%
  arrange(model2)
opt2.g <- loo.g %>% as_tibble(rownames="model") %>% 
  arrange(desc(as.numeric(elpd_diff)))

# saveRDS(out.g[[opt.g$model[1]]], "out/regr/opt_g.rds")
saveRDS(out.g$m8, "out/regr/opt_g.rds")
saveRDS(out.g[[opt2.g$model[1]]], "out/regr/opt2_g.rds")
saveRDS(out.g, "out/regr/all_g.rds")
saveRDS(out.g$m10, "out/regr/full_g.rds")








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


# *Abund refers to the NegBin component and *Zi refers to the zero-inflation
form.r <- list(
  m1=bf(Juvenile ~ IntAbund + AdultAbund,
        IntAbund ~ 1 + (1|Site),
        lf(AdultAbund ~ 0 + Adult, cmc=F),
        nlf(zi ~ IntZi),
        IntZi ~ 1 + (1|Site),
        nl=T),
  m2=bf(Juvenile ~ IntAbund + AdultAbund,
        IntAbund ~ 1 + (1|Site),
        lf(AdultAbund ~ 0 + Adult, cmc=F),
        nlf(zi ~ IntZi + AdultZi),
        IntZi ~ 1 + (1|Site),
        lf(AdultZi ~ 0 + Adult, cmc=F),
        nl=T),
  m3=bf(Juvenile ~ IntAbund + AdultAbund + MgmtAbund,
        IntAbund ~ 1 + (1|Site),
        lf(AdultAbund ~ 0 + Adult, cmc=F),
        lf(MgmtAbund ~ 0 + Management, cmc=F),
        nlf(zi ~ IntZi + AdultZi),
        IntZi ~ 1 + (1|Site),
        lf(AdultZi ~ 0 + Adult, cmc=F),
        nl=T),
  m4=bf(Juvenile ~ IntAbund + AdultAbund + MgmtAbund,
        IntAbund ~ 1 + (1|Site),
        lf(AdultAbund ~ 0 + Adult, cmc=F),
        lf(MgmtAbund ~ 0 + Management, cmc=F),
        nlf(zi ~ IntZi + AdultZi + MgmtZi),
        IntZi ~ 1 + (1|Site),
        lf(AdultZi ~ 0 + Adult, cmc=F),
        lf(MgmtZi ~ 0 + Management, cmc=F),
        nl=T),
  m5=bf(Juvenile ~ IntAbund + AdultAbund + UpwellAbund,
        IntAbund ~ 1 + (1|Site),
        lf(AdultAbund ~ 0 + Adult, cmc=F),
        lf(UpwellAbund ~ 0 + Upwelling, cmc=F),
        nlf(zi ~ IntZi + AdultZi),
        IntZi ~ 1 + (1|Site),
        lf(AdultZi ~ 0 + Adult, cmc=F),
        nl=T),
  m6=bf(Juvenile ~ IntAbund + AdultAbund + UpwellAbund,
        IntAbund ~ 1 + (1|Site),
        lf(AdultAbund ~ 0 + Adult, cmc=F),
        lf(UpwellAbund ~ 0 + Upwelling, cmc=F),
        nlf(zi ~ IntZi + AdultZi + UpwellZi),
        IntZi ~ 1 + (1|Site),
        lf(AdultZi ~ 0 + Adult, cmc=F),
        lf(UpwellZi ~ 0 + Upwelling, cmc=F),
        nl=T),
  m7=bf(Juvenile ~ IntAbund + AdultAbund + UpwellAbund + MgmtAbund,
        IntAbund ~ 1 + (1|Site),
        lf(AdultAbund ~ 0 + Adult, cmc=F),
        lf(UpwellAbund ~ 0 + Upwelling, cmc=F),
        lf(MgmtAbund ~ 0 + Management, cmc=F),
        nlf(zi ~ IntZi + AdultZi),
        IntZi ~ 1 + (1|Site),
        lf(AdultZi ~ 0 + Adult, cmc=F),
        nl=T),
  m8=bf(Juvenile ~ IntAbund + AdultAbund + UpwellAbund + MgmtAbund,
        IntAbund ~ 1 + (1|Site),
        lf(AdultAbund ~ 0 + Adult, cmc=F),
        lf(UpwellAbund ~ 0 + Upwelling, cmc=F),
        lf(MgmtAbund ~ 0 + Management, cmc=F),
        nlf(zi ~ IntZi + AdultZi + UpwellZi + MgmtZi),
        IntZi ~ 1 + (1|Site),
        lf(AdultZi ~ 0 + Adult, cmc=F),
        lf(UpwellZi ~ 0 + Upwelling, cmc=F),
        lf(MgmtZi ~ 0 + Management, cmc=F),
        nl=T),
  m9=bf(Juvenile ~ IntAbund + AdultAbund + UpwellAbund + MgmtAbund + FullInteractionAbund,
        IntAbund ~ 1 + (1|Site),
        lf(AdultAbund ~ 0 + Adult, cmc=F),
        lf(UpwellAbund ~ 0 + Upwelling, cmc=F),
        lf(MgmtAbund ~ 0 + Management, cmc=F),
        lf(FullInteractionAbund ~ 0 + Upwelling:Management, cmc=F),
        nlf(zi ~ IntZi + AdultZi),
        IntZi ~ 1 + (1|Site),
        lf(AdultZi ~ 0 + Adult, cmc=F),
        nl=T),
  m10=bf(Juvenile ~ IntAbund + AdultAbund + UpwellAbund + MgmtAbund + FullInteractionAbund,
        IntAbund ~ 1 + (1|Site),
        lf(AdultAbund ~ 0 + Adult, cmc=F),
        lf(UpwellAbund ~ 0 + Upwelling, cmc=F),
        lf(MgmtAbund ~ 0 + Management, cmc=F),
        lf(FullInteractionAbund ~ 0 + Upwelling:Management, cmc=F),
        nlf(zi ~ IntZi + AdultZi + UpwellZi + MgmtZi + FullInteractionZi),
        IntZi ~ 1 + (1|Site),
        lf(AdultZi ~ 0 + Adult, cmc=F),
        lf(UpwellZi ~ 0 + Upwelling, cmc=F),
        lf(MgmtZi ~ 0 + Management, cmc=F),
        lf(FullInteractionZi ~ 0 + Upwelling:Management, cmc=F),
        nl=T)
  # m9=bf(Juvenile ~ IntAbund + FullAbund,
  #       IntAbund ~ 1 + (1|Site),
  #       lf(FullAbund ~ 0 + Adult*Upwelling*Management, cmc=F),
  #       nlf(zi ~ IntZi + AdultZi),
  #       IntZi ~ 1 + (1|Site),
  #       lf(AdultZi ~ 0 + Adult, cmc=F),
  #       nl=T),
  # m10=bf(Juvenile ~ IntAbund + FullAbund,
  #       IntAbund ~ 1 + (1|Site),
  #       lf(FullAbund ~ 0 + Adult*Upwelling*Management, cmc=F),
  #       nlf(zi ~ IntZi + FullZi),
  #       IntZi ~ 1 + (1|Site),
  #       lf(FullZi ~ 0 + Adult*Upwelling*Management, cmc=F),
  #       nl=T)
)
priors.r <- vector("list", length(form.r))
for(i in seq_along(priors.r)) {
  priors.r[[i]] <- c(prior(normal(0, 0.1), class="sd", lb=0, nlpar="IntAbund"),
                   prior(normal(0.5, 1), class="b", nlpar="IntAbund"),
                   prior(normal(0, 0.1), class="sd", lb=0, nlpar="IntZi"),
                   prior(normal(0, 1), class="b", nlpar="IntZi"))
  if(grepl("AdultAbund", paste0(form.r[[i]], collapse=" "))) {
    priors.r[[i]] <- c(priors.r[[i]], prior(normal(-1, 1), class="b", nlpar="AdultAbund", ub=0))
  }
  if(grepl("AdultZi", paste0(form.r[[i]], collapse=" "))) {
    priors.r[[i]] <- c(priors.r[[i]], prior(normal(1, 1), class="b", nlpar="AdultZi", ub=0))
  }
  if(grepl("MgmtAbund", paste0(form.r[[i]], collapse=" "))) {
    priors.r[[i]] <- c(priors.r[[i]], prior(normal(0, 1), class="b", nlpar="MgmtAbund"))
  }
  if(grepl("MgmtZi", paste0(form.r[[i]], collapse=" "))) {
    priors.r[[i]] <- c(priors.r[[i]], prior(normal(0, 1), class="b", nlpar="MgmtZi"))
  }
  if(grepl("UpwellAbund", paste0(form.r[[i]], collapse=" "))) {
    priors.r[[i]] <- c(priors.r[[i]], prior(normal(0, 1), class="b", nlpar="UpwellAbund"))
  }
  if(grepl("UpwellZi", paste0(form.r[[i]], collapse=" "))) {
    priors.r[[i]] <- c(priors.r[[i]], prior(normal(0, 1), class="b", nlpar="UpwellZi"))
  }
  if(grepl("FullAbund", paste0(form.r[[i]], collapse=" "))) {
    priors.r[[i]] <- c(priors.r[[i]], prior(normal(0, 1), class="b", nlpar="FullAbund"))
  }
  if(grepl("FullZi", paste0(form.r[[i]], collapse=" "))) {
    priors.r[[i]] <- c(priors.r[[i]], prior(normal(0, 1), class="b", nlpar="FullZi"))
  }
  if(grepl("FullInteractionAbund", paste0(form.r[[i]], collapse=" "))) {
    priors.r[[i]] <- c(priors.r[[i]], prior(normal(0, 0.5), class="b", nlpar="FullInteractionAbund"))
  }
  if(grepl("FullInteractionZi", paste0(form.r[[i]], collapse=" "))) {
    priors.r[[i]] <- c(priors.r[[i]], prior(normal(0, 0.5), class="b", nlpar="FullInteractionZi"))
  }
}
out.r <- map2(form.r, priors.r,
             ~brm(.x, data=rcrDens.df, cores=4, save_pars=save_pars(all=T), 
                  family=zero_inflated_negbinomial(), control=list(adapt_delta=0.9),
                  prior=.y))
out.r <- map(out.r, ~add_criterion(.x, "loo", moment_match=T, cores=cores))
loo.r <- map(out.r, ~loo(.x, cores=cores)) %>% loo_compare()
saveRDS(loo.r, "out/regr/loo_r.rds")

opt.r <- loo.r %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -2) %>%
  mutate(model2=paste0("m", str_pad(str_sub(model, 2, -1), 2, "left", "0"))) %>%
  arrange(model2)
opt2.r <- loo.r %>% as_tibble(rownames="model") %>% 
  arrange(desc(as.numeric(elpd_diff)))

saveRDS(out.r[[opt.r$model[1]]], "out/regr/opt_r.rds")
saveRDS(out.r[[opt2.r$model[1]]], "out/regr/opt2_r.rds")
saveRDS(out.r, "out/regr/all_r.rds")
saveRDS(out.r$m10, "out/regr/full_r.rds")

# out.r <- brm(bf(Juvenile ~ Adult + (1|Site),
#                 zi ~ Adult + (1|Site)), 
#              family=zero_inflated_poisson(), data=rcrDens.df, cores=4,
#              prior=c(prior(normal(0, 0.1), class="sd", lb=0, group="Site"),
#                      prior(normal(0.5, 1), class="Intercept"),
#                      prior(normal(-1, 1), class="b", ub=0),
#                      prior(normal(0, 0.1), class="sd", lb=0, group="Site", dpar="zi"),
#                      prior(normal(0, 1), class="Intercept", dpar="zi"),
#                      prior(normal(0.5, 0.5), class="b", dpar="zi")), 
#              file="temp/out_r_ZIP.rds")
# pp_check(out.r, ndraws=100) + scale_x_continuous(trans="log1p")


rcrSize.df <- Growth_Surv %>% 
  select(Site, Management, Upwelling, Patch, Monitoring, 
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
  m4=bf(logHoldfast ~ Management + (1+Management|Site)),
  m5=bf(logHoldfast ~ Upwelling + (1|Site)),
  m6=bf(logHoldfast ~ Upwelling + (1+Upwelling|Site)),
  m7=bf(logHoldfast ~ Upwelling + Management + (1|Site)),
  m8=bf(logHoldfast ~ Upwelling + Management + (1+Upwelling+Management|Site)),
  m9=bf(logHoldfast ~ Upwelling * Management + (1|Site)),
  m10=bf(logHoldfast ~ Upwelling * Management + (1+Upwelling*Management|Site))
)

out.r_z <- map(form.r_z, 
             ~brm(.x, data=rcrSize.df, cores=4, save_pars=save_pars(all=T), 
                  control=list(adapt_delta=0.9)))
out.r_z <- map(out.r_z, ~add_criterion(.x, "loo", moment_match=T, cores=cores))
loo.r_z <- map(out.r_z, ~loo(.x, cores=cores)) %>% loo_compare()
saveRDS(loo.r_z, "out/regr/loo_r_z.rds")

opt.r_z <- loo.r_z %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -2) %>%
  mutate(model2=paste0("m", str_pad(str_sub(model, 2, -1), 2, "left", "0"))) %>%
  arrange(model2)
opt2.r_z <- loo.r_z %>% as_tibble(rownames="model") %>% 
  arrange(desc(as.numeric(elpd_diff)))

saveRDS(out.r_z[[opt.r_z$model[1]]], "out/regr/opt_r_z.rds")
saveRDS(out.r_z[[opt2.r_z$model[1]]], "out/regr/opt2_r_z.rds")
saveRDS(out.r_z, "out/regr/all_r_z.rds")
saveRDS(out.r_z$m4, "out/regr/full_r_z.rds")

# out.r_z <- brm(logHoldfast ~ Management + (1|Site), 
#                data=rcrSize.df, cores=4,
#                prior=c(prior(normal(0, 0.5), class="sd", lb=0),
#                        prior(normal(0, 1), class="b")),
#                file="temp/out_r_z.rds")








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
         Site2=factor(Site2, labels=janitor::make_clean_names(levels(Site2))),
         Upwelling=if_else(Zone=="North-Center", "Y", "N")) #%>%
  # filter(!is.na(Upwelling))



# From holdfast diameter to biomass

# candidate models of increasing complexity
form.hf_wt <- list(
  m0=bf(logWeight ~ logHoldfast),
  m1=bf(logWeight ~ logHoldfast + (1|Site2)),
  m2=bf(logWeight ~ logHoldfast + Management + (1|Site2)),
  m3=bf(logWeight ~ logHoldfast * Management + (1|Site2)),
  # m4=bf(logWeight ~ logHoldfast + Upwelling + (1|Site2)),
  # m5=bf(logWeight ~ logHoldfast * Upwelling + (1|Site2)),
  m6=bf(logWeight ~ logHoldfast + (1+logHoldfast|Site2)),
  m7=bf(logWeight ~ logHoldfast + Management + (1 + logHoldfast + Management|Site2)),
  m8=bf(logWeight ~ logHoldfast * Management + (1 + logHoldfast * Management|Site2))#,
  # m9=bf(logWeight ~ logHoldfast + Upwelling + (1 + logHoldfast + Upwelling|Site2)),
  # m10=bf(logWeight ~ logHoldfast * Upwelling + (1 + logHoldfast * Upwelling|Site2)),
  # m11=bf(logWeight ~ logHoldfast + Management + Upwelling + (1 |Site2)),
  # m12=bf(logWeight ~ logHoldfast * Management * Upwelling + (1 |Site2))
)

out.hf_wt <- map(form.hf_wt, 
                  ~brm(.x, data=allom.df, cores=4, save_pars=save_pars(all=T), 
                       control=list(adapt_delta=0.9)))
out.hf_wt <- map(out.hf_wt, ~add_criterion(.x, "loo", moment_match=T, cores=cores))
loo.hf_wt <- map(out.hf_wt, ~loo(.x, cores=cores)) %>% loo_compare()
saveRDS(loo.hf_wt, "out/regr/loo_allom_hw.rds")

opt.hf_wt <- loo.hf_wt %>% as_tibble(rownames="model") %>% 
  filter(elpd_diff > -2) %>%
  mutate(model2=paste0("m", str_pad(str_sub(model, 2, -1), 2, "left", "0"))) %>%
  arrange(model2)
opt2.hf_wt <- loo.hf_wt %>% as_tibble(rownames="model") %>% 
  arrange(desc(as.numeric(elpd_diff)))

# No support for upwelling, so re-fit using all data, where m1 is best

saveRDS(out.hf_wt[[opt.hf_wt$model[1]]], "out/regr/opt_allom_hw.rds")
saveRDS(out.hf_wt[[opt2.hf_wt$model[1]]], "out/regr/opt2_allom_hw.rds")
saveRDS(out.hf_wt, "out/regr/all_allom_hw.rds")
saveRDS(out.hf_wt$m6, "out/regr/full_allom_hw.rds")



# # From holdfast diameter to total length
# 
# form.hf_len <- list(
#   m0=bf(logLength ~ logHoldfast),
#   m1=bf(logLength ~ logHoldfast + (1|Site2)),
#   m2=bf(logLength ~ logHoldfast + I(logHoldfast^2) + (1|Site2)),
#   m3=bf(logLength ~ logHoldfast + I(logHoldfast^2) + Management + (1|Site2)),
#   m4=bf(logLength ~ logHoldfast * Management + I(logHoldfast^2) * Management + (1|Site2)),
#   m4=bf(logLength ~ logHoldfast * Management + I(logHoldfast^2) * Management + ( 1+ logHoldfast * Management + I(logHoldfast^2) * Management|Site2)),
# )
# 
# out.hf_len <- map(form.hf_len, 
#                   ~brm(.x, data=allom.df, cores=4, save_pars=save_pars(all=T), 
#                        control=list(adapt_delta=0.9)))
# out.hf_len <- map(out.hf_len, ~add_criterion(.x, "loo", moment_match=T, cores=cores))
# loo.hf_len <- map(out.hf_len, ~loo(.x, cores=cores)) %>% loo_compare()
# saveRDS(loo.hf_len, "out/regr/loo_allom_hl.rds")
# 
# opt.hf_len <- loo.hf_len %>% as_tibble(rownames="model") %>% 
#   filter(elpd_diff > -4) %>%
#   arrange(model)
# 
# saveRDS(out.hf_len[[opt.hf_len$model[1]]], "out/regr/opt_allom_hl.rds")
# saveRDS(out.hf_len, "out/regr/all_allom_hl.rds")






# summaries ---------------------------------------------------------------

opt.f <- dir("out/regr", "opt_")
opt.mods <- map(opt.f, ~readRDS(glue("out/regr/{.x}"))) %>%
  setNames(str_sub(opt.f, 5, -5))
map(opt.mods, ~.x$formula)

opt2.f <- dir("out/regr", "opt2_")
opt2.mods <- map(opt2.f, ~readRDS(glue("out/regr/{.x}"))) %>%
  setNames(str_sub(opt2.f, 6, -5))
map(opt2.mods, ~.x$formula)
