# KELPER2
# ipmr tests
# Tim Szewczyk



# IPM construction might be easier with ipmr, a recently publish R package



# setup -------------------------------------------------------------------

library(tidyverse); library(brms); library(ipmr)

sens <- function(ipm_obj, d_z) { 
  w <- right_ev(ipm_obj)[[1]] 
  v <- left_ev(ipm_obj)[[1]] 
  return( outer(v, w) / sum(v * w * d_z) ) 
}

elas <- function(ipm_obj, d_z) { 
  K <- make_iter_kernel(ipm_obj)$mega_matrix 
  sensitivity <- sens(ipm_obj, d_z) 
  lamb <- lambda(ipm_obj) 
  out <- sensitivity * (K / d_z) / lamb 
  return(out) 
}

R_nought <- function(ipm_obj) { 
  Pm <- ipm_obj$sub_kernels$P 
  Fm <- ipm_obj$sub_kernels$F 
  I <- diag(dim(Pm)[1]) 
  N <- solve(I - Pm) 
  R <- Fm %*% N 
  return( Re(eigen(R)$values)[1] ) 
}

gen_time <- function(ipm_obj) { 
  lamb <- unname(lambda(ipm_obj)) 
  r_nought <- R_nought(ipm_obj) 
  return(log(r_nought) / log(lamb)) 
}


# paper example -----------------------------------------------------------

data.ls <- list(
  s_i=-0.65,
  s_z=0.75,
  G_i=0.96,
  G_z=0.66,
  sd_G=0.67,
  mu_r=-0.08,
  sd_r=0.76,
  r_n_i=-1,
  r_n_z=0.3
)

test.ipm <- init_ipm("simple", "di", "det") %>%
  define_kernel(name="P",
                formula=surv*Grow,
                surv=plogis(s_i + s_z*z_1),
                Grow=dnorm(z_2, mu_G, sd_G),
                mu_G=G_i + G_z*z_1,
                data_list=data.ls,
                states=list(c("z"))) %>%
  define_kernel(name="F",
                formula=recr_number*recr_size,
                recr_number=exp(r_n_i + r_n_z*z_1),
                recr_size=dnorm(z_2, mu_r, sd_r),
                data_list=data.ls,
                states=list(c("z"))) %>%
  define_impl(
    list("P"=list(int_rule="midpoint", state_start="z", state_end="z"),
         "F"=list(int_rule="midpoint", state_start="z", state_end="z"))
  ) %>%
  define_domains(z=c(-2.65, 4.5, 250)) %>%
  define_pop_state(n_z=rep(1/250, 250)) %>%
  make_ipm()




sens <- function(ipm_obj, d_z) { 
  w <- right_ev(ipm_obj)[[1]] 
  v <- left_ev(ipm_obj)[[1]] 
  return( outer(v, w) / sum(v * w * d_z) ) 
}

elas <- function(ipm_obj, d_z) { 
  K <- make_iter_kernel(ipm_obj)$mega_matrix 
  sensitivity <- sens(ipm_obj, d_z) 
  lamb <- lambda(ipm_obj) 
  out <- sensitivity * (K / d_z) / lamb 
  return(out) 
}



R_nought <- function(ipm_obj) { 
  Pm <- ipm_obj$sub_kernels$P 
  Fm <- ipm_obj$sub_kernels$F 
  I <- diag(dim(Pm)[1]) 
  N <- solve(I - Pm) 
  R <- Fm %*% N 
  return( Re(eigen(R)$values)[1] ) 
}

gen_time <- function(ipm_obj) { 
  lamb <- unname(lambda(ipm_obj)) 
  r_nought <- R_nought(ipm_obj) 
  return(log(r_nought) / log(lamb)) 
}

mesh_info <- int_mesh(test.ipm)

sens_mat <- sens(test.ipm, mesh_info$d_z) 
elas_mat <- elas(test.ipm, mesh_info$d_z)
R0 <- R_nought(test.ipm) 
gen_T <- gen_time(test.ipm)



library(ggplot2) 
library(gridExtra) 
p_df <- ipm_to_df(test.ipm$sub_kernels$P) 
f_df <- ipm_to_df(test.ipm$sub_kernels$F) 
k_df <- ipm_to_df(K$mega_matrix)
sens_df <- ipm_to_df(sens_mat) 
elas_df <- ipm_to_df(elas_mat)


def_theme <-
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 16),
    axis.ticks = element_line(size = 1.5),
    axis.ticks.length = unit(0.08, "in"),
    axis.title.x = element_text(size = 20, margin = margin(
      t = 10,
      r = 0,
      l = 0,
      b = 2
    )),
    axis.title.y = element_text(size = 20, margin = margin(
      t = 0,
      r = 10,
      l = 2,
      b = 0
    )),
    legend.text = element_text(size = 16)
  )


p_plt <-
  ggplot(p_df) + geom_tile(aes(x = t, y = t_1, fill = value)) + geom_contour(
    aes(x = t, y = t_1, z = value),
    color = "black",
    size = 0.7,
    bins = 5
  ) + scale_fill_gradient("Value", low = "red", high = "yellow") + scale_x_continuous(name = "size (t)") + scale_y_continuous(name = "size (t + 1)") + def_theme + theme(legend.title = element_blank()) + ggtitle("P kernel")









# brms test ---------------------------------------------------------------


library(brms)


data(iceplant_ex) 
# grow_mod <- lm(log_size_next ~ log_size, data = iceplant_ex)
grow_mod <- brm(log_size_next ~ log_size, data=iceplant_ex, family=gaussian(), cores=4)
grow_sd <- sd(resid(grow_mod))
surv_mod <- brm(survival ~ log_size, data=iceplant_ex, family=bernoulli(), cores=4)
# surv_mod <- glm(survival ~ log_size, data = iceplant_ex, family = binomial())
repr_mod <- glm(repro ~ log_size, data = iceplant_ex, family = binomial())
flow_mod <- glm(flower_n ~ log_size, data = iceplant_ex, family = poisson())
recr_data <- subset(iceplant_ex, is.na(log_size))
recr_mu <- mean(recr_data$log_size_next) 
recr_sd <- sd(recr_data$log_size_next)
recr_n <- length(recr_data$log_size_next) 
flow_n <- sum(iceplant_ex$flower_n, na.rm = TRUE) 
recr_pr <- recr_n / flow_n

L <- min(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2 
U <- max(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2 
n_mesh_p <- 100


pred_par_list <-
  list(
    grow_mod = grow_mod,
    grow_sdv = grow_sd,
    surv_mod = surv_mod,
    repr_mod = repr_mod,
    flow_mod = flow_mod,
    recr_n = recr_n,
    flow_n = flow_n,
    recr_mu = recr_mu,
    recr_sd = recr_sd,
    recr_pr = recr_pr
  ) 
predict_method_carpobrotus <-
  init_ipm(sim_gen = "simple",
           di_dd = "di",
           det_stoch = "det") %>%
  define_kernel(
    name = "P",
    formula = s * G,
    family = "CC",
    G = dnorm(z_2, mu_g, grow_sdv),
    mu_g=colMeans(posterior_predict(grow_mod, newdata=tibble(log_size=z_1))),
    s=colMeans(posterior_predict(surv_mod, newdata=tibble(log_size=z_1))),
    data_list = pred_par_list,
    states = list(c('z')),
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "G")
  ) %>% 
  define_kernel(
    name = "F",
    formula = recr_pr * r_s * r_d * p_f,
    family = "CC",
    r_s = predict(flow_mod, newdata = data.frame(log_size = z_1), type = "response"),
    r_d = dnorm(z_2, recr_mu, recr_sd),
    p_f = predict(repr_mod, newdata = data.frame(log_size = z_1), type = "response"),
    data_list = pred_par_list,
    states = list(c("z")),
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "r_d")
  ) %>% 
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule = rep('midpoint', 2),
      state_start = rep('z', 2),
      state_end = rep('z', 2)
    )
  ) %>% 
  define_domains(z = c(L, U, n_mesh_p)) %>% 
  define_pop_state(n_z = rep(1 /100, n_mesh_p)) %>% 
  make_ipm(iterate = TRUE, iterations = 100)









# lessonia ----------------------------------------------------------------

library(tidyverse); library(glue); library(lubridate)
theme_set(theme_bw())

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
survey.i <- bind_rows(
  survey.i,
  survey.i %>% slice_head(n=1) %>%
    mutate(Monitoring=0,
           Date_Obs=Date_Harvest)
) %>%
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

surv.df <- Growth_Surv %>% filter(!is.na(doesSurvive))
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
             file="temp/out_s.rds")


grow.df <- Growth_Surv %>% 
  filter(doesSurvive) %>%
  mutate(annHoldfastGrowth=(nextHoldfast-Holdfast)/days_to_nextObs*365,
         annStipeGrowth=(nextStipes-Stipes)/days_to_nextObs*365,
         annLengthGrowth=(nextLength-Total_Length)/days_to_nextObs*365,
         nextHoldfastAnnual=Holdfast+annHoldfastGrowth,
         logNextHoldfastAnnual=log(nextHoldfastAnnual))
out.g <- brm(bf(logNextHoldfastAnnual ~ logHoldfast*Management + 
                  (1+logHoldfast*Management|Site/Patch/Plant)),
             family="gaussian", data=grow.df, cores=4,
             prior=c(prior(normal(0, 0.1), class="sd", lb=0, group="Site"),
                     prior(normal(0, 0.1), class="sd", lb=0, group="Site:Patch"),
                     prior(normal(0, 0.1), class="sd", lb=0, group="Site:Patch:Plant"),
                     prior(normal(0, 1), class="b"),
                     prior(normal(2, 1), class="Intercept")), 
             file="temp/out_g.rds")

rcrSize.df <- Growth_Surv %>% 
  select(Site, Management, Patch, Monitoring, 
         Plant, Holdfast, Stipes, Total_Length) %>%
  left_join(., survey.i) %>%
  arrange(Site, Patch, Plant, Monitoring) %>%
  mutate(prevSurv=Monitoring-1) %>%
  arrange(Site, Patch, Plant, Monitoring) %>%
  group_by(Site, Patch, Plant) %>%
  mutate(maxAge_days=cumsum(days_since_prevObs),
         maxAge_yrs=maxAge_days/365) %>%
  group_by(Site, Patch, Monitoring) %>%
  mutate(nPlants=sum(Holdfast>0),
         totHoldfast=sum(pi*(Holdfast/2)^2)) %>%
  ungroup %>%
  filter(Monitoring==1,
         Holdfast < 10)

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
  mutate(N=Juvenile+Adult)

out.r <- brm(bf(Juvenile ~ Adult + (1|Site),
                zi ~ Adult*Management + (1|Site)), 
             family=zero_inflated_negbinomial(), data=rcrDens.df, cores=4,
             prior=c(prior(normal(0, 0.1), class="sd", lb=0, group="Site"),
                     prior(normal(0.5, 1), class="Intercept"),
                     prior(normal(-1, 1), class="b", ub=0),
                     prior(normal(0, 0.1), class="sd", lb=0, group="Site", dpar="zi"),
                     prior(normal(0, 1), class="Intercept", dpar="zi"),
                     prior(normal(0.5, 0.5), class="b", dpar="zi")), 
             file="temp/out_r.rds")


L <- log(0.1) * 1.2 
U <- max(Growth_Surv$logHoldfast, na.rm = TRUE) * 1.2 
n_mesh_p <- 100
nSim <- 10
  
n_iter <- 20
pred_par_list <-  list(
  grow_mod = out.g,
  sd_g = sd(resid(out.g)),
  surv_mod = out.s,
  recr_mu = mean(log(rcrSize.df$Holdfast)),
  recr_sd = sd(log(rcrSize.df$Holdfast)),
  recr_mod = out.r,
  ndraws=10,
  re.form=NA,
  n_mesh_p=n_mesh_p
) 

for(i in 1:nSim) {
  OA.ipm <- init_ipm(sim_gen = "simple",
                     di_dd = "dd",
                     det_stoch = "det") %>%
    define_kernel(
      name = "P",
      formula = s * G,
      family = "CC",
      G = dnorm(z_2, mu_g, sd_g),
      mu_g=colMeans(posterior_epred(grow_mod, 
                                    newdata=tibble(logHoldfast=z_1,
                                                   Management=mgmt),
                                    re.form=re.form, ndraws=ndraws)),
      s=brms::inv_logit_scaled(colMeans(posterior_epred(surv_mod, 
                                                        newdata=tibble(logHoldfast=z_1, 
                                                                       Management=mgmt, 
                                                                       maxAge_yrs=1), 
                                                        re.form=re.form, ndraws=ndraws, nlpar="pSurvYr"))),
      data_list = c(pred_par_list, mgmt="OA"),
      states = list(c('z')),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "G")
    ) %>% 
    define_kernel(
      name = "F",
      formula = (r_s) * r_d,
      family = "CC",
      r_s=colMeans(posterior_epred(recr_mod,
                                   newdata=tibble(Adult=sum(n_z_t),
                                                  Management=mgmt),
                                   re.form=re.form, ndraws=ndraws)),
      r_d = dnorm(z_2, recr_mu, recr_sd),
      data_list = c(pred_par_list, mgmt="OA"),
      states = list(c("z")),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "r_d")
    ) %>% 
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "F"),
        int_rule = rep('midpoint', 2),
        state_start = rep('z', 2),
        state_end = rep('z', 2)
      )
    ) %>% 
    define_domains(z = c(L, U, n_mesh_p)) %>% 
    define_pop_state(n_z = c(1, rep(0, n_mesh_p-1))) %>% 
    make_ipm(iterate = TRUE, iterations = n_iter)
  
  TURF.ipm <- init_ipm(sim_gen = "simple",
                       di_dd = "dd",
                       det_stoch = "det") %>%
    define_kernel(
      name = "P",
      formula = s * G,
      family = "CC",
      G = dnorm(z_2, mu_g, sd_g),
      mu_g=colMeans(posterior_epred(grow_mod, 
                                    newdata=tibble(logHoldfast=z_1,
                                                   Management=mgmt),
                                    re.form=re.form, ndraws=ndraws)),
      s=brms::inv_logit_scaled(colMeans(posterior_epred(surv_mod, 
                                                        newdata=tibble(logHoldfast=z_1, 
                                                                       Management=mgmt, 
                                                                       maxAge_yrs=1), 
                                                        re.form=re.form, ndraws=ndraws, nlpar="pSurvYr"))),
      data_list = c(pred_par_list, mgmt="TURF"),
      states = list(c('z')),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "G")
    ) %>% 
    define_kernel(
      name = "F",
      formula = (r_s) * r_d,
      family = "CC",
      r_s=colMeans(posterior_epred(recr_mod,
                                   newdata=tibble(Adult=sum(n_z_t),
                                                  Management=mgmt),
                                   re.form=re.form, ndraws=ndraws)),
      r_d = dnorm(z_2, recr_mu, recr_sd),
      data_list = c(pred_par_list, mgmt="TURF"),
      states = list(c("z")),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "r_d")
    ) %>% 
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "F"),
        int_rule = rep('midpoint', 2),
        state_start = rep('z', 2),
        state_end = rep('z', 2)
      )
    ) %>% 
    define_domains(z = c(L, U, n_mesh_p)) %>% 
    define_pop_state(n_z = c(1, rep(0, n_mesh_p-1))) %>% 
    make_ipm(iterate = TRUE, iterations = n_iter)
  mesh_info <- int_mesh(OA.ipm)
  pop.df <- bind_rows(
    tibble(year=rep(1:(n_iter+1), each=n_mesh_p),
           size_ln=rep(mesh_info$z_1[1:n_mesh_p], times=n_iter+1),
           N=c(OA.ipm$pop_state$n_z),
           Management="OA"),
    tibble(year=rep(1:(n_iter+1), each=n_mesh_p),
           size_ln=rep(mesh_info$z_1[1:n_mesh_p], times=n_iter+1),
           N=c(TURF.ipm$pop_state$n_z),
           Management="TURF")) %>%
    mutate(size=exp(size_ln),
           stage=if_else(size>10, "Adult", "Juvenile"),
           sizeClass=case_when(size < 10 ~ "0-10",
                               between(size, 10, 20) ~ "10-20",
                               between(size, 20, 30) ~ "20-30",
                               between(size, 30, 40) ~ "30-40",
                               between(size, 40, 50) ~ "40-50",
                               size > 50 ~ "50+"),
           sim=i)
  saveRDS(pop.df, glue("temp/pop_{i}.rds"))
}




plot(lambda(OA.ipm), type="l", ylim=c(0,2.5), xlab="Year", ylab="lambda")
lines(1:(n_iter), lambda(TURF.ipm), col="steelblue")
abline(h=1, lty=3)


pop.df <- map_dfr(dir("temp", "pop_", full.names=T), readRDS)
ggplot(pop.df, aes(year, N, colour=size, group=size)) + 
  geom_line(alpha=0.5) + scale_colour_viridis_c() + facet_wrap(~Management)
pop.df %>%
  group_by(sim, Management, year, sizeClass) %>%
  summarise(N=sum(N)) %>%
  ggplot(aes(year, N, colour=sizeClass, group=paste(sim, sizeClass))) + 
  geom_line(size=1, alpha=0.5) + scale_colour_viridis_d() +
  facet_wrap(~Management)
pop.df %>%
  group_by(sim, Management, year, stage) %>%
  summarise(N=sum(N)) %>%
  ggplot(aes(year, N, colour=stage, group=paste(sim, stage))) + 
  geom_line(alpha=0.5) + 
  facet_wrap(~Management)

mesh_info <- int_mesh(lessonia.ipm)
sens_mat <- sens(lessonia.ipm, mesh_info$d_z) 
elas_mat <- elas(lessonia.ipm, mesh_info$d_z)
R0 <- R_nought(lessonia.ipm) 
gen_T <- gen_time(lessonia.ipm)


library(ggplot2) 
library(gridExtra) 
p_df <- ipm_to_df(lessonia.ipm$sub_kernels$P)
f_df <- ipm_to_df(lessonia.ipm$sub_kernels$F) 
K_df <- p_df %>% rename(P=value) %>% mutate(F=f_df$value, K=P+F)
# sens_df <- ipm_to_df(sens_mat)
# elas_df <- ipm_to_df(elas_mat)

ggplot(K_df, aes(x=t, y=t_1, fill=K)) + 
  geom_tile() + 
  scale_fill_gradient("Value", low="dodgerblue", high=scales::muted("red")) +
  scale_x_continuous(name = "size (t)") + 
  scale_y_continuous(name = "size (t + 1)") + 
  theme(legend.title = element_blank()) + ggtitle("K kernel")
