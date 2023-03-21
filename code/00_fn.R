# KELPER2
# helper functions
# Tim Szewczyk



# Miscellaneous helper functions



# aliases -----------------------------------------------------------------

invlogit <- brms::inv_logit_scaled
logit <- brms::logit_scaled




# Generate a sequence with log-scaled spacing
seq_ln <- function(from=1, to=1, length.out=10) {
  ln_from <- log(from)
  ln_to <- log(to)
  exp(seq(log(from), log(to), length.out=length.out))
}





# IPM wrapper -------------------------------------------------------------


# !! deprecated !! use build_IPM
run_ipmr <- function(param.ls, mesh.ls, N.init) {
  init_ipm(sim_gen="simple",
           di_dd="dd",
           det_stoch="det") %>%
    define_kernel(
      name="P",
      formula=s * G * s_densEffect,
      family="CC",
      G=dnorm(z_2, mu_g, sd_g),
      mu_g=posterior_epred(grow_mod, re.form=re.form, ndraws=ndraws,
                           newdata=tibble(logHoldfast=z_1, Management=mgmt)) %>%
        colMeans(),
      s=posterior_epred(surv_mod, re.form=re.form, ndraws=ndraws, nlpar="pSurvYr",
                        newdata=tibble(logHoldfast=z_1, Management=mgmt, maxAge_yrs=1)) %>%
        colMeans() %>% brms::inv_logit_scaled(),
      s_densEffect=1,
      # s_densEffect=max(0, min(1, 1 - ((sum(n_z_t * pi * (exp(z_1)/2)^2) - cumsum(n_z_t * pi * (exp(z_1)/2)^2))/K)^0.5 ) ),
      # s_densEffect=max(0, min(1, (1 - (sum(n_z_t) - cumsum(n_z_t))/K ))),
      data_list=param.ls,
      states=list(c('z')),
      evict_cor=TRUE,
      evict_fun=truncated_distributions("norm", "G")
    ) %>% 
    define_kernel(
      name="F",
      formula=r_s * r_d,
      family="CC",
      r_s=posterior_epred(recr_mod, re.form=re.form, ndraws=ndraws,
                          newdata=tibble(Adult=sum(n_z_t*(z_1>log(10))), Management=mgmt)) %>%
        colMeans(),
      r_d=dnorm(z_2, recr_mu, recr_sd),
      data_list=param.ls,
      states=list(c("z")),
      evict_cor=TRUE,
      evict_fun=truncated_distributions("norm", "r_d")
    ) %>% 
    define_impl(
      make_impl_args_list(
        kernel_names=c("P", "F"),
        int_rule=rep('midpoint', 2),
        state_start=rep('z', 2),
        state_end=rep('z', 2)
      )
    ) %>% 
    define_domains(z=c(mesh.ls$L, mesh.ls$U, mesh.ls$n_mesh_p)) %>% 
    define_pop_state(n_z=N.init) %>% 
    make_ipm(iterate=TRUE, iterations=mesh.ls$n_yr)
}


# !! deprecated !! use build_IPM
run_ipmr_regMan <- function(param.ls, mesh.ls, N.init) {
  
  IPM <- init_ipm(sim_gen="simple",
                  di_dd="dd",
                  det_stoch="det") %>%
    define_kernel(
      name="P",
      formula=s * G,
      family="CC",
      G=dnorm(z_2, mu_g, sd_g),
      mu_g=coef.g[1] + coef.g[2]*z_1,
      s=plogis(coef.s[1] + coef.s[2]*z_1 + 
                 coef.s[3]*(mgmt=="TURF") + coef.s[4]*z_1*(mgmt=="TURF") +
                 s_dd),
      s_dd=s_DD * sqrt((sum(n_z_t * pi * (exp(z_1)/2)^2) - cumsum(n_z_t * pi * (exp(z_1)/2)^2))),
      coef.s=colMeans(fixef(surv_mod, summary=F)[sample.int(4000, ndraws),]),
      coef.g=colMeans(fixef(grow_mod, summary=F)[sample.int(4000, ndraws),]),
      data_list=param.ls,
      states=list(c('z')),
      evict_cor=TRUE,
      evict_fun=truncated_distributions("norm", "G")
    ) %>% 
    define_kernel(
      name="F",
      formula=r_s * r_d,
      family="CC",
      r_s=posterior_epred(recr_mod, re.form=re.form, ndraws=ndraws,
                          newdata=tibble(Adult=sum(n_z_t*(z_1>log(adultThresh))), Management=mgmt)) %>%
        colMeans(),
      r_d=dnorm(z_2, recr_mu, recr_sd),
      data_list=param.ls,
      states=list(c("z")),
      evict_cor=TRUE,
      evict_fun=truncated_distributions("norm", "r_d")
    ) %>% 
    define_impl(
      make_impl_args_list(
        kernel_names=c("P", "F"),
        int_rule=rep('midpoint', 2),
        state_start=rep('z', 2),
        state_end=rep('z', 2)
      )
    ) %>% 
    define_domains(z=c(mesh.ls$L, mesh.ls$U, mesh.ls$n_mesh_p)) %>% 
    define_pop_state(n_z=N.init) %>% 
    make_ipm(iterate=TRUE, iterations=mesh.ls$n_yr)
}


build_IPM <- function(par.ls, mesh.ls) {
  
  init_ipm(sim_gen="simple",
           di_dd="dd",
           det_stoch="det") %>%
    define_kernel(
      name="P",
      formula=s*G,
      family="CC",
      G=dnorm(z_2, mu_g, sd_g),
      mu_g=b.g[1] + b.g[2]*z_1,
      b.g=fixef(out_g, summary=F)[draw,],
      sd_g=c(as.matrix(out_g, variable="sigma", draw=draw)),
      s=plogis(b.s[1] + b.s[2]*z_1 + b.s[3]*(mgmt=="TURF") + b.s[4]*z_1*(mgmt=="TURF") + s_dd),
      b.s=fixef(out_s, summary=F)[draw,],
      # Asymmetric competition
      # s_DD * (Area[z > z_i])^pow
      # Allows for tuning the shape of the impact by size
      s_dd=s_DD * (sum(n_z_t * pi * (exp(z_1)/2)^2) - cumsum(n_z_t * pi * (exp(z_1)/2)^2))^s_DD_pow,
      data_list=par.ls,
      states=list(c('z')),
      evict_cor=TRUE,
      evict_fun=truncated_distributions("norm", "G")
    ) %>% 
    define_kernel(
      name="F",
      formula=r_s*r_d,
      family="CC",
      r_s_nRcrPred=c(posterior_predict(out_r, re.form=NA, draw_ids=draw,
                                       newdata=tibble(Adult=sum(n_z_t*(z_1>log(adultThresh))), 
                                                      Management=mgmt,
                                                      Upwelling=upwell))),
      r_s=min(rcr_max, r_s_nRcrPred),
      r_d=dnorm(z_2, mu_rcr, sd_rcr),
      mu_rcr=c(posterior_linpred(out_r_z, re.form=NA, draw_ids=draw, newdata=tibble(Management=mgmt))),
      sd_rcr=c(as.matrix(out_r_z, variable="sigma", draw=draw)),
      data_list=par.ls,
      states=list(c("z")),
      evict_cor=TRUE,
      evict_fun=truncated_distributions("norm", "r_d")
    ) %>% 
    define_impl(
      make_impl_args_list(
        kernel_names=c("P", "F"),
        int_rule=rep('midpoint', 2),
        state_start=rep('z', 2),
        state_end=rep('z', 2)
      )
    ) %>% 
    define_domains(z=c(mesh.ls$L, mesh.ls$U, mesh.ls$n_mesh_p))
}





build_IPM2 <- function(par.ls, mesh.ls) {
  
  init_ipm(sim_gen="simple",
           di_dd="dd",
           det_stoch="det") %>%
    define_kernel(
      name="P",
      formula=s*G,
      family="CC",
      G=dnorm(z_2, mu_g, sd_g),
      # mu_g=b.g[1] + b.g[2]*z_1 + b.g[3]*(mgmt=="TURF") + b.g[4]*(upwell=="Y"),
      mu_g=b.g[1] + b.g[2]*z_1 + b.g[3]*(upwell=="Y") + b.g[4]*z_1*(upwell=="Y"),
      b.g=fixef(out_g, summary=F)[draw,],
      sd_g=c(as.matrix(out_g, variable="sigma", draw=draw)),
      s=plogis(b.s[1] + b.s[2]*z_1 + b.s[3]*(mgmt=="TURF") + b.s[4]*(upwell=="Y") + 
                 b.s[5]*z_1*(mgmt=="TURF") + b.s[6]*z_1*(upwell=="Y") +
                 b.s[7]*(mgmt=="TURF")*(upwell=="Y") + 
                 b.s[8]*z_1*(mgmt=="TURF")*(upwell=="Y") +
                 s_dd),
      b.s=fixef(out_s, summary=F)[draw,],
      # Asymmetric competition
      # s_DD * (Area[z > z_i])^pow
      # Allows for tuning the shape of the impact by size
      s_dd=s_DD * (sum(n_z_t * pi * (exp(z_1)/2)^2) - cumsum(n_z_t * pi * (exp(z_1)/2)^2))^s_DD_pow,
      data_list=par.ls,
      states=list(c('z')),
      evict_cor=TRUE,
      evict_fun=truncated_distributions("norm", "G")
    ) %>% 
    define_kernel(
      name="F",
      formula=r_s*r_d,
      family="CC",
      r_s=0,
      r_d=dnorm(z_2, 0, 0.1),
      data_list=par.ls,
      states=list(c("z")),
      evict_cor=TRUE,
      evict_fun=truncated_distributions("norm", "r_d")
    ) %>% 
    define_impl(
      make_impl_args_list(
        kernel_names=c("P", "F"),
        int_rule=rep('midpoint', 2),
        state_start=rep('z', 2),
        state_end=rep('z', 2)
      )
    ) %>% 
    define_domains(z=c(mesh.ls$L, mesh.ls$U, mesh.ls$n_mesh_p))
}




calc_new_recruits <- function(ipm, pars, z_2) {
  nRcr_k <- min(pars$rcr_max, 
                posterior_predict(
                  pars$out_r, re.form=NA, draw_ids=pars$draw,
                  newdata=tibble(Adult=sum(ipm$pop_state$n_z[adult_i,1]), 
                                 Management=pars$mgmt,
                                 Upwelling=pars$upwell)))
  mu_rcr <- c(posterior_linpred(pars$out_r_z, re.form=NA, draw_ids=pars$draw, 
                                newdata=tibble(Management=pars$mgmt)))
  sd_rcr <- c(as.matrix(pars$out_r_z, variable="sigma", draw=pars$draw))
  r_d <- dnorm(z_2, mu_rcr, sd_rcr)
  return(ipm$pop_state$n_z[,2] + nRcr_k*r_d/sum(r_d))
}



# IPM analysis ------------------------------------------------------------


