## -----------------------------------------------------------------------------
##
## Runner to build a Gadget3 model for Blue Ling 
##
## -----------------------------------------------------------------------------

library(mfdb)
library(gadget3)
library(gadgetutils)
library(gadgetplots)
library(tidyverse)
library(g3experiments)

## Local version of g3l_bound_penalty that allows bounding of random effects
## bounding random effects crashes model atm
#source("exploratory_models/gadget3/src/g3l_bounds_penalty.R")

## Model directory
base_dir <- 'exploratory_models/gadget3'

## Model version
vers <- 'models/04-ranexpclamprec_newton50it_tol10'

## -----------------------------------------------------------------------------
## OPTIONS 
## -----------------------------------------------------------------------------

## Whether or not to call the setup-data scripts
read_data <- FALSE

## Different options for running g3 and its diagnostics
run_iterative <- FALSE
run_jitter <- FALSE
run_retro <- FALSE
run_leaveout <- FALSE
run_bootstrap <- FALSE

## Parameter option:
include_bound_penalty <- TRUE
init_rec_scalar <- FALSE       # single scalar for initial conditions and recruitment?
timevarying_K <- FALSE        # time-varying K (vonb growth parameter)
exponentiate_bbin <- TRUE     # exponentiate the beta-binomial parameter
exponentiate_recruitment <- TRUE   # exponentiate recruitment parameters (necessary to include bounds for recruitment as a random effect)

## Iterative re-weighting:
cv_floor <- 0
iter_group_SI <- FALSE         # whether or not to group all SI's together for optimisation
iter_group_aut <- FALSE        # group together autumn likelihood components

## Should age readings from bottom trawl be included?
bmt_age <- FALSE
stratified <- FALSE
slope_early_SIs <- TRUE

## Model options:
single_stock_model <- TRUE 
single_fleet <- TRUE            # Single commercial fleet?
dome_comm <- FALSE             # Only applies if single_fleet == TRUE
dome_bmt <- FALSE             # Only applies if single_fleet == FALSE
dome_lln <- FALSE             # Only applies if single_fleet == FALSE

## Random effects
random_recruitment <- TRUE
random_initial <- FALSE
penalise_recruitment <- 0
penalise_initial <- 0
#bound_random_effects <- TRUE   
wgts_vers <- '04-BASELINE-single_stock_expbbin_tol10it3'

fix_initial_conditions <- TRUE


maxlengthgroupgrowth <- 5  # Maximum length group growth
lencv <- 0.1               # CV for initial conditions standard deviations

## Initial abundance mode
initial_abund_mode <- 2
init_sd_parametric <- FALSE

## Optimisation controls:
usenlminb <- FALSE
nlminb_control <- list(trace=1, eval.max=2000, iter.max=1000, rel.tol=1e-10)
  
# For g3_optim, BFGS in iterative, retro, jitter, leaveout
control_list <- list(maxit = 3000, reltol = 1e-10)#.Machine$double.eps^2)
random_control <- list(maxit = 1000, reltol = 1e-10)#.Machine$double.eps^2)#1e-10)

# TMB::newton optimisations
newton_control <- list(maxit = 50, tol = 1e-10, tol10 = 0)


## -----------------------------------------------------------------------------

tyr <- lubridate::year(Sys.Date())

year_range <- { if (bmt_age) 1970 else 1982 }:(tyr - 1)
#year_range <- 2000:2010
peel <- 0

species_name <- 'bli' 

## -----------------------------------------------------------------------------

reitmapping <- 
  read.table(
    system.file("demo-data", "reitmapping.tsv", package="mfdb"),
    header=TRUE,
    as.is=TRUE)

defaults <- list(
  area = mfdb_group("1" = unique(reitmapping$SUBDIVISION)),
  timestep = mfdb_timestep_quarterly,
  year = year_range,
  species = toupper(species_name))

# Map area names to integer area numbers (in this case only "1" ==> 1, but could be anything)
areas <- structure(
    seq_along(defaults$area),
    names = names(defaults$area))

# Timekeeping for the model, i.e. how long we run for
time_actions <- list(
  g3a_time(start_year = min(defaults$year), 
           end_year = max(defaults$year),
           defaults$timestep),
  list())

## Data and model folders
fs::dir_create(file.path(base_dir, c('data')))
fs::dir_create(file.path(base_dir, vers))
fs::dir_create(file.path(base_dir, vers, c('RETRO', 'WGTS', 'JITTER')))

## Read in a params_final to estimate weights init consitions
if (fix_initial_conditions | random_initial | random_recruitment){
  load(file.path(base_dir, 'models', wgts_vers, 'WGTS/params_final.Rdata'))
}

## ------------------------------------------------------------------------------------

source(file.path(base_dir, '00-setup', 'setup-stocks.R'))  # Generates stock objects

if (single_stock_model){
  source(file.path(base_dir, '00-setup', 'setup-model-simple-single-stock.R'))  # Generates mat_stock_actions / imm_stock_actions  
}else{
  source(file.path(base_dir, '00-setup', 'setup-model-simple.R'))  # Generates mat_stock_actions / imm_stock_actions  
}


## Load data objects ----------------------------------------------------------
if(read_data){
  mdb <- mfdb('Iceland', db_params = list(host = 'mfdb.hafro.is'))
  #mdb <- mfdb("../../mfdb/copy/iceland.duckdb")
  source(file.path(base_dir, '00-setup', 'setup-fleet-data.R'))
  source(file.path(base_dir, '00-setup', 'setup-catchdistribution.R'))
  source(file.path(base_dir, '00-setup', 'setup-indices.R'))
  source(file.path(base_dir, '00-setup', 'setup-initial_parameters.R'))
} else {
  fs::dir_ls(file.path(base_dir, 'data')) %>%
    stringr::str_subset('.Rdata') %>%
    lapply(load,.GlobalEnv)
}

## Configure model actions ------------------------------------------------------------


source(file.path(base_dir, '00-setup', 'setup-fleets.R'))  # Generates fleet_actions
source(file.path(base_dir, '00-setup', 'setup-likelihood.R'))  # Generates likelihood_actions
source(file.path(base_dir, '00-setup', 'setup-randomeffects.R'))  # Generates random actions

##### Compile the r- and tmb-based models ######################################


## Collate actions
actions <- c(
  stock_actions,
  fleet_actions,
  likelihood_actions,
  if (random_recruitment || random_initial) random_actions else NULL,
  time_actions,
  NULL
  )

# Turn actions into an R function
model <- g3_to_r(actions)#, strict = TRUE, trace = TRUE)

# Turn actions into C++ objective function code
tmb_model <- g3_to_tmb(actions)

## Fill in the parameter template
tmb_param <- 
  attr(tmb_model, 'parameter_template') %>% 
  g3_init_guess('\\.rec', 100, 0.001, 200, ifelse(random_recruitment, penalise_recruitment, 1)) %>% 
  g3_init_guess('\\.init.[0-9]', 100, 0.001, 200, ifelse(random_initial, penalise_initial, 1)) %>% 
  g3_init_guess('recl', 30, 10, 50, 1) %>% 
  g3_init_guess('reclcv', lencv, lencv*0.9, lencv*1.1, 0) %>% 
  g3_init_guess('lencv', 0.2, 0.1, 0.3, 0) %>% 
  g3_init_guess('rec.sd', 5, 4, 20, 1) %>% 
  g3_init_guess('init.scalar', 50, 1, 100, 1) %>% 
  g3_init_guess('rec.scalar', 5, 1, 10, 1) %>% 
  g3_init_guess('init.rec.scalar', 50, 1, 100, 1) %>% 
  #g3_init_guess('\\.rec.sigma', 0.2, -1, 1, 0) %>% 
  g3_init_guess('Linf', 140, 100, 200, 1) %>% 
  g3_init_guess('\\.K', 0.09, 0.04, 0.2, 1) %>%
  g3_init_guess('\\.t0', 1, -1, 5, 0) %>%
  g3_init_guess('bbin', 100, 1e-05, 1000, 1) %>% 
  g3_init_guess('\\.alpha', 0.5, 0.01, 3, 1) %>% 
  g3_init_guess('\\.l50', 50, 10, 100, 1) %>% 
  g3_init_guess('init.F', 0, 0.1, 1, 0) %>% 
  g3_init_guess('\\.M$', 0.15, 0.001, 1, 0) %>% 
  g3_init_guess('^M$', 0.15, 0.001, 1, 0) %>%  
  g3_init_guess('mat_initial_alpha', 1, 0.5, 2, 1) %>% 
  g3_init_guess('mat_initial_a50', 7, 3, 15, 0) %>% 
  g3_init_guess('mat.alpha', 0.07, 0.01, 0.2, 1) %>%
  g3_init_guess('mat.l50', mat.l50$l50, 0.75*mat.l50$l50, 1.25*mat.l50$l50, 1) %>%
  g3_init_guess('sigma_alpha', init.sigma.coef[['alpha']]) %>%
  g3_init_guess('sigma_beta', init.sigma.coef[['beta']]) %>%
  g3_init_guess('sigma_gamma', init.sigma.coef[['gamma']]) %>% 
  
  ## Fleet selection parameters
  g3_init_guess('andersen.p0$', 0, NA, NA, 0) %>%
  g3_init_guess('andersen.p2$', 1, NA, NA, 0) %>%
  g3_init_guess('andersen.L$', max(unlist(lapply(stocks, gadget3::g3_stock_def, 'minlen'))), NA, NA, 0) %>%
  g3_init_guess('\\.p1$', log(2), 0, 3, 1) %>%
  g3_init_guess('\\.p3$', 0.1, 0, 10, 1) %>%
  g3_init_guess('\\.p4$', 0.02, 0, 1e3, 1) %>%
  
  ## Random effects/penalities
  g3_init_guess('recruitment_sigma', 0.2, 0.01, 10, ifelse(penalise_recruitment == 0, 1, 0)) %>% 
  g3_init_guess('initial_sigma', 0.2, 0.01, 10, ifelse(penalise_initial == 0, 1, 0)) %>% 
  g3_init_guess('rnd_recruitment_weight', 1, 0.01, 100, 0) %>% 
  g3_init_guess('rnd_initial_weight', 1, 0.01, 100, 0) %>% 
  g3_init_guess('zero', 0, -1, 1, 0) %>% 
  
  ## Weight-length
  g3_init_guess('walpha', lw.constants$estimate[1]) %>% 
  g3_init_guess('wbeta', lw.constants$estimate[2])

## Create narrower bounds for time-varying parameters
if (timevarying_K) tmb_param <- g3_init_guess(tmb_param, '\\.K', 0.12, 0.1, 0.14, 1)


## Initial sd's
if (any(grepl('\\.init\\.sd', tmb_param$switch))){
  
  if (!single_stock_model){
    tmb_param[grepl('mat\\.init\\.sd', tmb_param$switch), 'value'] <-
      init.sigma %>% filter(age %in%
                              gadget3::g3_stock_def(mat_stock, 'minage'):
                              gadget3::g3_stock_def(mat_stock, 'maxage')) %>% .$ms
    
    tmb_param[grepl('imm\\.init\\.sd', tmb_param$switch), 'value'] <-
      init.sigma %>% filter(age %in%
                              gadget3::g3_stock_def(imm_stock, 'minage'):
                              gadget3::g3_stock_def(imm_stock, 'maxage')) %>% .$ms  
  }else{
    tmb_param[grepl('\\.init\\.sd', tmb_param$switch), 'value'] <-
      init.sigma %>% filter(age %in%
                              gadget3::g3_stock_def(single_stock, 'minage'):
                              gadget3::g3_stock_def(single_stock, 'maxage')) %>% .$ms
  }
  
  ## Turn off optimisation
  tmb_param <-
    tmb_param %>%
    mutate(optimise = case_when(grepl('init.sd', switch) ~ FALSE,
                                grepl('.M.[\\.[0-9]', switch) ~ FALSE,
                                TRUE~optimise))

}

## Turn off the preliminary and terminal year random effects, still estimate the terminal one
if (random_recruitment){
  tmb_param[tmb_param$switch %in% c('bling.rec.1982', 'bling.rec_exp.1982'), 'random'] <- FALSE
  tmb_param[tmb_param$switch %in% c('bling.rec.2022', 'bling.rec_exp.2022'), c('optimise', 'random')] <- c(TRUE, FALSE)
}

if (fix_initial_conditions){
  tmb_param$value[paste0('bling.init.', 3:20)] <- params_final$value[paste0('bling.init.', 3:20)]
  tmb_param$value['bling.init.scalar'] <- params_final$value['bling.init.scalar']
  tmb_param[tmb_param$switch %in% paste0('bling.init.', 3:20), 'optimise'] <- FALSE
  tmb_param[tmb_param$switch %in% 'bling.init.scalar', 'optimise'] <- FALSE
}

## --------------------------------------------------------------------------

if (include_bound_penalty){
  actions <- c(actions, list(g3l_bounds_penalty(tmb_param %>% 
                                                  filter(!grepl('_sigma$|_sigma_exp$',switch)))))#, include_random = exponentiate_recruitment && bound_random_effects)))
  model <- g3_to_r(actions)
  tmb_model <- g3_to_tmb(actions)
}

## Run the R-model
result <- model(tmb_param$value)
result[[1]]

# List all available reports
print(names(attributes(result)))

## Write out parameters and both models
save(tmb_param, file = file.path(base_dir, vers, 'tmb_param.Rdata'))
save(model, file = file.path(base_dir, vers, 'r_model.Rdata'))
save(tmb_model, file = file.path(base_dir, vers, 'tmb_model.Rdata'))

## -----------------------------------------------------------------------------

if (TRUE){
  
  tmp <- params_final[grepl('weight$', params_final$switch), 'switch']
  tmp <- tmp[tmp %in% tmb_param$switch]
  tmb_param$value[tmp] <- params_final$value[tmp]
  rm(tmp)
  
  # Compile and generate TMB ADFun (see ?TMB::MakeADFun)
  obj.fun <- g3_tmb_adfun(tmb_model, 
                          jitter_params(tmb_param), 
                          inner.control = newton_control)
  
  obj.fun$env$tracepar <- TRUE
  
  if (usenlminb){
    out <- nlminb(start = obj.fun$par, 
                  objective = obj.fun$fn, 
                  gradient = obj.fun$gr, 
                  control = nlminb_control)
  }else{
    out <- optim(par = obj.fun$par, 
                 fn = obj.fun$fn, 
                 gr = obj.fun$gr,
                 method = 'BFGS',
                 control = c(random_control, 
                             list(parscale = g3_tmb_parscale(tmb_param)))
    )  
  }
  
  sdout <- TMB::sdreport(obj.fun, out$par)
  
  save(tmb_model, tmb_param, 
       obj.fun,
       out, sdout, file = file.path(base_dir, vers, 'randomout.Rdata'))
  
  ## Merge estimated parameters into template
  pars <- summary(sdout, 'all')[,1]
  tmbpars <- g3_tmb_relist(tmb_param, pars[match(names(g3_tmb_par(tmb_param)), names(pars))])
  modcpp <- g3_to_tmb(c(attr(tmb_model, 'actions'), list(g3a_report_detail(attr(tmb_model, 'actions')))))
  newpars <- attr(modcpp, 'parameter_template')
  newpars$value[names(tmbpars)] <- tmbpars
  newpars$value$report_detail <- 1L
  
  ## Fit
  fit <- g3_fit(g3_to_r(attr(modcpp, 'actions')), newpars)
  save(fit, file = file.path(base_dir, vers, 'fit.Rdata'))
  gadget_plots(fit, file.path(base_dir, vers, 'figs'), file_type = 'html')
  
}
  
## -----------------------------------------------------------------------------

## Setup grouping for iterative re-weighting 
if (iter_group_SI){
  
  grouping <- list(sind = c('log_si_aut_1',
                            'log_si_aut_2a',
                            'log_si_aut_2b',
                            'log_si_aut_3a',
                            'log_si_aut_3b',
                            'log_si_aut_3c',
                            'log_si_aut_3d'))
  if (!single_fleet){
    grouping <- c(grouping,
                  list(bmt = c('ldist_bmt', if (bmt_age) 'aldist_bmt' else NULL),
                       aut = c('ldist_aut', if (iter_group_aut) 'matp_aut' else NULL)))
  }else{
    grouping <- c(grouping, list(comm = c('ldist_comm')))
  }
}else{
  
  grouping <- list(sind1 = c('log_si_aut_1',
                             'log_si_aut_2a',
                             'log_si_aut_2b'),
                   sind2 = c('log_si_aut_3a',
                             'log_si_aut_3b',
                             'log_si_aut_3c',
                             'log_si_aut_3d'))
  if (!single_fleet){
    grouping <- c(grouping,
                  list(bmt = c('ldist_bmt', if (bmt_age) 'aldist_bmt' else NULL),
                       aut = c('ldist_aut', if (iter_group_aut) 'matp_aut' else NULL)))
  }else{
    grouping <- c(grouping, list(comm = c('ldist_comm')))
  } 
}

## Iterative re-weighting
if (run_iterative){
 
  # Compile and generate TMB ADFun (see ?TMB::MakeADFun)
  obj.fun <- g3_tmb_adfun(tmb_model, tmb_param)
   
  ## Run iterative re-weighting
  params.out <- g3_iterative(file.path(base_dir, vers),
                             wgts = 'WGTS',
                             model = tmb_model,
                             params.in = tmb_param,
                             grouping = grouping,
                             method = 'BFGS',
                             control = control_list,
                             use_parscale = TRUE,
                             shortcut = FALSE,
                             cv_floor = cv_floor,
                             resume_final = FALSE)
  
  
  ## Get the model fit
  fit <- g3_fit(model, params.out)
  save(fit, file = file.path(base_dir, vers, 'fit.Rdata'))
  gadget_plots(fit, file.path(base_dir, vers, 'figs'), file_type = 'html')
  
  sdout <- TMB::sdreport(obj.fun, g3_tmb_par(params.out))
  save(sdout, file = file.path(base_dir, vers, 'sdout.Rdata'))
  
  
}else{
  load(file = file.path(base_dir, vers, 'fit.Rdata'))
}

## Run the retro
if (run_retro){
  
  npeels <- 6
  
  load(file = file.path(base_dir, vers, 'WGTS/params_final.Rdata'))
  
  tmp <- params_final$value[grepl('_weight$', params_final$switch)]
  retro_pars <- tmb_param
  retro_pars$value[names(tmp)] <- tmp
  
  retro_model <- list()
  retro_params <- list()
  
  for(peel in 1:npeels){
    
    source(file.path(base_dir, '00-setup', 'setup-likelihood.R'))  # Generates likelihood_actions
    
    retro_actions <- 
      c(stock_actions,
        fleet_actions,
        likelihood_actions,
        if (penalise_recruitment) random_actions else NULL,
        time_actions,
        list(g3l_bounds_penalty(tmb_param %>% 
                                  filter(!grepl('_sigma$|_sigma_exp$',switch, optimise))))# filter(optimise, !grepl('\\.init\\.[0-9]|\\.rec\\.[0-9]', switch)))))filter(optimise, !grepl('\\.init\\.[0-9]|\\.rec\\.[0-9]', switch))))
      )
    retro_model[[peel]] <- g3_to_tmb(retro_actions)
    retro_params[[peel]] <- params_final# retro_pars #gadgetutils::jitter_params(params_final)
    retro_params[[peel]]$value$retro_years <- peel
  } 
  
  peel <- 0
  
  retro <- 
    parallel::mclapply(1:npeels,function(x){
      g3_optim(retro_model[[x]],
               retro_params[[x]],
               control = control_list)
    }, 
    mc.cores = parallel::detectCores())
  ## Collate fit
  retro_fit <- 
    1:npeels %>%
    set_names(paste0('r',1:npeels)) %>% 
    purrr::map(function(x) g3_fit(model = retro_model[[x]], params = retro[[x]])) 
  
  save(retro_fit, file = file.path(base_dir, vers,'RETRO', 'retro_fit.Rdata'))
  
  #gadget_plots(fit, file.path(base_dir, vers, 'figs'), 'html', retrofit = retro_fit)
  
  
  # retro <- g3_retro(file.path(base_dir, vers),
  #                   tmb_model,
  #                   params.out,
  #                   num.years = 10)
}


## Iterative re-weighting step-by-step
if (FALSE){
  res1 <-  g3_lik_out(model, tmb_param) 
  res2 <-  g3_iterative_setup(res1, grouping = grouping)
  
  res3 <-  parallel::mclapply(res2$params, function(x) g3_iterative_run(x, 
                                                                        tmb_model), mc.cores = parallel::detectCores())
  res4 <-  parallel::mclapply(res3, function(x) g3_lik_out(model,x), mc.cores = parallel::detectCores())
  res5 <-  g3_iterative_final(res4)
  res6 <-  parallel::mclapply(res5, function(x) g3_iterative_run(x, tmb_model), mc.cores = parallel::detectCores())
  res7 <-  parallel::mclapply(res6, function(x) g3_lik_out(model,x), mc.cores = parallel::detectCores())
  
  lapply(res6, function(x) attr(x, 'summary')) %>% dplyr::bind_rows(.id = 'group')
  
}




