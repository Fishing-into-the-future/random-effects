## -----------------------------------------------------------------------------
##
## Runner to set up stocks parameters and actions using g3 defaults as far as possible
##
## -----------------------------------------------------------------------------

## Parametric relationship for initial sd's
g3a_initial_sigma <- function(alpha = gadget3::g3_parameterized('initial_sigma_alpha', by_stock = 'species'), 
                              beta = gadget3::g3_parameterized('initial_sigma_beta', by_stock = 'species'), 
                              gamma = gadget3::g3_parameterized('initial_sigma_gamma', by_stock = 'species'), 
                              mean_l = gadget3::g3a_renewal_vonb_t0(),
                              by_stock = FALSE){
  
  gadget3:::f_substitute(
    ~mean_l * ( alpha + beta/age + gamma * age),
    list(alpha = alpha,
         beta = beta,
         gamma = gamma,
         mean_l = mean_l)
  )
}

## Initial sd's
# One per age group per stock
initsd <- gadget3::g3_parameterized('init.sd', by_stock = TRUE, by_age = TRUE, value = 1)

# Parameteric model
if (init_sd_parametric){
  initsd <- g3a_initial_sigma(mean_l = initvonb)
}

if (fix_initial_conditions){
  mean_l <- gadget3::g3a_renewal_vonb_recl(Linf = params_final$value$bling.Linf,
                                          K = params_final$value$bling.K,
                                          recl = params_final$value$bling.recl,
                                          recage = gadget3::g3_stock_def(single_stock, 'minage'))
  
}else{
  mean_l <- gadget3::g3a_renewal_vonb_recl(K = gadget3::g3_parameterized('K', by_year = timevarying_K, by_stock = TRUE),
                                          recage = gadget3::g3_stock_def(single_stock, 'minage'))
}

initsd <- gadget3:::f_substitute(~mean_l * lencv,
                                 list(mean_l = mean_l, 
                                      lencv = g3_parameterized('lencv', by_stock = TRUE)))

## -----------------------------------------------------------------------------
## Setup model actions
## -----------------------------------------------------------------------------

stock_actions <- 
  list(
    ## Initial conditions
    g3a_initialconditions_normalparam(single_stock,
                                      factor_f = 
                                        gadget3::g3a_renewal_initabund(scalar = 
                                                                         gadget3::g3_parameterized(
                                                                           ifelse(init_rec_scalar, 'init.rec.scalar', 'init.scalar'),
                                                                           value = 1, 
                                                                           by_stock = TRUE),
                                                                       recage = gadget3::g3_stock_def(single_stock, 'minage')),
                                      mean_f = mean_l,
                                      stddev_f = initsd),
    ## Natural mortality
    gadget3::g3a_naturalmortality(single_stock, gadget3::g3a_naturalmortality_exp()),
    ## Ageing
    gadget3::g3a_age(single_stock),
    ## Growth 
    gadget3::g3a_growmature(single_stock,
                            ## Growth
                            impl_f = gadget3::g3a_grow_impl_bbinom(
                              delta_len_f = gadget3::g3a_grow_lengthvbsimple(
                                kappa_f = g3_parameterized('K', by_year = timevarying_K, by_stock = TRUE)
                                ),
                              delta_wgt_f = gadget3::g3a_grow_weightsimple(),
                              beta_f = gadget3::g3_parameterized('bbin', by_stock = TRUE, exponentiate = exponentiate_bbin),
                              maxlengthgroupgrowth = maxlengthgroupgrowth
                              )
                            ),
    # Renewal
    gadget3::g3a_renewal_normalparam(single_stock,
                                     factor_f = g3_parameterized('rec',
                                                                 by_stock = TRUE,
                                                                 by_year = TRUE,
                                                                 scale = g3_parameterized(
                                                                   name = ifelse(init_rec_scalar, 'init.rec.scalar', 'rec.scalar'),
                                                                   by_stock = TRUE),
                                                                 exponentiate = exponentiate_recruitment, 
                                                                 ifmissing = NaN),
                                     mean_f = gadget3::g3a_renewal_vonb_recl(K = gadget3::g3_parameterized('K', by_year = timevarying_K, by_stock = TRUE),
                                                                             recage = gadget3::g3_stock_def(single_stock, 'minage'))),
    
    list()
    
  )

