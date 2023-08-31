## -----------------------------------------------------------------------------
##
## Runner to set up stocks parameters and actions
##
## -----------------------------------------------------------------------------

## Common id for parameters that correspond to all stocks
comp_id <- 'species'

## -----------------------------------------------------------------------------
## Parameters
## -----------------------------------------------------------------------------

model_pars <- 
  list(
   
    ## Weight-length
    walpha = gadget3::g3_parameterized('walpha', by_stock = comp_id),
    wbeta = gadget3::g3_parameterized('wbeta', by_stock = comp_id),
    
    ## Growth
    Linf = gadget3::g3_parameterized('Linf', by_stock = comp_id, by_year = FALSE),
    K = gadget3::g3_parameterized('K', by_stock = comp_id, by_year = timevarying_K),
    bbin = gadget3::g3_parameterized('bbin', by_stock = comp_id),
    maxlengthgroupgrowth = maxlengthgroupgrowth,
    
    ## Maturity
    mat_alpha = gadget3::g3_parameterized('mat_alpha', by_stock = comp_id, by_year = FALSE),
    mat_l50 = gadget3::g3_parameterized('mat_l50', by_stock = comp_id, by_year = FALSE),
    
    ## Recruitment
    renewal = gadget3::g3_parameterized(name = 'rec',
                                        by_stock = stocks,
                                        by_year = TRUE,
                                        scale = 
                                          gadget3::g3_parameterized(name = 'rec.scalar',
                                                                    by_stock = stocks),
                                        ifmissing = NaN),
    
    recl = gadget3::g3_parameterized('recl', by_stock = stocks),
    recsd = gadget3::g3_parameterized('recl', 
                                      by_stock = stocks, 
                                      scale = 
                                        gadget3::g3_parameterized('reclcv', by_stock = stocks)),
    
    ## Initial conditions
    initsd = gadget3::g3_parameterized('init.sd', by_stock = stocks, by_age = TRUE),
    
    ## Natural mortality
    natm = gadget3::g3_parameterized('M', by_stock = FALSE, by_age = FALSE)
    
    
    
    
  )

## Initial vonb
initvonb <- gadget3::g3a_renewal_vonb(Linf = model_pars$Linf,
                                      K = model_pars$K,
                                      recl = model_pars$recl,
                                      recage = gadget3::g3_stock_def(imm_stock, 'minage'))

## -----------------------------------------------------------------------------
## Setup model actions
## -----------------------------------------------------------------------------

## IMMATURE ACTIONS

imm_actions <- 
  list(
    
    ## INITIAL CONDITIONS
    g3a_initialconditions_normalparam(imm_stock,
                                      factor_f = 
                                        gadgetutils::init_abund(imm_stock,
                                                                mat_stock,
                                                                comp_id,
                                                                mature = FALSE,
                                                                init_mode = initial_abund_mode,
                                                                naturalmortality = model_pars$natm),
                                      mean_f = initvonb,
                                      stddev_f = gadgetutils::init_sd(stock = imm_stock, 
                                                                      id = comp_id, 
                                                                      parametric = TRUE,
                                                                      mean_len = initvonb),
                                      alpha_f = model_pars$walpha,
                                      beta_f = model_pars$wbeta),
    
    ## NATURAL MORTALITY
    gadget3::g3a_naturalmortality(imm_stock, 
                                  gadget3::g3a_naturalmortality_exp(param_f = model_pars$natm)),
    
    ## AGEING
    gadget3::g3a_age(imm_stock, output_stocks = list(mat_stock)),
    
    ## GROWTH AND MATURITY
    gadget3::g3a_growmature(imm_stock,
                            
                            ## Growth
                            impl_f = gadget3::g3a_grow_impl_bbinom(
                              delta_len_f =
                                gadget3::g3a_grow_lengthvbsimple(
                                  linf_f = model_pars$Linf,
                                  kappa_f = model_pars$K
                                ),
                              delta_wgt_f = gadget3::g3a_grow_weightsimple(
                                alpha_f = model_pars$walpha,
                                beta_f = model_pars$wbeta
                              ),
                              beta_f = model_pars$bbin,
                              maxlengthgroupgrowth = model_pars$maxlengthgroupgrowth
                              ),
                            
                            ## Maturity
                            maturity_f =
                              gadget3::g3a_mature_continuous(
                                alpha = model_pars$mat_alpha,
                                l50 = model_pars$mat_l50
                              ),
                            output_stocks = list(mat_stock),
                            transition_f = TRUE
                            ),
    
    # RENEWAL
    gadget3::g3a_renewal_normalparam(imm_stock,
                                     factor_f = model_pars$renewal,
                                     mean_f = initvonb,
                                     stddev_f = model_pars$recsd,
                                     alpha_f = model_pars$walpha,
                                     beta_f = model_pars$wbeta,
                                     run_f = gadget3:::f_substitute(
                                       ~cur_step == 1 &&
                                         age == minage &&
                                         cur_time > 0 &&
                                         !cur_year_projection,
                                       list(minage = gadget3::g3_stock_def(imm_stock, 'minage')))),
    
    list()
    
  )

## MATURE ACTIONS

mat_actions <- 
  list(
    
    ## INITIAL CONDITIONS
    g3a_initialconditions_normalparam(mat_stock,
                                      factor_f = 
                                        gadgetutils::init_abund(imm_stock,
                                                                mat_stock,
                                                                comp_id,
                                                                mature = TRUE,
                                                                init_mode = initial_abund_mode,
                                                                naturalmortality = model_pars$natm),
                                      mean_f = initvonb,
                                      stddev_f = gadgetutils::init_sd(stock = mat_stock, 
                                                                      id = comp_id, 
                                                                      parametric = TRUE,
                                                                      mean_len = initvonb),
                                      alpha_f = model_pars$walpha,
                                      beta_f = model_pars$wbeta),
    
    ## NATURAL MORTALITY
    gadget3::g3a_naturalmortality(mat_stock, 
                                  gadget3::g3a_naturalmortality_exp(param_f = model_pars$natm)),
    
    ## AGEING
    gadget3::g3a_age(mat_stock, output_stocks = list()),
    
    ## GROWTH
    gadget3::g3a_growmature(mat_stock,
                            impl_f = gadget3::g3a_grow_impl_bbinom(
                              delta_len_f =
                                gadget3::g3a_grow_lengthvbsimple(
                                  linf_f = model_pars$Linf,
                                  kappa_f = model_pars$K
                                  ),
                              delta_wgt_f = gadget3::g3a_grow_weightsimple(
                                alpha_f = model_pars$walpha,
                                beta_f = model_pars$wbeta
                                ),
                              beta_f = model_pars$bbin,
                              maxlengthgroupgrowth = model_pars$maxlengthgroupgrowth
                              )
                            ),

    list()
    
  )

## Compile stock actions
stock_actions <- c(imm_actions, mat_actions, list())

