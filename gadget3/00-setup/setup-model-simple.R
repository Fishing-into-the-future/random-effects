## -----------------------------------------------------------------------------
##
## Runner to set up stocks parameters and actions using g3 defaults as far as possible
##
## -----------------------------------------------------------------------------

## Setup formulas that depend on setup options

## Two internal functions (from gadgetutils) that should not be made visible
# (1) ogive for splitting ages between stocks
g3a_initial_ageprop <- function(alpha = gadget3::g3_parameterized('mat_initial_alpha', by_stock = 'species'), 
                                a50 = gadget3::g3_parameterized('mat_initial_a50', by_stock = 'species'),
                                by_stock = FALSE){  
  gadget3:::f_substitute(~bounded(-1 * alpha * (age - a50), 0, 1),
                         list(alpha = alpha, 
                              a50 = a50)
  )
}

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

initvonb2 <- 
  gadget3::g3a_renewal_vonb_t0(K = g3_parameterized('K', 
                                                    by_stock = stocks, 
                                                    by_year = timevarying_K),
                                 by_stock = stocks)


initvonb <- 
  gadget3::g3a_renewal_vonb_recl(K = g3_parameterized('K', 
                                                    by_stock = stocks, 
                                                    by_year = timevarying_K),
                               by_stock = stocks,
                               recage = 3)


## Initial numbers
# One parameter per age group per stock
initabun <- gadget3::g3a_renewal_initabund(by_stock = TRUE, by_stock_f = stocks, recage = gadget3::g3_stock_def(imm_stock, 'minage'))

# One parameter per age group across stocks
if (initial_abund_mode == 2){
  initabun <- gadget3::g3a_renewal_initabund(by_stock = stocks, proportion_f = g3a_initial_ageprop(), recage = gadget3::g3_stock_def(imm_stock, 'minage'))
}

## Initial sd's
# One per age group per stock
initsd <- gadget3::g3_parameterized('init.sd', by_stock = TRUE, by_age = TRUE)

# Parameteric model
if (init_sd_parametric){
  initsd <- g3a_initial_sigma(mean_l = initvonb)
}


## -----------------------------------------------------------------------------
## Setup model actions
## -----------------------------------------------------------------------------

## IMMATURE ACTIONS

imm_actions <- 
  list(
    
    ## INITIAL CONDITIONS
    g3a_initialconditions_normalparam(imm_stock,
                                      factor_f = initabun,
                                      mean_f = initvonb,
                                      stddev_f = initsd),
    
    ## NATURAL MORTALITY
    gadget3::g3a_naturalmortality(imm_stock, gadget3::g3a_naturalmortality_exp()),

    ## AGEING
    gadget3::g3a_age(imm_stock, output_stocks = list(mat_stock)),

    ## GROWTH AND MATURITY
    gadget3::g3a_growmature(imm_stock,

                            ## Growth
                            impl_f = gadget3::g3a_grow_impl_bbinom(
                              delta_len_f = gadget3::g3a_grow_lengthvbsimple(kappa_f = g3_parameterized('K',
                                                                                                        by_stock = stocks,
                                                                                                        by_year = timevarying_K),
                                                                             by_stock = stocks),
                              delta_wgt_f = gadget3::g3a_grow_weightsimple(),
                              maxlengthgroupgrowth = maxlengthgroupgrowth,
                              by_stock = stocks ## Ensures the bbin parameter across stocks within species
                            ),

                            ## Maturity
                            maturity_f = gadget3::g3a_mature_continuous(),
                            output_stocks = list(mat_stock),
                            transition_f = ~TRUE
    ),

    # RENEWAL
    gadget3::g3a_renewal_normalparam(imm_stock, mean_f = initvonb),
    
    list()
    
  )

## MATURE ACTIONS

mat_actions <- 
  list(
    
    ## INITIAL CONDITIONS
    g3a_initialconditions_normalparam(mat_stock,
                                      factor_f = initabun,
                                      mean_f = initvonb,
                                      stddev_f = initsd),
    
    ## NATURAL MORTALITY
    gadget3::g3a_naturalmortality(mat_stock, gadget3::g3a_naturalmortality_exp()),

    ## AGEING
    gadget3::g3a_age(mat_stock),

    # GROWTH
    gadget3::g3a_growmature(mat_stock,
                            impl_f = gadget3::g3a_grow_impl_bbinom(
                              delta_len_f =
                                gadget3::g3a_grow_lengthvbsimple(kappa_f = g3_parameterized('K',
                                                                                            by_stock = stocks,
                                                                                            by_year = timevarying_K),
                                                                 by_stock = stocks),
                              delta_wgt_f = gadget3::g3a_grow_weightsimple(),
                              maxlengthgroupgrowth = maxlengthgroupgrowth,
                              by_stock = stocks)
    ),
    
    list()
    
  )

## Compile stock actions
stock_actions <- c(imm_actions, mat_actions, list())

