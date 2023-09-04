## Setup a random walk for recruitment...
random_actions <-
  c(
    list(
      if (random_recruitment){
        g3l_random_walk('recruitment',
                        substitute(x,
                                   list(x = 
                                          {if (!exponentiate_recruitment){
                                            substitute(log(avoid_zero(x2)),
                                                       list(x2 = 
                                                              g3_parameterized(paste0(
                                                                if (single_stock_model) gadgetutils::g3_stock_name(single_stock)
                                                                else gadgetutils::g3_stock_name(imm_stock, 'species'), '.rec',
                                                                ifelse(exponentiate_recruitment, '_exp', '')
                                                                ),
                                                                by_year = TRUE,
                                                                random = penalise_recruitment == 0)))
                                            }else{
                                              g3_parameterized(paste0(
                                                if (single_stock_model) gadgetutils::g3_stock_name(single_stock)
                                                else gadgetutils::g3_stock_name(imm_stock, 'species'), '.rec',
                                                ifelse(exponentiate_recruitment, '_exp', '')
                                                ),
                                                by_year = TRUE,
                                                random = penalise_recruitment == 0)
                                              }}
                                        )),
                        weight = g3_parameterized('rnd_recruitment_weight',scale = -1),
                        sigma_f = g3_parameterized('bling_recruitment_sigma'))  
        } else NULL  
  ),
  
  if (random_initial){
    gadget3::g3_stock_def(single_stock, 'minage'):gadget3::g3_stock_def(single_stock, 'maxage') %>% 
      map(function(age) 
        g3l_random_dnorm(sprintf('initial_condititions_age%s',age),
                         substitute(log(avoid_zero(x)), 
                                    list(x = 
                                           g3_parameterized(sprintf('bling.init.%s',age),
                                                            random = penalise_initial == 0))),
                         mean_f = g3_parameterized('zero', scale = 0),
                         weight = g3_parameterized('rnd_initial_weight', scale = -1),
                         sigma_f = g3_parameterized('bling_initial_sigma')))  
    } else NULL 
  ,
  list()
  )