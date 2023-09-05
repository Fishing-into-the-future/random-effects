## -------------------------------------------------------------------------------
##
## Runner to set up likelihoods
##
## -------------------------------------------------------------------------------

nll_breakdown <- TRUE  # Turn to TRUE to get per-step nll
lik_report <- TRUE

likelihood_actions <- list(
  g3l_understocking(stocks, nll_breakdown = nll_breakdown, weight = 100),
  
  if (single_fleet){
    
    g3l_catchdistribution(
      'ldist_comm',
      obs_data = ldist.comm[[boot_repl]],
      fleets = list(comm),
      stocks = stocks,
      g3l_distribution_sumofsquares(),
      nll_breakdown = nll_breakdown,
      report = lik_report)
    
  }else{
    
    c(
      g3l_catchdistribution(
      'ldist_lln',
      obs_data = ldist.lln[[boot_repl]],
      fleets = list(lln),
      stocks = stocks,
      g3l_distribution_sumofsquares(),
      nll_breakdown = nll_breakdown,
      report = lik_report),
    
    g3l_catchdistribution(
      'ldist_bmt',
      obs_data = ldist.bmt[[boot_repl]] %>%
        filter(!(year == 2020 & step == 2)),
      fleets = list(bmt),
      stocks = stocks,
      g3l_distribution_sumofsquares(),
      nll_breakdown = nll_breakdown,
      report = lik_report),
    
    if (bmt_age){
      g3l_catchdistribution(
        'aldist_bmt',
        (aldist.bmt[[boot_repl]]),
        fleets = list(bmt),
        stocks = stocks,
        g3l_distribution_sumofsquares(over = {if (stratified) c('area', 'length') else 'area' }),
        nll_breakdown = nll_breakdown,
        report = lik_report)
    }else NULL
    )
  },
  
  g3l_catchdistribution(
    'ldist_aut',
    obs_data = ldist.aut[[boot_repl]] %>% filter(year > 1999, year != 2011),
    fleets = list(aut),
    stocks = stocks,
    g3l_distribution_sumofsquares(),
    nll_breakdown = nll_breakdown,
    report = lik_report),
  
  if (!single_stock_model){
    
    g3l_catchdistribution(
      'matp_aut',
      (matp.aut[[boot_repl]] %>% filter(year > 1999, year != 2011) %>% 
         rename(stock = maturity_stage) %>%
         mutate(stock = recode(as.factor(stock), blingimm = 'bling_imm', blingmat = 'bling_mat'))),
      fleets = list(aut),
      stocks = stocks,
      g3l_distribution_sumofsquares(),
      nll_breakdown = nll_breakdown,
      report = lik_report) 
  },
  
  g3l_abundancedistribution(
    'si_aut_1',
    (aut.SI[[boot_repl]]$len20 %>% filter(year > 1999, year != 2011)) %>% filter(year < tyr - peel),
    fleets = list(),
    stocks = stocks,
    g3l_distribution_surveyindices_log(beta = if (slope_early_SIs) NULL else 1), 
    nll_breakdown = nll_breakdown,
    report = lik_report),
  
  g3l_abundancedistribution(
    'si_aut_2a',
    (aut.SI[[boot_repl]]$len52 %>% filter(year > 1999, year != 2011)) %>% filter(year < tyr - peel),
    fleets = list(),
    stocks = stocks,
    g3l_distribution_surveyindices_log(beta = if (slope_early_SIs) NULL else 1), 
    nll_breakdown = nll_breakdown,
    report = lik_report),
  
  g3l_abundancedistribution(
    'si_aut_2b',
    (aut.SI[[boot_repl]]$len60 %>% filter(year > 1999, year != 2011)) %>% filter(year < tyr - peel),
    fleets = list(),
    stocks = stocks,
    g3l_distribution_surveyindices_log(beta = 1),
    nll_breakdown = nll_breakdown,
    report = lik_report),
  
  g3l_abundancedistribution(
    'si_aut_3a',
    (aut.SI[[boot_repl]]$len72 %>% filter(year > 1999, year != 2011)) %>% filter(year < tyr - peel),
    fleets = list(),
    stocks = stocks,
    g3l_distribution_surveyindices_log(beta = 1),
    nll_breakdown = nll_breakdown,
    report = lik_report),
  
  g3l_abundancedistribution(
    'si_aut_3b',
    (aut.SI[[boot_repl]]$len80 %>% filter(year > 1999, year != 2011)) %>% filter(year < tyr - peel),
    fleets = list(),
    stocks = stocks,
    g3l_distribution_surveyindices_log(beta = 1),
    nll_breakdown = nll_breakdown,
    report = lik_report),
  
  g3l_abundancedistribution(
    'si_aut_3c',
    (aut.SI[[boot_repl]]$len92 %>% filter(year > 1999, year != 2011)) %>% filter(year < tyr - peel),
    fleets = list(),
    stocks = stocks,
    g3l_distribution_surveyindices_log(beta = 1),
    nll_breakdown = nll_breakdown,
    report = lik_report),
  
  g3l_abundancedistribution(
    'si_aut_3d',
    (aut.SI[[boot_repl]]$len100 %>% filter(year > 1999, year != 2011)) %>% filter(year < tyr - peel),
    fleets = list(),
    stocks = stocks,
    g3l_distribution_surveyindices_log(beta = 1),
    nll_breakdown = nll_breakdown,
    report = lik_report),
  
  list()
  
)
