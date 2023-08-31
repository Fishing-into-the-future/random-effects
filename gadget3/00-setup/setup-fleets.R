## -----------------------------------------------------------------------------
##
## Configure fleets and associated actions
##
## -----------------------------------------------------------------------------

## FLEET STOCKS
## Longline
lln <- 
  g3_fleet('lln') %>% 
  g3s_livesonareas(areas[c('1')])

## Bottom trawl
bmt <- 
  g3_fleet('bmt') %>% 
  g3s_livesonareas(areas[c('1')])

## Foreign
foreign <- 
  g3_fleet('foreign') %>% 
  g3s_livesonareas(areas[c('1')])

## Autumn survey
aut <- 
  g3_fleet("aut") %>%
  g3s_livesonareas(areas[c("1")])

## For single fleet option
comm <- 
  g3_fleet("comm") %>% 
  g3s_livesonareas(areas[c("1")])

## ------------------------------------------------------------------------------

## Fleet actions
fleet_actions <- 
  list(
    if (single_fleet){
      comm %>% 
        g3a_predate_fleet(stocks,
                          suitabilities = 
                            stocks %>% 
                            set_names(.,map(.,'name')) %>% 
                            map(function(x){
                              if (dome_comm){
                                g3_suitability_andersen(
                                  p0 = g3_parameterized('andersen.p0', by_stock = 'species'),
                                  p1 = g3_parameterized('comm.p1', by_stock = 'species'),
                                  p2 = g3_parameterized('andersen.p2', by_stock = 'species'),
                                  p3 = g3_parameterized('comm.p3', by_stock = 'species'),
                                  p4 = g3_parameterized('comm.p4', by_stock = 'species'),
                                  p5 = g3_parameterized('andersen.L', by_stock = 'species')
                                )
                              }else{
                                g3_suitability_exponentiall50(g3_parameterized('comm.alpha', by_stock = 'species'),
                                                              g3_parameterized('comm.l50', by_stock = 'species')) 
                              }
                            }),
                          catchability_f = g3a_predate_catchability_totalfleet(g3_timeareadata('comm_landings', comm_landings[[1]] %>%
                                                                                                 mutate(area = as.numeric(area),
                                                                                                        step = as.numeric(step),
                                                                                                        year = as.numeric(year))))) 
    }else{
      c(
        lln %>% 
          g3a_predate_fleet(stocks,
                            suitabilities = 
                              stocks %>% 
                              set_names(.,map(.,'name')) %>% 
                              map(function(x){
                                if (dome_lln){
                                  g3_suitability_andersen(
                                    p0 = g3_parameterized('andersen.p0', by_stock = 'species'),
                                    p1 = g3_parameterized('lln.p1', by_stock = 'species'),
                                    p2 = g3_parameterized('andersen.p2', by_stock = 'species'),
                                    p3 = g3_parameterized('lln.p3', by_stock = 'species'),
                                    p4 = g3_parameterized('lln.p4', by_stock = 'species'),
                                    p5 = g3_parameterized('andersen.L', by_stock = 'species')
                                  )
                                }else{
                                  g3_suitability_exponentiall50(g3_parameterized('lln.alpha', by_stock = 'species'),
                                                                g3_parameterized('lln.l50', by_stock = 'species')) 
                                }
                              }),
                            catchability_f = g3a_predate_catchability_totalfleet(g3_timeareadata('lln_landings', lln_landings[[1]] %>%
                                                                                                   mutate(area = as.numeric(area),
                                                                                                          step = as.numeric(step),
                                                                                                          year = as.numeric(year))))),
        bmt %>% 
          g3a_predate_fleet(stocks,
                            suitabilities = 
                              stocks %>% 
                              set_names(.,map(.,'name')) %>% 
                              map(function(x){
                                if (dome_lln){
                                  g3_suitability_andersen(
                                    p0 = g3_parameterized('andersen.p0', by_stock = 'species'),
                                    p1 = g3_parameterized('bmt.p1', by_stock = 'species'),
                                    p2 = g3_parameterized('andersen.p2', by_stock = 'species'),
                                    p3 = g3_parameterized('bmt.p3', by_stock = 'species'),
                                    p4 = g3_parameterized('bmt.p4', by_stock = 'species'),
                                    p5 = g3_parameterized('andersen.L', by_stock = 'species')
                                  )
                                }else{
                                  g3_suitability_exponentiall50(g3_parameterized('bmt.alpha', by_stock = 'species'),
                                                                g3_parameterized('bmt.l50', by_stock = 'species')) 
                                }
                              }),
                            catchability_f = g3a_predate_catchability_totalfleet(g3_timeareadata('bmt_landings', bmt_landings[[1]] %>%
                                                                                                   mutate(area = as.numeric(area),
                                                                                                          step = as.numeric(step),
                                                                                                          year = as.numeric(year)))))
        
      )
    },
    
    foreign  %>% 
      g3a_predate_fleet(stocks,
                        suitabilities = 
                          stocks %>% 
                          set_names(.,map(.,'name')) %>% 
                          map(function(x){
                            if (dome_lln){
                              g3_suitability_andersen(
                                p0 = g3_parameterized('andersen.p0', by_stock = 'species'),
                                p1 = g3_parameterized(ifelse(single_fleet, 'comm.p1', 'lln.p1'), by_stock = 'species'),
                                p2 = g3_parameterized('andersen.p2', by_stock = 'species'),
                                p3 = g3_parameterized(ifelse(single_fleet, 'comm.p3', 'lln.p3'), by_stock = 'species'),
                                p4 = g3_parameterized(ifelse(single_fleet, 'comm.p4', 'lln.p4'), by_stock = 'species'),
                                p5 = g3_parameterized('andersen.L', by_stock = 'species')
                              )
                            }else{
                              g3_suitability_exponentiall50(g3_parameterized(ifelse(single_fleet, 'comm.alpha', 'lln.alpha'), by_stock = 'species'),
                                                            g3_parameterized(ifelse(single_fleet, 'comm.l50', 'lln.l50'), by_stock = 'species')) 
                            }
                          }),
                        catchability_f = g3a_predate_catchability_totalfleet(g3_timeareadata('foreign_landings', foreign_landings[[1]] %>%
                                                                                               mutate(area = as.numeric(area),
                                                                                                      step = as.numeric(step),
                                                                                                      year = as.numeric(year))))),
    aut %>% 
      g3a_predate_fleet(stocks,
                        suitabilities = 
                          stocks %>% 
                          set_names(.,map(.,'name')) %>% 
                          map(function(x) g3_suitability_exponentiall50(g3_parameterized('aut.alpha', by_stock = 'species'),
                                                                        g3_parameterized('aut.l50', by_stock = 'species'))),
                        catchability_f = g3a_predate_catchability_totalfleet(g3_timeareadata('aut_landings', aut_landings %>%
                                                                                               mutate(area = as.numeric(area),
                                                                                                      step = as.numeric(step),
                                                                                                      year = as.numeric(year))))),
    list()
  )
     
  
    
    
    