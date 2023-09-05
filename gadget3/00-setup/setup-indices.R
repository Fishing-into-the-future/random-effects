## -----------------------------------------------------------------------------
## Survey indices
## -----------------------------------------------------------------------------

## Autumn survey indices
aut.SI <- 
  mfdb_sample_count(mdb, 
                    c('length'), 
                    c(list(
                      data_source = 'iceland-ldist',
                      sampling_type = 'AUT',
                      length = mfdb_interval('len', c(20,52,60,72,80,92,100,140), open_ended = c('lower', 'upper'))),
                      defaults)) %>% 
  #.[[1]] %>%
  purrr::map(function(y){
    y %>% 
      split(.,~length) %>% 
      purrr::map(function(x){
        structure(x, length = attr(x, 'length')[unique(x$length)])
      })
  }) 


## OUTPUT
if (run_bootstrap){
  save(aut.SI, file = file.path(base_dir, 'data', 'bootstrap_indices.Rdata'))
}else{
  save(aut.SI, file = file.path(base_dir, 'data', 'indices.Rdata'))  
}


