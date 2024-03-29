## -----------------------------------------------------------------------------
##
## Catch distributions
##
## -----------------------------------------------------------------------------

minage <- lapply(stocks, gadget3::g3_stock_def, 'minage') |> unlist() |> min()
maxage <- lapply(stocks, gadget3::g3_stock_def, 'maxage') |> unlist() |> max()
minlen <- lapply(stocks, gadget3::g3_stock_def, 'minlen') |> unlist() |> min()
maxlen <- 
  lapply(stocks, gadget3::g3_stock_def, 'maxlen') |> 
  unlist() |> (\(.) replace(., is.infinite(.), NA))() |> max(na.rm = TRUE)
dl <- lapply(stocks, function(x) diff(gadget3::g3_stock_def(x, 'minlen'))) |> unlist() |> min()
  
## Query length data to create catch distribution components
## AUTUMN SURVEY
ldist.aut <-
  mfdb_sample_count(mdb, 
                    c('age', 'length'), 
                    c(list(
                      data_source = 'iceland-ldist',
                      sampling_type = 'AUT',
                      age = mfdb_interval("all",c(minage,maxage),
                                       open_ended = c("upper","lower")),
                      length = mfdb_interval("len", 
                                             seq(minlen, maxlen, by = dl),
                                             open_ended = c("upper","lower"))),
                      defaults))

## Maturity
matp.aut <- 
  mfdb_sample_count(mdb, c('maturity_stage','age','length'),
                    append(defaults,
                           list(sampling_type='AUT',
                                #age=mfdb_group(mat_ages=minage:maxage),
                                length = mfdb_interval('len',
                                                       seq(minlen, maxlen, by = 2*dl),
                                                       open_ended = c('lower','upper')),              
                                maturity_stage = mfdb_group(blingimm = 1, blingmat = 2:5))))

for(i in seq_along(matp.aut)){
  attributes(matp.aut[[i]])$age$all <- minage:maxage
}

ldist.aut[[1]]$step <- 4
matp.aut[[1]]$step <- 4

## LONGLINE
ldist.lln <- 
  mfdb_sample_count(mdb, 
                    c('age', 'length'), 
                    c(list(
                      sampling_type = 'SEA',
                      data_source = 'iceland-ldist',
                      gear = c('LLN','HLN'),
                      age = mfdb_interval("all",c(minage,maxage),
                                          open_ended = c("upper","lower")),
                      length = mfdb_interval("len", 
                                             seq(minlen, maxlen, by = dl),
                                             open_ended = c("upper","lower"))),
                      defaults))

## BOTTOM TRAWL LENGTH AND AGE
ldist.bmt <- 
  mfdb_sample_count(mdb,
                    c('age', 'length'), 
                    c(list(
                      sampling_type = 'SEA',
                      data_source = 'iceland-ldist',
                      gear=c('BMT','NPT','DSE','PSE','PGT','SHT','GIL'),
                      age = mfdb_interval("all",c(minage,maxage),
                                          open_ended = c("upper","lower")),
                      length = mfdb_interval("len", 
                                             seq(minlen, maxlen, by = dl),
                                             open_ended = c("upper","lower"))),
                      defaults))

ldist.comm <- 
  mfdb_sample_count(mdb,
                    c('age', 'length'), 
                    c(list(
                      sampling_type = 'SEA',
                      data_source = 'iceland-ldist',
                      gear=c('BMT','NPT','DSE','PSE','PGT','SHT','GIL','HLN','LLN'),
                      age = mfdb_interval("all",c(minage,maxage),
                                          open_ended = c("upper","lower")),
                      length = mfdb_interval("len", 
                                             seq(minlen, maxlen, by = dl),
                                             open_ended = c("upper","lower"))),
                      defaults))

aldist.bmt <-
  mfdb_sample_count(mdb,
                    c('age', 'length'),
                    c(list(sampling_type = 'SEA',
                           data_source = 'iceland-aldist',
                           gear=c('BMT','NPT','DSE','PSE','PGT','SHT','GIL'),
                           age = mfdb_interval('age',c(minage:maxage),open_ended = c('upper')),
                           length = mfdb_interval("len",
                                                  seq(minlen, maxlen, by = dl),
                                                  open_ended = c("upper","lower"))),
                      defaults))

## OUTPUT
if (run_bootstrap){
  save(ldist.aut, matp.aut,
       ldist.lln, ldist.bmt, aldist.bmt, ldist.comm,
       file = file.path(base_dir, 'data', 'catchdistribution_bootstrap_data.Rdata'))
}else{
  save(ldist.aut, matp.aut,
       ldist.lln, ldist.bmt, aldist.bmt, ldist.comm,
       file = file.path(base_dir, 'data', 'catchdistribution_data.Rdata'))
}

