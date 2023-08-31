## -----------------------------------------------------------------------------
##
## Catches by fleet
##
## -----------------------------------------------------------------------------

## Survey
aut_landings <- 
  structure(data.frame(year = defaults$year, step = 4, area = 1, total_weight = 1),
            area_group = mfdb_group(`1` = 1))

## Long- & Handline
lln_landings <- mfdb_sample_totalweight(mdb, NULL, c(list(
  gear=c('HLN','LLN'),
  sampling_type = 'LND',
  species = defaults$species),
  defaults))
  
## Bottom trawl etc (Gilnet included here)
bmt_landings <- mfdb_sample_totalweight(mdb, NULL, c(list(
  gear=c('BMT','NPT','DSE','PSE','PGT','SHT','GIL'),
  sampling_type = 'LND',
  species = defaults$species),
  defaults))

foreign_landings <-
  mfdb_sample_totalweight(mdb, NULL,
                          c(list(
                            sampling_type = 'FLND',
                            data_source = c('lods.foreign.landings','statlant.foreign.landings'),
                            species = defaults$species),
                            defaults))

## For single fleet option:
comm_landings <- mfdb_sample_totalweight(mdb, NULL, c(list(
  gear=c('BMT','NPT','DSE','PSE','PGT','SHT','GIL','HLN', 'LLN'),
  sampling_type = 'LND',
  species = defaults$species),
  defaults)) 
  

if (TRUE){
  save(lln_landings, bmt_landings, foreign_landings, aut_landings, comm_landings,
       file = file.path(base_dir, "data", "fleet_data.Rdata"))
}
