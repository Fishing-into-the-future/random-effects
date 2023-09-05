## Useful constants
## Weight-Length relationship (Autumn survey)
lw.constants <- 
  mfdb_dplyr_sample(mdb) %>% 
  filter(species == local(defaults$species),
         sampling_type == 'AUT',
         !is.na(weight),
         !is.na(length),
         weight > 0,
         length > 10) %>% 
  select(length,weight) %>% 
  collect(n=Inf) %>% 
  lm(log(weight/1e3)~log(length),.) %>% 
  broom::tidy() %>% 
  select(estimate)
## transport back to right dimension
lw.constants$estimate[1] <- exp(lw.constants$estimate[1])

## initial conditions sigma
init.sigma <- 
  mfdb_dplyr_sample(mdb) %>% 
  dplyr::filter(species == local(defaults$species), 
                !is.na(length),
                !is.na(age),
                year > 1950,
                length > 10)  %>% 
  dplyr::select(age,length) %>% 
  dplyr::collect(n=Inf) %>% 
  dplyr::group_by(age) %>% 
  dplyr::summarise(ml=mean(length,na.rm=TRUE),ms=sd(length,na.rm=TRUE), n=length(na.exclude(length)))

## Initial coefficients for sd
init.sigma.coef <- 
  init.sigma %>% 
  filter(age >= min(unlist(lapply(stocks, gadget3::g3_stock_def, 'minage'))) & 
         age <= max(unlist(lapply(stocks, gadget3::g3_stock_def, 'minage')))) %>% 
  lm(I(ms/ml)~I(1/age) + age, data = .) %>% 
  coefficients() %>% 
  setNames(c('alpha', 'beta', 'gamma'))


## initial guess for the maturity ogive:
mat.l50 <- 
  mfdb_dplyr_sample(mdb) %>% 
  filter(species == local(defaults$species),
         sampling_type == 'AUT',
         !is.na(maturity_stage)) %>% 
  select(length,maturity_stage) %>% 
  group_by(length,maturity_stage) %>% 
  dplyr::summarise(n=n()) %>% 
  group_by(length) %>% 
  dplyr::mutate(p=n/sum(n)) %>% 
  ungroup() %>% 
  filter(maturity_stage=='2',p>0.50,length > 20) %>% 
  dplyr::summarise(l50=min(length)) %>% 
  collect(n=Inf)

## No data so taken from the Faroese stock annex 
## Magnusson and Magnusson 95 estimated and L and A50 for Icelandic males (9) and females (11)
mat.a50 <- 10

## OUTPUT
save(lw.constants, mat.l50, init.sigma, mat.a50, init.sigma.coef,
     file = file.path(base_dir, 'data', 'init_param.Rdata'))

