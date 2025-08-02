
compile_zvl <- function(n_sims) {
  pars_ve_zvl <- local({
    load(here::here("pars", "pars_ve_zvl_rw.rdata"))
    ve_zvl_age <- read_csv(here::here("pars", "pars_ve_zvl_agp.csv"))

    crossing(Age = 50:100, Key = 1:max(pars_ve_zvl$Key), Yr = 1:50) %>% 
      filter(Age - Yr >= 50) %>% 
      left_join(pars_ve_zvl, by = c("Key", "Yr")) %>% 
      left_join(ve_zvl_age, by = "Age") %>% 
      rename(VE0 = VE) %>% 
      mutate(
        oddVE = log(VE0 / (1 - VE0)) + or70spline,
        VE = 1 / (1 + exp(-oddVE)),
        Vaccine = "ZVL"
      ) %>% 
      select(Key, Age, Yr, Vaccine, VE, IC) %>% 
      sample_table(n_sims = n_sims)
    
  }) 
  
  f <- here::here("pars", "pars_ve_zvl_rwa_zlg.rdata")
  save(pars_ve_zvl, file = f)
  return(f) 
}

