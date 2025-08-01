library(tidyverse)
library(tidybayes)


source(here::here("scripts_premodelling", "fn_ors.R"))


n_sims <- 2e3
load(here::here("pars", "fitted_ve_offset.rdata"))


## ZVL: zi-gamma version selected ----- 
ve_zvl <- local({load(here::here("pars", "fitted_ve_zvl_zig.rdata")); sel})


pars_ve_zvl <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_zvl) %>% 
  mutate(
    VE = p0 * (1 - pgamma(Yr, alpha, beta)),
    VE = align_to(VE, mean(VE[Yr == 1]), tar = offset_zvl$ve0),
    Vaccine = "ZVL",
    IC = F
  ) %>% 
  select(Key, Yr, Vaccine, VE, IC)

save(pars_ve_zvl, file = here::here("pars", "pars_ve_zvl_rw.rdata"))


pars_ve_zvl <- pars_ve_zvl %>% 
  crossing(Age = 50:100) %>% 
  filter(Age - Yr >= 50) %>% 
  left_join(ve_zvl) %>% 
  left_join(offset_zvl$agp) %>% 
  rename(VE0 = VE) %>% 
  mutate(
    oddVE = log(VE0 / (1 - VE0)) + or70spline,
    VE = 1 / (1 + exp(-oddVE)),
  ) %>% 
  select(Key, Age, Yr, Vaccine, VE0, VE, IC)

save(pars_ve_zvl, file = here::here("pars", "pars_ve_zvl_rwa.rdata"))


## RZV: zi-exponential version selected ----- 

tag <- glue::as_glue("y11m_zie")
ve_rzv <- local({load(here::here("pars", "fitted_ve_rzv_" + tag + ".rdata")); sel})


pars_ve_tr <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_rzv) %>% 
  mutate(
    VE = p0 * (1 - pgamma(Yr, alpha, beta)),
    Vaccine = "RZV",
    Variant = "TR",
    IC = F
  ) %>% 
  select(Key, Vaccine, Variant, Yr, VE, IC)


pars_ve_rzv <- pars_ve_tr
save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_tr.rdata"))
save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_tr_zie.rdata"))


pars_ve_rw <- pars_ve_tr %>% 
  mutate(
    VE = align_to(VE, mean(VE[Yr == 1]), tar = offset_rzv$ve0),
    Variant = "RW"
  )


pars_ve_rzv <- pars_ve_rw
save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_rw.rdata"))
save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_rw_zie.rdata"))


pars_ve_rzv <- pars_ve_rw %>%
  mutate(
    VE = apply_lor(VE, offset_rzv$re),
    Variant = "ReRZV"
  )

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rerzv_rw.rdata"))


pars_ve_rzv <- pars_ve_rw %>%
  mutate(
    VE = apply_lor(VE, offset_rzv$single),
    Variant = "RZV_Single"
  )

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_single_rw.rdata"))


pars_ve_rzv <- pars_ve_rw %>%
  mutate(
    VE = apply_lor(VE, offset_rzv$re + offset_rzv$single),
    Variant = "ReRZV_Single"
  )

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_resingle_rw.rdata"))


pars_ve_rzv <- pars_ve_rw %>% 
  group_by(Key) %>% 
  mutate(
    VE = case_when(
      Yr <= 11 ~ VE,
      T ~ VE[Yr == 11]
    )
  )

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_long_rw_zie.rdata"))


pars_ve_rzv <- pars_ve_rw %>% 
  group_by(Key) %>% 
  mutate(
    VE = case_when(
      Yr <= 11 ~ VE,
      T ~ 0
    )
  )

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_short_rw_zie.rdata"))


## RZV gamma -----

tag <- glue::as_glue("y11m_zig")
ve_rzv <- local({load(here::here("pars", "fitted_ve_rzv_" + tag + ".rdata")); sel})


pars_ve_tr <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_rzv) %>% 
  mutate(
    VE = p0 * (1 - pgamma(Yr, alpha, beta)),
    Vaccine = "RZV",
    Variant = "TR",
    IC = F
  ) %>% 
  select(Key, Vaccine, Variant, Yr, VE, IC)

pars_ve_rzv <- pars_ve_tr 
save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_tr_zig.rdata"))


pars_ve_rw <- pars_ve_tr %>% 
  mutate(
    VE = align_to(VE, mean(VE[Yr == 1]), tar = offset_rzv$ve0),
    Variant = "RW",
  )

pars_ve_rzv <- pars_ve_tr 
save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_rw_zig.rdata"))


