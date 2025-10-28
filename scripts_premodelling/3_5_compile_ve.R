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
  left_join(offset_zvl$agp) %>% 
  rename(VE0 = VE) %>% 
  mutate(
    oddVE = log(VE0 / (1 - VE0)) + or70spline,
    VE = 1 / (1 + exp(-oddVE)),
  ) %>% 
  select(Key, Age, Yr, Vaccine, VE0, VE, IC)

save(pars_ve_zvl, file = here::here("pars", "pars_ve_zvl_rwa.rdata"))


## RZV: zi-exponential version selected ----- 

for(tag in c("y10_zig", "y10_zie", "y11_zig", "y11_zie", "y11m_zig", "y11m_zie")) {
  tag <- glue::as_glue(tag)
  ve_rzv <- local({load(here::here("pars", "fitted_ve_rzv_" + tag + ".rdata")); sel})
  
  
  ## Trial-based VE ----
  pars_ve_rzv <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
    left_join(ve_rzv) %>% 
    mutate(
      p_im = p0,
      VE_t = (1 - pgamma(Yr, alpha, beta)),
      Variant = "TR",
      IC = F,
      RZV_2d = p_im * VE_t,
      RZV_1d = apply_lor(p_im, offset_rzv$single) * VE_t,
      ReRZV_2d = apply_lor(p_im, offset_rzv$re) * VE_t,
      ReRZV_1d = apply_lor(p_im, offset_rzv$re + offset_rzv$single) * VE_t
    ) %>% 
    select(Key, Yr, p_im, VE_t, Variant, IC, starts_with(c("RZV", "ReRZV")))

  save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_tr_" + tag + ".rdata"))
  

  
  ## Realworld VE ----
  pars_ve_rzv <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
    left_join(ve_rzv) %>% 
    mutate(
      p_im = p0,
      lor = find_lor(mean(p_im), offset_rzv$ve0),
      p_im = apply_lor(p_im, lor),
      VE_t = (1 - pgamma(Yr, alpha, beta)),
      Variant = "RW",
      IC = F,
      RZV_2d = p_im * VE_t,
      RZV_1d = apply_lor(p_im, offset_rzv$single) * VE_t,
      ReRZV_2d = apply_lor(p_im, offset_rzv$re) * VE_t,
      ReRZV_1d = apply_lor(p_im, offset_rzv$re + offset_rzv$single) * VE_t
    )
  
  save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_rw_" + tag + ".rdata"))
}


## Select meta
tag <- glue::as_glue("y11_zie")

for(variant in c("tr", "rw")) {
  variant <- glue::as_glue(variant)
  
  file.copy(from = here::here("pars", "pars_ve_rzv_" + variant + "_" + tag + ".rdata"), 
            to = here::here("pars", "pars_ve_rzv_" + variant + ".rdata"), overwrite = T)
}

