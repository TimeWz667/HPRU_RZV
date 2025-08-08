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
apply_ve_offset <- function(ve, off) {
  prop <- ve %>% 
    filter(Yr == 1) %>% 
    mutate(
      prop = apply_lor(VE, off) / VE 
    ) %>% 
    select(Key, prop)
  
  
  ve %>%
    left_join(prop) %>%
    mutate(
      VE = VE * prop,
    ) %>% select(-prop)
}


for(tag in c("y10_zig", "y10_zie", "y11_zig", "y11_zie", "y11m_zig", "y11m_zie")) {
  tag <- glue::as_glue(tag)
  ve_rzv <- local({load(here::here("pars", "fitted_ve_rzv_" + tag + ".rdata")); sel})
  
  
  ## Trial-based VE ----
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
  save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_uv2_tr_" + tag + ".rdata"))
  
  
  pars_ve_rzv <- pars_ve_tr %>%
    apply_ve_offset(offset_rzv$single) %>% 
    mutate(Variant = "RZV_Single")
  
  save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_uv1_tr_" + tag + ".rdata"))
  
  
  pars_ve_rzv <- pars_ve_tr %>%
    apply_ve_offset(offset_rzv$re) %>% 
    mutate(Variant = "ReRZV")
  
  save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_re2_tr_" + tag + ".rdata"))
  
  
  pars_ve_rzv <- pars_ve_tr %>%
    apply_ve_offset(offset_rzv$re + offset_rzv$single) %>% 
    mutate(Variant = "ReRZV_Single")
  
  save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_re1_tr_" + tag + ".rdata"))
  
  
  
  ## Realworld VE ----
  prop <- pars_ve_tr %>% 
    filter(Yr == 1) %>% 
    mutate(
      lor = find_lor(mean(VE), offset_rzv$ve0),
      prop = apply_lor(VE, lor) / VE 
    ) %>% 
    select(Key, prop)
  
  
  pars_ve_rw <- pars_ve_tr %>% 
    left_join(prop) %>% 
    mutate(
      VE = VE * prop,
      Variant = "RW"
    ) %>% select(-prop)
  
  
  pars_ve_rzv <- pars_ve_rw
  save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_uv2_rw_" + tag + ".rdata"))
  
  
  pars_ve_rzv <- pars_ve_rw %>%
    apply_ve_offset(offset_rzv$single) %>% 
    mutate(Variant = "RZV_Single")
  
  save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_uv1_rw_" + tag + ".rdata"))
  
  
  pars_ve_rzv <- pars_ve_rw %>%
    apply_ve_offset(offset_rzv$re) %>% 
    mutate(Variant = "ReRZV")
  
  save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_re2_rw_" + tag + ".rdata"))
  
  
  pars_ve_rzv <- pars_ve_rw %>%
    apply_ve_offset(offset_rzv$re + offset_rzv$single) %>% 
    mutate(Variant = "ReRZV_Single")
  
  save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_re1_rw_" + tag + ".rdata"))
  
  
  # pars_ve_rzv <- pars_ve_rw %>% 
  #   group_by(Key) %>% 
  #   mutate(
  #     VE = case_when(
  #       Yr <= 11 ~ VE,
  #       T ~ VE[Yr == 11]
  #     )
  #   )
  # 
  # save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_uv2long_rw_" + tag + ".rdata"))
  # 
  # 
  # pars_ve_rzv <- pars_ve_rw %>% 
  #   group_by(Key) %>% 
  #   mutate(
  #     VE = case_when(
  #       Yr <= 11 ~ VE,
  #       T ~ 0
  #     )
  #   )
  # 
  # save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_uv2short_rw_" + tag + ".rdata"))
}


## Select meta
tag <- glue::as_glue("y11_zie")

for(variant in c("uv2_tr", "re2_tr", "uv1_tr", "re1_tr", 
                 "uv2_rw", "re2_rw", "uv1_rw", "re1_rw")) {
  variant <- glue::as_glue(variant)
  
  file.copy(from = here::here("pars", "pars_ve_rzv_" + variant + "_" + tag + ".rdata"), 
            to = here::here("pars", "pars_ve_rzv_" + variant + ".rdata"), overwrite = T)
}

