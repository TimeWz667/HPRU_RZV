library(tidyverse)
library(tidybayes)



n_sims <- 2e3

## ZVL -----
ve_zvl <- local({load(here::here("pars", "fitted_ve_zvl_zlg.rdata")); sel})
ve_zvl_agp <- read_csv(here::here("pars", "fitted_ve_zvl_agp.csv"))


pars_ve_zvl <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_zvl) %>% 
  mutate(
    VE = p0 * (1 - pgamma(Yr, alpha, beta)),
    Vaccine = "ZVL",
    IC = F
  ) %>% 
  select(Key, Yr, Vaccine, VE, IC)

save(pars_ve_zvl, file = here::here("pars", "pars_ve_zvl_rw.rdata"))


pars_ve_zvl <- crossing(Age = 50:100, Key = 1:n_sims, Yr = 1:50) %>% 
  filter(Age - Yr >= 50) %>% 
  left_join(ve_zvl) %>% 
  left_join(ve_zvl_agp) %>% 
  mutate(
    VE0 = p0 * (1 - pgamma(Yr, alpha, beta)),
    oddVE = log(VE0 / (1 - VE0)) + or70spline,
    VE = 1 / (1 + exp(-oddVE)),
    Vaccine = "ZVL",
    IC = F
  ) %>% 
  select(Key, Age, Yr, Vaccine, VE0, VE, IC)

save(pars_ve_zvl, file = here::here("pars", "pars_ve_zvl_rwa.rdata"))

pars_ve_zvl %>% 
  group_by(Age, Yr) %>% 
  summarise(across(c(VE, VE0), mean)) %>% 
  ggplot() +
  geom_point(aes(x = Age, y = VE, colour = Yr)) +
  scale_y_continuous("Vaccine effectiveness, %", labels = scales::percent) +
  expand_limits(y = 0:1)


ve_zvl <- local({load(here::here("pars", "pars_ve_zvl_zle.rdata")); sel})
ve_zvl_agp <- read_csv(here::here("pars", "pars_ve_zvl_agp.csv"))


pars_ve_zvl <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_zvl) %>% 
  mutate(
    VE = p0 * (1 - pgamma(Yr, alpha, beta)),
    Vaccine = "ZVL",
    IC = F
  ) %>% 
  select(Key, Yr, Vaccine, VE, IC)

save(pars_ve_zvl, file = here::here("pars", "pars_ve_zvl_rw_zle.rdata"))



pars_ve_zvl <- crossing(Age = 50:100, Key = 1:n_sims, Yr = 1:50) %>% 
  filter(Age - Yr >= 50) %>% 
  left_join(ve_zvl) %>% 
  left_join(ve_zvl_agp) %>% 
  mutate(
    VE0 = p0 * (1 - pgamma(Yr, alpha, beta)),
    oddVE = log(VE0 / (1 - VE0)) + or70spline,
    VE = 1 / (1 + exp(-oddVE)),
    Vaccine = "ZVL",
    IC = F
  ) %>% 
  select(Key, Age, Yr, Vaccine, VE0, VE, IC)

save(pars_ve_zvl, file = here::here("pars", "pars_ve_zvl_rwa_zle.rdata"))



## RZV gamma ----- 

find_lor <- function(p0, p1) log(p1 / (1 - p1)) - log(p0 / (1 - p0))
apply_lor <- function(p0, lor) 1 / (1 + exp(-log(p0 / (1 - p0)) - lor))


tag <- glue::as_glue("y11m")
ve_rzv <- local({load(here::here("pars", "pars_ve_rzv_" + tag + "_zl_gamma.rdata")); sel})

ve_y1 <- ve_rzv %>% 
  mutate(
    VE = p0 * (1 - pgamma(1, alpha, beta))
  ) %>% 
  pull(VE) %>% mean


lor_rw <- find_lor(ve_y1, 0.79)
lor_single <- find_lor(0.701, 0.569) # Lzurieta HS 2021

lor_re <- find_lor(0.79, 0.75)

lor_ic <- find_lor(0.705, 0.641)
lor_ic_single <- find_lor(0.705, 0.37)


save(lor_rw, lor_re, lor_single, lor_ic, lor_ic_single, file = here::here("pars", "pars_ve_lor.rdata"))



### RZV gamma -----
pars_ve_rzv <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_rzv) %>% 
  mutate(
    VE = p0 * (1 - pgamma(Yr, alpha, beta)),
    VE = apply_lor(VE, lor_rw),
    Vaccine = "RZV",
    Variant = "RW",
    IC = F
  ) %>% 
  select(Key, Vaccine, Yr, VE, IC)

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_rw_zlg.rdata"))


pars_ve_rzv_booster <- pars_ve_rzv %>% 
  mutate(Variant = "RW_reRZV")

save(pars_ve_rzv_booster, file = here::here("pars", "pars_ve_rerzv_rw_zlg.rdata"))


pars_ve_rzv_booster <- pars_ve_rzv %>% 
  mutate(
    VE = apply_lor(VE, lor_single),
    Variant = "RW_reRZV(S)"
  )

save(pars_ve_rzv_booster, file = here::here("pars", "pars_ve_rerzv_single_rw_zlg.rdata"))


pars_ve_rzv <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_rzv) %>% 
  mutate(
    VE = p0 * (1 - pgamma(pmin(Yr, 10), alpha, beta)) ,
    VE = apply_lor(VE, lor_rw),
    Vaccine = "RZV",
    Variant = "RW",
    IC = F
  ) %>% 
  select(Key, Vaccine, Yr, VE, IC)

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_rw_long_zlg.rdata"))


pars_ve_rzv <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_rzv) %>% 
  mutate(
    VE = p0 * (1 - pgamma(pmin(Yr, 10), alpha, beta)) ,
    VE = apply_lor(VE, lor_rw),
    VE = ifelse(Yr <= 10, VE, 0),
    Vaccine = "RZV",
    Variant = "RW",
    IC = F
  ) %>% 
  select(Key, Vaccine, Yr, VE, IC)

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_rw_short_zlg.rdata"))


## RZV exp -----

ve_rzv <- local({load(here::here("pars", "pars_ve_rzv_" + tag + "_zl_exp.rdata")); sel})


pars_ve_rzv <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_rzv) %>% 
  mutate(
    VE = p0 * (1 - pgamma(Yr, alpha, beta)) * 0.91,
    Vaccine = "RZV",
    Variant = "RW",
    IC = F
  ) %>% 
  select(Key, Vaccine, Yr, VE, IC)

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_rw_zle.rdata"))


pars_ve_rzv_booster <- pars_ve_rzv %>% 
  mutate(Variant = "RW_reRZV")

save(pars_ve_rzv_booster, file = here::here("pars", "pars_ve_rerzv_rw_zle.rdata"))


pars_ve_rzv_booster <- pars_ve_rzv %>% 
  mutate(
    VE = apply_lor(VE, lor_single),
    Variant = "RW_reRZV(S)"
  )

save(pars_ve_rzv_booster, file = here::here("pars", "pars_ve_rerzv_single_rw_zle.rdata"))


bind_rows(
  pars_ve_rzv %>% mutate(Variant = "RW"),
  pars_ve_rzv_booster
) %>% 
  group_by(Yr, Variant) %>% 
  summarise(
    M = mean(VE),
    L = quantile(VE, 0.025),
    U = quantile(VE, 0.975)
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Yr, ymin = L, ymax = U, fill = Variant), alpha = 0.3) +
  geom_line(aes(x = Yr, y = M, colour = Variant)) + 
  scale_y_continuous("Vaccine effectiveness, %", labels = scales::percent) +
  expand_limits(y = 0:1)


pars_ve_rzv <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_rzv) %>% 
  mutate(
    VE = p0 * (1 - pgamma(pmin(Yr, 10), alpha, beta)) ,
    VE = apply_lor(VE, lor_rw),
    Vaccine = "RZV",
    Variant = "RW",
    IC = F
  ) %>% 
  select(Key, Vaccine, Yr, VE, IC)

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_rw_long_zle.rdata"))


pars_ve_rzv <- crossing(Key = 1:n_sims, Yr = 1:50) %>% 
  left_join(ve_rzv) %>% 
  mutate(
    VE = p0 * (1 - pgamma(pmin(Yr, 10), alpha, beta)) ,
    VE = apply_lor(VE, lor_rw),
    VE = ifelse(Yr <= 10, VE, 0),
    Vaccine = "RZV",
    Variant = "RW",
    IC = F
  ) %>% 
  select(Key, Vaccine, Yr, VE, IC)

save(pars_ve_rzv, file = here::here("pars", "pars_ve_rzv_rw_short_zle.rdata"))


