library(tidyverse)
library(readxl)

theme_set(theme_bw())


source(here::here("scripts_premodelling", "fn_ors.R"))

## Data loading
dat <- read_xlsx(here::here("data", "processed_vaccine", "VE.xlsx"))


dat %>% 
  filter(Source == "Izurieta 2021")


dat %>% 
  filter(Source == "Sun")


offset_rzv <- list(
  ve0 = 0.701,
  single = find_lor(0.701, 0.569), # Izurieta HS 2021,
  re = find_lor(0.79, 0.75), # Mbinta JF 2022,
  ic = find_lor(0.705, 0.641),
  ic_single = find_lor(0.705, 0.37)
)

### ZVL
dat %>% 
  filter(Source == "Walker 2018")


dat_ve_agp <- dat %>% 
  filter(Use) %>% 
  filter(Type == "HZ") %>% 
  filter(Vaccine == "Zostavax") %>% 
  filter(Source == "Mbinta 2022") %>% 
  filter(!is.na(`Age-group`)) %>% 
  mutate(
    Agp = c(55, 65, 75, 85),
    M = M / 100,
    odd = log(M / (1 - M)),
    or = (odd - odd[2])
  ) %>% 
  select(Agp, or, M)


d_agp <- tibble(Age = 50:100) %>% 
  mutate(
    Agp = case_when(
      Age < 60 ~ 55,
      Age < 70 ~ 65, 
      Age < 80 ~ 75,
      T ~ 85
    )
  ) %>% 
  left_join(dat_ve_agp) %>% 
  mutate(
    odd = log(M / (1 - M)),
    or = odd - odd[Age == 70],
    or_spline = splinefun(unique(Agp), unique(or), method = "natural")(Age)
  )


g_test_agp <- d_agp %>% 
  mutate(
    M2 = splinefun(unique(Agp), unique(M), method = "natural")(Age),
    odd3 = splinefun(unique(Agp), unique(odd), method = "natural")(Age),
    M3 = 1 / (1 + exp(-odd3)),
    or4 = splinefun(unique(Agp), unique(or), method = "natural")(Age),
    odd4 = or4 + odd[Age == 70],
    M4 = 1 / (1 + exp(-odd4))
  ) %>%   
  ggplot() +
  geom_line(aes(x = Age, y = M)) +
  geom_line(aes(x = Age, y = M2)) +
  geom_line(aes(x = Age, y = M3)) +
  geom_line(aes(x = Age, y = M4)) +
  expand_limits(y = 0:1)

g_test_agp



offset_zvl <- list(
  ve0 = 0.69,
  agp = d_agp %>% 
    select(Age, or70 = or, or70spline = or_spline)
)  


save(offset_rzv, offset_zvl, file = here::here("pars", "fitted_ve_offset.rdata"))




