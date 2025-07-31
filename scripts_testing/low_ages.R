library(targets)
library(tidyverse)

theme_set(theme_bw())


for (d in dir("R")) source(here::here("R", d))






pars_proj <- tar_read(pars_proj, 1)

prefix <- gsub("Year0", "", names(pars_proj)[1])

names(pars_proj) <- gsub(prefix, "", pars_proj %>% names())









strategy <- strategy_changeonly

a_proj_ce <- function(pars, strategy, year0 = 2013, year_vac = 2033, year1 = 2050, age0 = 50, age1 = 100) {
  p_uptake <- pars$Uptake
  p_demo <- pars$Demography
  
  population <- with(model_proj, {
    yr <- year0
    pop <- tibble(Year = yr, Age = age0:age1, Vaccine = "None", TimeVac = -1, Prop = 1) 
    pop <- pop %>% uptaking(yr = yr, strategy = strategy, pars_uptake = p_uptake)
    
    collector <- pop
    
    while (yr < year1) {
      yr <- yr + 1
      pop <- pop %>% 
        ageing(yr, age0, age1) %>% 
        uptaking(yr = yr, strategy = strategy, pars_uptake = p_uptake)
      collector <- collector %>% bind_rows(pop)
      # print(nrow(collector))
    }
    
    collector
  }) %>% 
    left_join(p_demo$N, by = c("Year", "Age")) %>%
    mutate(
      N = N * Prop
    ) %>% 
    select(-Prop)
  
  
  ves_rzv <- bind_rows(
    pars$VE_RZV_2d,
    pars$VE_RZV_1d,
    pars$VE_ReRZV_2d,
    pars$VE_ReRZV_1d,
  )
  
  ves_zvl <- pars$VE_ZVL
  
  
  population <- population %>% crossing(Key = 1:pars$N_Sims)
  
  vaccinated <- bind_rows(
    population %>% 
      filter(Vaccine == "ZVL") %>% 
      left_join(ves_zvl, by = c("Key", "Age", "Vaccine", "TimeVac")),
    population %>% 
      filter(Vaccine == "None") %>% 
      mutate(Protection = 0),
    population %>% 
      filter(!(Vaccine %in% c("ZVL", "None"))) %>% 
      left_join(ves_rzv, by = c("Key", "Vaccine", "TimeVac"))
  ) %>% 
    mutate(
      N_Uptake_ZVL = ifelse((TimeVac == 1) * (Vaccine == "ZVL"), N, 0),
      N_Uptake_RZV = ifelse((TimeVac == 1) * endsWith(Vaccine, "RZV_1d"), N, 0) +
        ifelse((TimeVac == 1) * endsWith(Vaccine, "RZV_2d"), 2 * N, 0),
      N_Covered_ZVL = ifelse((TimeVac == 1) * (Vaccine == "ZVL"), N, 0),
      N_Covered_RZV = ifelse((TimeVac == 1) * !(Vaccine %in% c("None", "ZVL")), N, 0),
      
    ) %>% 
    group_by(Key, Year, Age) %>% 
    summarise(
      Protection = weighted.mean(Protection, w = N),
      across(starts_with("N"), sum)
    ) %>% 
    ungroup()
  
  
  sims <- vaccinated %>% 
    left_join(p_demo$DeathIm %>% select(Year, Age, r_mor_bg = r_death), by = c("Year", "Age")) %>% 
    inner_join(pars$Epidemiology, by = c("Key", "Age")) %>% 
    mutate(
      r_mor = (1 - Protection) * r_mor_hz + r_mor_bg,
      p_mor_hz = (1 - Protection) * r_mor_hz * (1 - exp(- r_mor)) / r_mor,
      N_HZ = (1 - Protection) * r_hz * N,
      N_HZ_GP = p_gp * N_HZ,
      N_HZ_Hosp = N_HZ - N_HZ_GP,
      N_HZ_PHN = p_phn * N_HZ,
      N_HZ_PHN_GP = p_gp * N_HZ_PHN
    ) %>% 
    select(-starts_with(c("r_", "p_"))) %>% 
    left_join(pars$CostEff, by = c("Key", "Age")) %>% 
    filter(Year >= 2023) %>% 
    mutate(
      Q_Life = N,
      Q_HZ = - N_HZ * QL_ph,
      Q_HZ_Norm = - N_HZ * QL_pn,
      Q_All = Q_Life + Q_HZ,
      C_Hosp = N_HZ_Hosp * cost_hosp_pp_inf,
      C_GP_NonPHN = (N_HZ_GP - N_HZ_PHN_GP) * cost_GP_pp_non_PHN_HZ_inf,
      C_GP_PHN = N_HZ_PHN_GP * cost_GP_pp_PHN_inf,
      C_GP = C_GP_NonPHN + C_GP_PHN,
      C_Med = C_GP + C_Hosp,
      across(starts_with("C_Vac"), \(x) ifelse(is.na(x), 0, x)),
      across(starts_with(c("C_", "Q_", "N_")), \(x) ifelse(is.na(x), 0, x))
    ) %>% 
    filter(Age >= 60) %>% 
    group_by(Key, Year) %>% 
    summarise(
      Protection = sum(N * Protection) / sum(N),
      N_Vac = sum(N_Uptake_RZV),
      across(c(N, N_HZ, N_HZ_Hosp, Q_All, C_Med), sum)
    )
  
}

strategy_null()
sims0 <- a_proj_ce(pars_proj, strategy_changeonly)
sims1 <- a_proj_ce(pars_proj, strategy_scheduled65)
sims2 <- a_proj_ce(pars_proj, strategy_scheduled)


dis <- 0.035

sims0
sims1
sims2


sims0 %>% 
  left_join(sims1, by = c("Key", "Year"), suffix = c("_0", "_1")) %>% 
  group_by(Key) %>% 
  mutate(
    d = 1 / (1 + dis) ^ (Year - 2024), 
    dN_Vac = N_Vac_1 - N_Vac_0,
    dN_HZ = N_HZ_1 - N_HZ_0,
    dN_HZ_Hosp = N_HZ_Hosp_1 - N_HZ_Hosp_0,
    dQ_All = Q_All_1 - Q_All_0,
    dC_Med = C_Med_1 - C_Med_0,
    across(c(dN_Vac, dN_HZ, dN_HZ_Hosp, dQ_All, dC_Med), \(x) cumsum(x * d)),
    Thres20 = (dQ_All * 20000 - dC_Med) / dN_Vac,
    Thres30 = (dQ_All * 30000 - dC_Med) / dN_Vac
  ) %>% 
  group_by(Year) %>% 
  summarise(
    across(starts_with("d"), mean),
    Thres = pmin(median(Thres20, na.rm = T), quantile(Thres30, 0.1, na.rm = T))
  ) %>% 
  ggplot() +
  geom_line(aes(x = Year, y = Thres))


sims1 %>% 
  left_join(sims2, by = c("Key", "Year"), suffix = c("_0", "_1")) %>% 
  group_by(Key) %>% 
  mutate(
    d = 1 / (1 + dis) ^ (Year - 2024), 
    dN_Vac = N_Vac_1 - N_Vac_0,
    dN_HZ = N_HZ_1 - N_HZ_0,
    dN_HZ_Hosp = N_HZ_Hosp_1 - N_HZ_Hosp_0,
    dQ_All = Q_All_1 - Q_All_0,
    dC_Med = C_Med_1 - C_Med_0,
    across(c(dN_Vac, dN_HZ, dN_HZ_Hosp, dQ_All, dC_Med), \(x) cumsum(x * d)),
    Thres20 = (dQ_All * 20000 - dC_Med) / dN_Vac,
    Thres30 = (dQ_All * 30000 - dC_Med) / dN_Vac
  ) %>% 
  group_by(Year) %>% 
  summarise(
    across(starts_with("d"), mean),
    Thres = pmin(median(Thres20, na.rm = T), quantile(Thres30, 0.1, na.rm = T))
  ) %>% 
  ggplot() +
  geom_line(aes(x = Year, y = Thres))



