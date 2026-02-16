library(targets)
library(tidyverse)


strategy_changeonly <- function(df, year) {
  require(tidyverse)
  
  df <- df %>% 
    mutate(
      eli = case_when(
        Vaccine != "None" ~ NA,
        Age < 70 ~ NA,
        Age < 80 ~ ifelse(year < 2023, "ZVL", "RZV_2d"),
        T ~ NA
      ),
      tp_uptake = case_when(
        is.na(eli) ~ NA,
        Age == 70 ~ "Ini70",
        T ~ "Cat"
      )
    )
  
  return(df)
}


strategy_scheduled <- function(df, year) {
  require(tidyverse)
  
  if (year < 2023) {
    df <- df %>% 
      mutate(
        eli = case_when(
          Vaccine != "None" ~ NA,
          Age < 70 ~ NA,
          Age < 80 ~ "ZVL",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age == 70 ~ "Ini70",
          T ~ "Cat"
        )
      )
  } else if (year < 2028) {
    df <- df %>% 
      mutate(
        eli = case_when(
          Vaccine != "None" ~ NA,
          Age < 65 ~ NA,
          Age >= 80 ~ NA,
          Age >= 70 ~ "RZV_2d",
          Age < (65 + year - 2023) ~ "RZV_2d",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age == 65 ~ "Ini65",
          Age == 70 ~ "Ini70",
          T ~ "Cat"
        )
      )
  } else if (year < 2033) {
    df <- df %>% 
      mutate(
        eli = case_when(
          Vaccine != "None" ~ NA,
          Age < 60 ~ NA,
          Age >= 80 ~ NA,
          Age >= 65 ~ "RZV_2d",
          Age < (60 + year - 2028) ~ "RZV_2d",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age == 60 ~ "Ini60",
          Age == 65 ~ "Ini65",
          T ~ "Cat"
        )
      )
  } else {
    df <- df %>% 
      mutate(
        eli = case_when(
          Vaccine != "None" ~ NA,
          Age < 60 ~ NA,
          Age < 80 ~ "RZV_2d",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age == 60 ~ "Ini60",
          T ~ "Cat"
        )
      )
  }
}


strategy_scheduled65 <- function(df, year) strategy_scheduled(df, min(year, 2027))



model_proj <- list()
class(model_proj) <- "model_proj"


model_proj$uptaking <- function(df, yr, strategy) {
  df <- df %>% strategy(yr)
  if ((df %>% filter(!is.na(eli)) %>% nrow()) <= 0) {
    return(df)
  }
  df <- bind_rows(
    df %>% filter(is.na(eli)),
    df %>% filter(!is.na(eli)) %>% 
      mutate(
        p_uptake = case_when(
          tp_uptake == "Cat" ~ p_catch,
          tp_uptake == "Ini60" ~ p_initial_60,
          tp_uptake == "Ini65" ~ p_initial_65,
          tp_uptake == "Ini70" ~ p_initial_70,
          T ~ 0
        ),
        Prop_n = Prop * (1 - p_uptake), 
        Prop_v = Prop * p_uptake
      ) %>% 
      pivot_longer(c(Prop_n, Prop_v)) %>% 
      mutate(
        Prop = value,
        TimeVac = case_when(
          name == "Prop_n" ~ TimeVac,
          Vaccine == eli ~ TimeVac,
          T ~ 1
        ),
        Vaccine = ifelse(name == "Prop_n", Vaccine, eli)
      ) %>% 
      select(-c(name, value))
  ) %>% 
    select(-c(eli, p_uptake, tp_uptake)) %>% 
    arrange(Year, Age)
  
  return(df)
}


model_proj$ageing <- function(df, yr, age0, age1) {
  ks <- df %>% pull(Key) %>% unique()
  
  new_in <- df %>% filter(Age == age0) %>% 
    group_by(Key) %>%
    filter(row_number()==1) %>% 
    mutate(Year = yr, Vaccine = "None", TimeVac = -1, Prop = 1) %>% 
    ungroup()
  
  
  df <- df %>% 
    mutate(
      Age = Age + 1,
      TimeVac = ifelse(TimeVac > 0, TimeVac + 1, TimeVac),
      Year = yr
    ) %>% 
    filter(Age <= age1) %>%       
    bind_rows(new_in) %>% 
    arrange(Age)
  
  return(df)
}


a_projection <- function(pars, strategy, year0 = 2013, year1 = 2050, age0 = 60, age1 = 100, n_sims = 5000) {
  require(tidyverse)
  
  # pars <- tar_read(pars_proj)
  # prefix <- gsub("Year0", "", names(pars)[1])
  # names(pars) <- gsub(prefix, "", names(pars))
  
  p_uptake <- pars$Uptake
  p_demo <- pars$Demography
  p_ce <- pars$CostEff
  
  ves_rzv <- bind_rows(
    pars$VE_RZV_2d,
    pars$VE_RZV_1d,
    pars$VE_ReRZV_2d,
    pars$VE_ReRZV_1d,
  )
  
  ves_zvl <- pars$VE_ZVL
  
  
  fn_vac <- function(df) {
    df <- df %>% 
      left_join(p_demo$N, by = c("Year", "Age")) %>%
      mutate(
        N = N * Prop
      ) %>% 
      select(-Prop)
    
    bind_rows(
      df %>% 
        filter(Vaccine == "ZVL") %>% 
        left_join(ves_zvl, by = c("Key", "Age", "Vaccine", "TimeVac")),
      df %>% 
        filter(Vaccine == "None") %>% 
        mutate(Protection = 0),
      df %>% 
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
        across(starts_with("N"), sum), .groups="keep"
      ) %>% 
      ungroup()
  }
  
  
  n_sims <- min(pars$N_Sims, n_sims)
  
  pb <- txtProgressBar(min = year0, max = year1, style = 3,  width = 50, char = "=") 

  
  population <- with(model_proj, {
    yr <- year0
    pop <- tibble(Year = yr, Age = age0:age1, Vaccine = "None", TimeVac = -1, Prop = 1) %>% 
      crossing(Key = 1:n_sims) %>% 
      left_join(p_uptake, by = "Key")
    
    pop <- pop %>% uptaking(yr = yr, strategy = strategy)
    
    collector <- pop
    
    while (yr < year1) {
      yr <- yr + 1
      setTxtProgressBar(pb, yr)
      pop <- pop %>% 
        ageing(yr, age0, age1) %>% 
        uptaking(yr = yr, strategy = strategy)
      collector <- collector %>% bind_rows(pop)
    }
    collector
  }) %>% fn_vac()
  
  
  sims <- population %>% 
    left_join(p_demo$DeathIm %>% select(Year, Age, r_mor_bg = r_death), by = c("Year", "Age")) %>% 
    inner_join(pars$Epidemiology, by = c("Key", "Age")) %>% 
    mutate(
      r_mor = (1 - Protection) * r_mor_hz + r_mor_bg,
      p_mor_hz = (1 - Protection) * r_mor_hz * (1 - exp(- r_mor)) / r_mor,
      N_HZ = (1 - Protection) * r_hz * N,
      N_HZ_GP = p_gp * N_HZ,
      N_HZ_Hosp = N_HZ - N_HZ_GP,
      N_HZ_PHN = p_phn * N_HZ,
      N_HZ_PHN_GP = p_gp * N_HZ_PHN,
      N_HZ_Death = p_mor_hz * N_HZ,
      QALY = 0
    ) %>% 
    left_join(pars$CostEff, by = c("Key", "Age")) %>% 
    mutate(
      Q_Life = QALY,
      Q_HZ = - N_HZ * QL_ph,
      Q_All = Q_Life + Q_HZ,
      C_Hosp = N_HZ_Hosp * cost_hosp_pp_inf,
      C_GP_NonPHN = (N_HZ_GP - N_HZ_PHN_GP) * cost_GP_pp_non_PHN_HZ_inf,
      C_GP_PHN = N_HZ_PHN_GP * cost_GP_pp_PHN_inf,
      C_GP = C_GP_NonPHN + C_GP_PHN,
      C_Med = C_GP + C_Hosp,
      across(starts_with("C_Vac"), \(x) ifelse(is.na(x), 0, x)),
      across(starts_with(c("C_", "Q_", "N_")), \(x) ifelse(is.na(x), 0, x))
    ) %>% 
    select(-starts_with(c("r_", "p_")))
  
  
  yss_agp <- sims %>% 
    filter(N > 0) %>% 
    filter(Age >= 60 & Age < 100) %>% 
    mutate(
      Agp = cut(Age, seq(60, 100, 5), right = F)
    ) %>% 
    group_by(Key, Year, Agp) %>% 
    summarise(
      Protection = weighted.mean(Protection, N),
      across(starts_with(c("C_", "Q_", "N_")), sum),
      N = sum(N)
    ) %>% 
    mutate(
      Coverage_ZVL = N_Covered_ZVL / N,
      Coverage_RZV = N_Covered_RZV / N,
      Coverage = Coverage_ZVL + Coverage_RZV,
      IncR_HZ = N_HZ / N,
      IncR_HZ_Hosp = N_HZ_Hosp  / N,
      IncR_HZ_PHN = N_HZ_PHN / N,
      MorR_HZ = N_HZ_Death / N
    ) %>% 
    filter(Year >= 2023) %>% 
    group_by(Year, Agp) %>% 
    summarise(across(everything(), mean)) %>% 
    ungroup()
  
  return(list(
    vtype = pars$vtype,
    Yss_Agp = yss_agp
  ))
  
}


exec_projection <- function(pars, year1 = 2050, n_sims = 5000) {
  age0 <- 60 
  res <- list(
    "Fr70" = strategy_changeonly,
    "Fr65" = strategy_scheduled65,
    "Fr60" = strategy_scheduled
  )
  
  res <- lapply(res, \(strategy) a_projection(pars, strategy, year0 = 2013, year1 = year1, age0 = age0, age1 = 100, n_sims = n_sims))
  
  yss_agp = lapply(names(res), \(k) res[[k]]$Yss_Agp %>% mutate(Scenario = k)) %>% bind_rows() 
  
  return(list(
    vtype = pars$vtype,
    Year1 = year1,
    Yss_Agp = yss_agp
  ))
}



pars_proj <- tar_read(pars_proj, 1)
tag <- gsub("Year0", "", names(pars_proj)[1])
names(pars_proj) <- gsub(tag, "", names(pars_proj))


pars_proj$Uptake <- pars_proj$Uptake %>% 
  mutate(
    p_initial_60 = 0.301,
    p_initial_65 = 0.339,
    p_initial_70 = 0.405
  )


proj <- exec_projection(pars_proj, n_sims = 2000)


write_csv(proj$Yss_Agp, file = here::here("out", "proj_agp.csv"))



proj$Yss_Agp %>% 
  filter(Year >= 2023) %>%
  ggplot() +
  geom_bar(aes(x = Year, y = N_Uptake_RZV, fill = Agp), position = "stack", stat = "identity", width = 1) +
  geom_vline(xintercept = c(2023, 2028) + 0.5, linetype = 2) + 
  facet_wrap(.~Scenario)


