

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
          tp_uptake == "Ini" ~ p_initial,
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
  
  # population <- with(model_proj, {
  #   yr <- year0
  #   pop <- tibble(Year = yr, Age = age0:age1, Vaccine = "None", TimeVac = -1, Prop = 1) %>% 
  #     crossing(Key = 1:n_sims) %>% 
  #     left_join(p_uptake, by = "Key")
  #   
  #   pop <- pop %>% uptaking(yr = yr, strategy = strategy)
  #   
  #   collector <- pop %>% fn_vac()
  #   
  #   while (yr < year1) {
  #     yr <- yr + 1
  #     setTxtProgressBar(pb, yr)
  #     pop <- pop %>% 
  #       ageing(yr, age0, age1) %>% 
  #       uptaking(yr = yr, strategy = strategy)
  #     collector <- collector %>% bind_rows(pop %>% fn_vac())
  #   }
  #   collector
  # })
  # 
  
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
      N_HZ_Death = p_mor_hz * N_HZ
    ) %>% 
    select(-starts_with(c("r_", "p_")))
  
  
  yss_agp <- sims %>% 
    filter(Age >= 60 & Age < 100) %>% 
    mutate(
      Agp = cut(Age, seq(60, 100, 5), right = F)
    ) %>% 
    group_by(Key, Year, Agp) %>% 
    summarise(
      Protection = weighted.mean(Protection, N),
      across(starts_with("N"), sum)
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
    ungroup()
  
  yss_68 <- sims %>% 
    filter(Age >= 60 & Age < 100) %>% 
    mutate(
      Agp = ifelse(Age < 80, "60_80", "80+")
    ) %>% 
    group_by(Key, Year, Agp) %>% 
    summarise(
      Protection = weighted.mean(Protection, N),
      across(starts_with("N"), sum)
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
    ungroup()
  
  yss_all <- sims %>% 
    filter(Age >= 60 & Age < 100) %>% 
    group_by(Key, Year) %>% 
    summarise(
      Protection = weighted.mean(Protection, N),
      across(starts_with("N"), sum)
    ) %>% 
    mutate(
      Agp = "All",
      Coverage_ZVL = N_Covered_ZVL / N,
      Coverage_RZV = N_Covered_RZV / N,
      Coverage = Coverage_ZVL + Coverage_RZV,
      IncR_HZ = N_HZ / N,
      IncR_HZ_Hosp = N_HZ_Hosp  / N,
      IncR_HZ_PHN = N_HZ_PHN / N,
      MorR_HZ = N_HZ_Death / N
    ) %>% 
    ungroup()
  
  return(list(
    vtype = pars$vtype,
    Yss_Agp = yss_agp,
    Yss_68 = yss_68,
    Yss_All = yss_all
  ))
  
}


exec_projection <- function(pars, year1 = 2050, n_sims = 5000) {
  age0 <- 60 
  res <- list(
    "Null" = strategy_null,
    "Stay" = strategy_zvl,
    "ToRZV" = strategy_changeonly,
    "Sch65" = strategy_scheduled65,
    "Sch" = strategy_scheduled,
    "Sch1d85" = strategy_scheduled_1d85,
    "Sch1d95" = strategy_scheduled_1d95,
    "Sch2d85" = strategy_scheduled_2d85,
    "Sch2d95" = strategy_scheduled_2d95
  )
  
  res <- lapply(res, \(strategy) a_projection(pars, strategy, year0 = 2013, year1 = year1, age0 = age0, age1 = 100, n_sims = n_sims))
  
  yss_agp = lapply(names(res), \(k) res[[k]]$Yss_Agp %>% mutate(Scenario = k)) %>% bind_rows() 
  yss_68 = lapply(names(res), \(k) res[[k]]$Yss_68 %>% mutate(Scenario = k)) %>% bind_rows() 
  yss_all = lapply(names(res), \(k) res[[k]]$Yss_All %>% mutate(Scenario = k)) %>% bind_rows() 
  
  return(list(
    vtype = pars$vtype,
    Year1 = year1,
    Yss_Agp = yss_agp,
    Yss_68 = yss_68,
    Yss_All = yss_all
  ))
}



