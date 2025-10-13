library(tidyverse)
library(targets)

source(here::here("R", "misc.R"))


ages <- seq(60, 95, 5)


sts <- tar_read(stats_uv, 1)


res_epi <- bind_rows(
  sts[[1]] %>% 
    filter(Age0 >= 80 | Arm != "RZV_1d") %>% 
    filter(Age0 %in% ages) %>% 
    pivot_longer(A:U) %>% 
    pivot_wider(names_from = Index) %>% 
    select(Age = Age0, Arm, Stat = name, N0, N_Vac = N_VacRZV_d, starts_with("Risk")),
  sts[[2]] %>% 
    filter(Age0 >= 80 | Arm == "RZV_2d") %>% 
    filter(Age0 %in% ages) %>% 
    pivot_longer(A:U) %>%
    pivot_wider(names_from = Index) %>% 
    mutate(
      Arm = ifelse(Arm == "RZV_1d", "RZV_1d - SOC", "RZV_2d - SOC")
    ) %>% 
    rename_if(startsWith(names(.), "d"), \(x) gsub("^d", "", x)) %>% 
    select(Age = Age0, Arm, Stat = name, N0, N_Vac = N_VacRZV_d, starts_with("Risk"))
)



res_ce <- lapply(1:3, \(i) {
  sts <- tar_read(stats_uv, i)
  
  bind_rows(
    sts[[1]] %>% 
      filter(Age0 >= 80 | Arm != "RZV_1d") %>% 
      filter(Age0 %in% ages) %>% 
      pivot_longer(A:U) %>% 
      pivot_wider(names_from = Index) %>% 
      mutate(
        Thres = NA
      ) %>% 
      select(Age = Age0, Arm, Stat = name, Q_HZ_d, Q_Life_d, Q_All_d, C_Hosp_d, C_GP_d, C_Med_d, Thres),
    sts[[2]] %>% 
      filter(Age0 >= 80 | Arm == "RZV_2d") %>% 
      filter(Age0 %in% ages) %>% 
      pivot_longer(A:U) %>%
      pivot_wider(names_from = Index) %>% 
      mutate(
        Thres = (dQ_All_d * 2e4 - dC_Med_d) / dN_VacRZV_d
      ) %>% 
      rename_if(startsWith(names(.), "d"), \(x) gsub("^d", "", x)) %>% 
      mutate(
        Arm = ifelse(Arm == "RZV_1d", "RZV_1d - SOC", "RZV_2d - SOC")
      ) %>% 
      select(Age = Age0, Arm, Stat = name, Q_HZ_d, Q_Life_d, Q_All_d, C_Hosp_d, C_GP_d, C_Med_d, Thres)
  )
})
  
  

sts <- tar_read(stats_re, 1)


res_epi_r <- bind_rows(
  sts[[1]] %>%
    filter(Scenario != "Overall" & Arm != "SOC") %>% 
    group_by(Scenario) %>% mutate(Age1 = mean(Age1, na.rm = T)) %>% ungroup() %>% 
    filter(Age0 == 70) %>% 
    filter(Age1 %in% ages) %>% 
    pivot_longer(A:U) %>% 
    pivot_wider(names_from = Index) %>% 
    mutate(
      Arm = case_when(
        Arm == "ReRZV_1d" ~ "ReRZV_1d",
        Arm == "ReRZV_2d" ~ "ReRZV_2d",
        T ~ "ZVL only"
      )
    ) %>% 
    select(Age = Age1, Arm, Stat = name, N0, N_Vac = N_VacRZV_d, starts_with("Risk")),
  sts[[2]] %>% 
    select(- N0) %>% 
    filter(Scenario != "Overall") %>%
    filter(Age0 == 70) %>% 
    filter(Age1 %in% ages) %>% 
    pivot_longer(A:U) %>%
    pivot_wider(names_from = Index) %>% 
    rename_if(startsWith(names(.), "d"), \(x) gsub("^d", "", x)) %>% 
    mutate(
      Arm = ifelse(Arm == "ReRZV_1d", "ReRZV_1d - ZVL only", "ReRZV_2d - ZVL only")
    ) %>% 
    select(Age = Age1, Arm, Stat = name, N0, N_Vac = N_VacRZV_d, starts_with("Risk"))
)


res_ce_r <- lapply(1:3, \(i) {
  sts <- tar_read(stats_re, i)
  
  bind_rows(
    sts[[1]] %>%
      filter(Scenario != "Overall" & Arm != "SOC") %>% 
      group_by(Scenario) %>% mutate(Age1 = mean(Age1, na.rm = T)) %>% ungroup() %>% 
      filter(Age0 == 70) %>% 
      filter(Age1 %in% ages) %>% 
      pivot_longer(A:U) %>% 
      pivot_wider(names_from = Index) %>% 
      mutate(
        Thres = NA,
        Arm = case_when(
          Arm == "ReRZV_1d" ~ "ReRZV_1d",
          Arm == "ReRZV_2d" ~ "ReRZV_2d",
          T ~ "ZVL only"
        )
      ) %>% 
      select(Age = Age1, Arm, Stat = name, Q_HZ_d, Q_Life_d, Q_All_d, C_Hosp_d, C_GP_d, C_Med_d, Thres),
    sts[[2]] %>% 
      select(- N0) %>% 
      filter(Scenario != "Overall") %>%
      filter(Age0 == 70) %>% 
      filter(Age1 %in% ages) %>% 
      pivot_longer(A:U) %>%
      pivot_wider(names_from = Index) %>% 
      mutate(
        Thres = (dQ_All_d * 2e4 - dC_Med_d) / dN_VacRZV_d
      ) %>% 
      rename_if(startsWith(names(.), "d"), \(x) gsub("^d", "", x)) %>% 
      mutate(
        Arm = ifelse(Arm == "ReRZV_1d", "ReRZV_1d - ZVL only", "ReRZV_2d - ZVL only")
      ) %>% 
      select(Age = Age1, Arm, Stat = name, Q_HZ_d, Q_Life_d, Q_All_d, C_Hosp_d, C_GP_d, C_Med_d, Thres)
  )
})


res <- bind_rows(
  res_epi %>% 
    left_join(res_ce[[1]], by = c("Age", "Arm", "Stat")) %>% 
    left_join(res_ce[[3]], by = c("Age", "Arm", "Stat"), suffix = c("", "_3.5%")) %>% 
    left_join(res_ce[[2]], by = c("Age", "Arm", "Stat"), suffix = c("", "_1.5%")),
  res_epi_r %>% 
    left_join(res_ce_r[[1]], by = c("Age", "Arm", "Stat")) %>% 
    left_join(res_ce_r[[3]], by = c("Age", "Arm", "Stat"), suffix = c("", "_3.5%")) %>% 
    left_join(res_ce_r[[2]], by = c("Age", "Arm", "Stat"), suffix = c("", "_1.5%"))
)



dir.create("out", showWarnings = F)


res %>% pull(Arm) %>% unique()

res %>% 
  mutate(
    Arm = factor(Arm, c("SOC", "RZV_2d", "RZV_2d - SOC", "RZV_1d", "RZV_1d - SOC", 
                        "ZVL only", "ReRZV_2d", "ReRZV_2d - ZVL only", "ReRZV_1d", "ReRZV_1d - ZVL only"))
  ) %>% 
  arrange(Age, Arm) %>% 
  write_csv(here::here("out", "DHSC.csv"))



