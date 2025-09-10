library(tidyverse)
library(targets)





ages <- seq(60, 95, 5)

sts <- tar_read(stats_uv, 1)

res00 <- bind_rows(
  sts[[1]] %>% 
    filter(Age0 >= 80 | Arm != "RZV_1d") %>% 
    filter(Age0 %in% ages) %>% 
    pivot_longer(A:U) %>% 
    pivot_wider(names_from = Index) %>% 
    mutate(
      QL_HZ_d = - Q_HZ_d,
      QL_Life_d = - Q_Life_d,
      QL_All_d = QL_HZ_d + QL_Life_d,
      Thres = NA
    ) %>% 
    select(Age = Age0, Arm, Stat = name, N0, N_Vac = N_VacRZV_d, starts_with("Risk"), starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres),
  sts[[2]] %>% 
    filter(Age0 >= 80 | Arm == "RZV_2d") %>% 
    filter(Age0 %in% ages) %>% 
    pivot_longer(A:U) %>%
    pivot_wider(names_from = Index) %>% 
    mutate(
      dQL_HZ_d = - dQ_HZ_d,
      dQL_Life_d = - dQ_Life_d,
      dQL_All_d = dQL_HZ_d + dQL_Life_d,
      Thres = (dQ_All_d * 2e4 - dC_Med_d) / dN_VacRZV_d
    ) %>% 
    rename_if(startsWith(names(.), "d"), \(x) gsub("^d", "", x)) %>% 
    mutate(
      Arm = ifelse(Arm == "RZV_1d", "RZV_1d - SOC", "RZV_2d - SOC")
    ) %>% 
    select(Age = Age0, Arm, Stat = name, N0, N_Vac = N_VacRZV_d, starts_with("Risk"), starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres)
)


sts <- tar_read(stats_uv, 2)

res15 <- bind_rows(
  sts[[1]] %>% 
    filter(Age0 >= 80 | Arm != "RZV_1d") %>% 
    filter(Age0 %in% ages) %>% 
    pivot_longer(A:U) %>% 
    pivot_wider(names_from = Index) %>% 
    mutate(
      QL_HZ_d = - Q_HZ_d,
      QL_Life_d = - Q_Life_d,
      QL_All_d = QL_HZ_d + QL_Life_d,
      Thres = NA
    ) %>% 
    select(Age = Age0, Arm, Stat = name, starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres),
  sts[[2]] %>% 
    filter(Age0 >= 80 | Arm == "RZV_2d") %>% 
    filter(Age0 %in% ages) %>% 
    pivot_longer(A:U) %>%
    pivot_wider(names_from = Index) %>% 
    mutate(
      dQL_HZ_d = - dQ_HZ_d,
      dQL_Life_d = - dQ_Life_d,
      dQL_All_d = dQL_HZ_d + dQL_Life_d,
      Thres = (dQ_All_d * 2e4 - dC_Med_d) / dN_VacRZV_d
    ) %>% 
    rename_if(startsWith(names(.), "d"), \(x) gsub("^d", "", x)) %>% 
    mutate(
      Arm = ifelse(Arm == "RZV_1d", "RZV_1d - SOC", "RZV_2d - SOC")
    ) %>% 
    select(Age = Age0, Arm, Stat = name, starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres)
)


sts <- tar_read(stats_uv, 3)

res35 <- bind_rows(
  sts[[1]] %>% 
    filter(Age0 >= 80 | Arm != "RZV_1d") %>% 
    filter(Age0 %in% ages) %>% 
    pivot_longer(A:U) %>% 
    pivot_wider(names_from = Index) %>% 
    mutate(
      QL_HZ_d = - Q_HZ_d,
      QL_Life_d = - Q_Life_d,
      QL_All_d = QL_HZ_d + QL_Life_d,
      Thres = NA
    ) %>% 
    select(Age = Age0, Arm, Stat = name, starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres),
  sts[[2]] %>% 
    filter(Age0 >= 80 | Arm == "RZV_2d") %>% 
    filter(Age0 %in% ages) %>% 
    pivot_longer(A:U) %>%
    pivot_wider(names_from = Index) %>% 
    mutate(
      dQL_HZ_d = - dQ_HZ_d,
      dQL_Life_d = - dQ_Life_d,
      dQL_All_d = dQL_HZ_d + dQL_Life_d,
      Thres = (dQ_All_d * 2e4 - dC_Med_d) / dN_VacRZV_d
    ) %>% 
    rename_if(startsWith(names(.), "d"), \(x) gsub("^d", "", x)) %>% 
    mutate(
      Arm = ifelse(Arm == "RZV_1d", "RZV_1d - SOC", "RZV_2d - SOC")
    ) %>% 
    select(Age = Age0, Arm, Stat = name, starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres)
)


sts <- tar_read(stats_re, 1)

res00r <- bind_rows(
  sts[[1]] %>%
    filter(Scenario != "Overall" & Arm != "SOC") %>% 
    group_by(Scenario) %>% mutate(Age1 = mean(Age1, na.rm = T)) %>% ungroup() %>% 
    filter(Age0 == 70) %>% 
    filter(Age1 %in% ages) %>% 
    pivot_longer(A:U) %>% 
    pivot_wider(names_from = Index) %>% 
    mutate(
      QL_HZ_d = - Q_HZ_d,
      QL_Life_d = - Q_Life_d,
      QL_All_d = QL_HZ_d + QL_Life_d,
      Thres = NA,
      Arm = case_when(
        Arm == "ReRZV_1d" ~ "ReRZV_1d",
        Arm == "ReRZV_2d" ~ "ReRZV_2d",
        T ~ "ZVL only"
      )
    ) %>% 
    select(Age = Age1, Arm, Stat = name, N0, N_Vac = N_VacRZV_d, starts_with("Risk"), starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres),
  sts[[2]] %>% 
    select(- N0) %>% 
    filter(Scenario != "Overall") %>%
    filter(Age0 == 70) %>% 
    filter(Age1 %in% ages) %>% 
    pivot_longer(A:U) %>%
    pivot_wider(names_from = Index) %>% 
    filter(name != "A") %>% 
    mutate(
      dQL_HZ_d = - dQ_HZ_d,
      dQL_Life_d = - dQ_HZ_Norm_d,
      dQL_All_d = dQL_HZ_d + dQL_Life_d,
      Thres = (dQ_All_d * 2e4 - dC_Med_d) / dN_VacRZV_d
    ) %>% 
    rename_if(startsWith(names(.), "d"), \(x) gsub("^d", "", x)) %>% 
    mutate(
      Arm = ifelse(Arm == "ReRZV_1d", "ReRZV_1d - ZVL only", "ReRZV_2d - ZVL only")
    ) %>% 
    select(Age = Age1, Arm, Stat = name, N0, N_Vac = N_VacRZV_d, starts_with("Risk"), starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres)
)


sts <- tar_read(stats_re, 2)

res15r <- bind_rows(
  sts[[1]] %>%
    filter(Scenario != "Overall" & Arm != "SOC") %>% 
    group_by(Scenario) %>% mutate(Age1 = mean(Age1, na.rm = T)) %>% ungroup() %>% 
    filter(Age0 == 70) %>% 
    filter(Age1 %in% ages) %>% 
    pivot_longer(A:U) %>% 
    pivot_wider(names_from = Index) %>% 
    mutate(
      QL_HZ_d = - Q_HZ_d,
      QL_Life_d = - Q_Life_d,
      QL_All_d = QL_HZ_d + QL_Life_d,
      Thres = NA,
      Arm = case_when(
        Arm == "ReRZV_1d" ~ "ReRZV_1d",
        Arm == "ReRZV_2d" ~ "ReRZV_2d",
        T ~ "ZVL only"
      )
    ) %>% 
    select(Age = Age1, Arm, Stat = name, N_Vac = N_VacRZV_d, starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres),
  sts[[2]] %>% 
    select(- N0) %>% 
    filter(Scenario != "Overall") %>%
    filter(Age0 == 70) %>% 
    filter(Age1 %in% ages) %>% 
    pivot_longer(A:U) %>%
    pivot_wider(names_from = Index) %>% 
    mutate(
      dQL_HZ_d = - dQ_HZ_d,
      dQL_Life_d = - dQ_Life_d,
      dQL_All_d = dQL_HZ_d + dQL_Life_d,
      Thres = (dQ_All_d * 2e4 - dC_Med_d) / dN_VacRZV_d
    ) %>% 
    rename_if(startsWith(names(.), "d"), \(x) gsub("^d", "", x)) %>% 
    mutate(
      Arm = ifelse(Arm == "ReRZV_1d", "ReRZV_1d - ZVL only", "ReRZV_2d - ZVL only")
    ) %>% 
    select(Age = Age1, Arm, Stat = name, N_Vac = N_VacRZV_d, starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres)
)


sts <- tar_read(stats_re, 3)

res35r <- bind_rows(
  sts[[1]] %>%
    filter(Scenario != "Overall" & Arm != "SOC") %>% 
    group_by(Scenario) %>% mutate(Age1 = mean(Age1, na.rm = T)) %>% ungroup() %>% 
    filter(Age0 == 70) %>% 
    filter(Age1 %in% ages | Arm == "Vac") %>% 
    pivot_longer(A:U) %>% 
    pivot_wider(names_from = Index) %>% 
    mutate(
      QL_HZ_d = - Q_HZ_d,
      QL_Life_d = - Q_Life_d,
      QL_All_d = QL_HZ_d + QL_Life_d,
      Thres = NA,
      Arm = case_when(
        Arm == "ReRZV_1d" ~ "ReRZV_1d",
        Arm == "ReRZV_2d" ~ "ReRZV_2d",
        T ~ "ZVL only"
      )
    )  %>% 
    select(Age = Age1, Arm, Stat = name, N_Vac = N_VacRZV_d, starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres),
  sts[[2]] %>% 
    select(- N0) %>% 
    filter(Scenario != "Overall") %>%
    filter(Age0 == 70) %>% 
    filter(Age1 %in% ages | Arm == "Vac") %>% 
    pivot_longer(A:U) %>%
    pivot_wider(names_from = Index) %>% 
    mutate(
      dQL_HZ_d = - dQ_HZ_d,
      dQL_Life_d = - dQ_Life_d,
      dQL_All_d = dQL_HZ_d + dQL_Life_d,
      Thres = (dQ_All_d * 2e4 - dC_Med_d) / dN_VacRZV_d
    ) %>% 
    rename_if(startsWith(names(.), "d"), \(x) gsub("^d", "", x)) %>% 
    mutate(
      Arm = case_when(
        Arm == "ReRZV_1d" ~ "ReRZV_1d - ZVL only",
        Arm == "ReRZV_2d" ~ "ReRZV_2d - ZVL only",
        T ~ ""
      )
    ) %>% 
    select(Age = Age1, Arm, Stat = name, N_Vac = N_VacRZV_d, starts_with("QL"), C_Hosp_d, C_GP_d, C_Med_d, Thres)
)




res <- bind_rows(
  res00 %>%
    left_join(res35, by = c("Age", "Arm", "Stat"), suffix = c("_0%", "_3.5%")) %>%
    left_join(res15, by = c("Age", "Arm", "Stat"), suffix = c("", "_1.5%")),
  res00r %>%
    left_join(res35r, by = c("Age", "Arm", "Stat"), suffix = c("_0%", "_3.5%")) %>%
    left_join(res15r, by = c("Age", "Arm", "Stat"), suffix = c("", "_1.5%"))
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



