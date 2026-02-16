

summarise_cohort_rev <- function(yss, wtp_m = 2e4, wtp_h = 3e4) {
  require(tidyverse)
  
  yss_diff <- local({
    temp <- yss %>% 
      pivot_longer(-c(Scenario, Age0, Arm, Key, N0, Year0), names_to = "Index")
    
    temp %>% 
      filter(Arm != "SOC") %>% 
      left_join(
        temp %>% 
          filter(Arm == "SOC") %>% 
          select(Scenario, Age0, Key, Index, value0 = value),
        relationship = "many-to-many"
      ) %>% 
      mutate(
        Index = paste0("d", Index),
        Diff = value - value0
      ) %>% 
      select(-value, -value0) %>% 
      pivot_wider(names_from = Index, values_from = "Diff") %>% 
      mutate(
        ICER = dC_All_d / dQ_All_d,
        Price0 = dC_VacRZV_d / dN_VacRZV_d,
        Thres20 = (dQ_All_d * wtp_m - dC_Med_d) / dN_VacRZV_d,
        Thres30 = (dQ_All_d * wtp_h - dC_Med_d) / dN_VacRZV_d,
      )
    
  })
  
  ss <- list()
  
  ss$stats_uv_ys <- yss %>% 
    group_by(Scenario, Age0, Arm) %>% 
    select(-Key) %>% 
    summarise_all(amlu) %>% 
    pivot_longer(ends_with(c("_A", "_M", "_L", "_U")), names_pattern = "(\\S+)_(A|M|L|U)", names_to = c("Index", "name")) %>% 
    pivot_wider() %>% 
    ungroup()
  
  
  
  ss$stats_uv_ce <- yss_diff %>% 
    group_by(Scenario, Age0, Arm, Year0, N0) %>% 
    select(-Key) %>% 
    summarise(
      across(everything(), amlu),
      Thres20_50_M = median(Thres20),
      Thres30_90_M = quantile(Thres30, 0.1)
    ) %>% 
    mutate(
      Thres_M = pmin(Thres20_50_M, Thres30_90_M)
    ) %>% 
    pivot_longer(ends_with(c("_A", "_M", "_L", "_U")), names_pattern = "(\\S+)_(A|M|L|U)", names_to = c("Index", "name")) %>% 
    pivot_wider() %>% 
    ungroup()
  
  return(ss)
}


summarise_cohort_rev_re <- function(yss, wtp_m = 2e4, wtp_h = 3e4) {
  require(tidyverse)
  
  
  ### Comparison -----
  yss_diff <- local({
    temp <- yss %>% 
      pivot_longer(-c(Scenario, Age0, Age1, Arm, Key, N0, Year0), names_to = "Index")
    
    dy0 <- temp %>% 
      filter(Scenario == "Overall") %>%
      filter(Arm != "SOC") %>% 
      left_join(
        temp %>% 
          filter(Arm == "SOC") %>% 
          filter(Scenario == "Overall") %>%
          select(Scenario, Age0, Key, Index, value0 = value),
        relationship = "many-to-many"
      )
    
    dy1 <- temp %>% 
      filter(Scenario != "Overall") %>%
      filter(Arm %in% c("ReRZV_1d", "ReRZV_2d")) %>% 
      left_join(
        temp %>% 
          filter(Arm == "Vac") %>% 
          filter(Scenario != "Overall") %>%
          select(Scenario, Age0, Key, Index, value0 = value),
        relationship = "many-to-many"
      )
    
    bind_rows(dy0, dy1) %>% 
      mutate(
        Index = paste0("d", Index),
        Diff = value - value0
      ) %>% 
      select(-value, -value0) %>% 
      pivot_wider(names_from = Index, values_from = "Diff") %>% 
      mutate(
        across(starts_with("d"), \(x) x/N0),
        ICER = dC_All_d / dQ_All_d,
        Price0 = dC_VacRZV_d / dN_VacRZV_d,
        Thres20 = (dQ_All_d * wtp_m - dC_Med_d) / dN_VacRZV_d,
        Thres30 = (dQ_All_d * wtp_h - dC_Med_d) / dN_VacRZV_d,
      )
  })
  
  
  ss <- list()
  
  ss$stats_re_ys <- yss %>% 
    group_by(Scenario, Arm, Age0, Age1) %>% 
    select(-Key) %>% 
    summarise_all(amlu) %>% 
    pivot_longer(ends_with(c("_A", "_M", "_L", "_U")), names_pattern = "(\\S+)_(A|M|L|U)", names_to = c("Index", "name")) %>% 
    pivot_wider() %>% 
    ungroup()
  
  
  ss$stats_re_ce <- yss_diff %>% 
    #select(-N0) %>% 
    group_by(Scenario, Age0, Age1, Arm, Year0) %>% 
    select(-Key) %>% 
    summarise(
      N0 = mean(N0),
      across(everything(), amlu),
      Thres20_50_M = median(Thres20, na.rm = T),
      Thres30_90_M = quantile(Thres30, 0.1, na.rm = T)
    ) %>% 
    mutate(
      Thres_M = pmin(Thres20_50_M, Thres30_90_M)
    ) %>% 
    pivot_longer(ends_with(c("_A", "_M", "_L", "_U")), names_pattern = "(\\S+)_(A|M|L|U)", names_to = c("Index", "name")) %>% 
    pivot_wider() %>% 
    ungroup()
  
  return(ss)
}
