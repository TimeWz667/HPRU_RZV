sens_ce <- function(yss_uv) { 
  require(tidyverse)
  
  fn_thres <- function(df) {
    df %>% 
      mutate(
        dC_Med_d = dC_Hosp_d + dC_GP_NonPHN_d + dC_GP_PHN_d,
        Thres20 = ((dQ_HZ_d + dQ_Life_d) * 2e4 - dC_Med_d) / dN_VacRZV_d,
        Thres30 = ((dQ_HZ_d + dQ_Life_d) * 3e4 - dC_Med_d) / dN_VacRZV_d,
      ) %>% 
      group_by(Age0) %>% 
      summarise(
        T20_50 = quantile(Thres20, 0.5),
        T30_90 = quantile(Thres30, 0.9), 
        Thres = pmin(T20_50, T30_90)
      )
  }
  
  
  yss <- local({
    temp <- yss_uv %>% 
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
        Thres20 = (dQ_All_d * 2e4 - dC_Med_d) / dN_VacRZV_d,
        Thres30 = (dQ_All_d * 3e4 - dC_Med_d) / dN_VacRZV_d,
      )
    
  }) %>% 
    filter(Arm == "RZV_2d")
  
  
  
  case0 <- yss %>% fn_thres
  
  
  case1 <- bind_rows(
    yss %>% mutate(
      dC_Hosp_d = dC_Hosp_d * (1 + 0.1)
    ) %>% fn_thres %>% mutate(Pars = "C_Hosp", Direction = +0.1),
    yss %>% mutate(
      dC_Hosp_d = dC_Hosp_d * (1 - 0.1)
    ) %>% fn_thres %>% mutate(Pars = "C_Hosp", Direction = -0.1),
    
    yss %>% mutate(
      dC_GP_NonPHN_d = dC_GP_NonPHN_d * (1 + 0.1)
    ) %>% fn_thres %>% mutate(Pars = "C_GP_NonPHN", Direction = +0.1),
    yss %>% mutate(
      dC_GP_NonPHN_d = dC_GP_NonPHN_d * (1 - 0.1)
    ) %>% fn_thres %>% mutate(Pars = "C_GP_NonPHN", Direction = -0.1),
    
    yss %>% mutate(
      dC_GP_PHN_d = dC_GP_PHN_d * (1 + 0.1)
    ) %>% fn_thres %>% mutate(Pars = "C_GP_PHN", Direction = +0.1),
    yss %>% mutate(
      dC_GP_PHN_d = dC_GP_PHN_d * (1 - 0.1)
    ) %>% fn_thres %>% mutate(Pars = "C_GP_PHN", Direction = -0.1),
    
    yss %>% mutate(
      dQ_HZ_d = dQ_HZ_d * (1 + 0.1)
    ) %>% fn_thres %>% mutate(Pars = "Q_HZ", Direction = +0.1),
    yss %>% mutate(
      dQ_HZ_d = dQ_HZ_d * (1 - 0.1)
    ) %>% fn_thres %>% mutate(Pars = "Q_HZ", Direction = -0.1),
    
    yss %>% mutate(
      dQ_Life_d = dQ_Life_d * (1 + 0.1)
    ) %>% fn_thres %>% mutate(Pars = "Q_Life", Direction = +0.1),
    yss %>% mutate(
      dQ_Life_d = dQ_Life_d * (1 - 0.1)
    ) %>% fn_thres %>% mutate(Pars = "Q_Life", Direction = -0.1)
  )
  
  
  lvs <- case1 %>% 
    filter(Age0 == 80) %>% 
    group_by(Pars) %>% 
    summarise(R = diff(range(Thres))) %>% 
    arrange(R) %>% 
    pull(Pars)
  
  sens_ce <- case1 %>% 
    select(Pars, Direction, Thres, Age0) %>% 
    left_join(case0 %>% select(Age0, Thres0 = Thres)) 

  return(sens_ce)
    
}


summarise_sens_ce <- function(sens_ce, prefix = "", folder = NA, ext = ".pdf") {
  require(tidyverse)
  
  if (!is.na(folder)) {
    root_tab <- here::here("docs", "tabs", folder)
    root_fig <- here::here("docs", "figs", folder)
    dir.create(root_tab, showWarnings = F)
    dir.create(root_fig, showWarnings = F)
  } else {
    root_tab <- here::here("docs", "tabs")
    root_fig <- here::here("docs", "figs")
  }
  
  if (prefix != "") {
    prefix_fig <- glue::as_glue("g_") + prefix + "_"
    prefix_tab <- glue::as_glue(prefix) + "_"
  } else {
    prefix_fig <- glue::as_glue("g_")
    prefix_tab <- glue::as_glue("")
  }

  write_csv(sens_ce, here::here(root_tab, prefix_tab + "sens_ce1d.csv"))
  
  
  labs_lvs <- c(
    "C_GP_PHN" = "Cost, GP, PHN",
    "C_GP_NonPHN" = "Cost, GP, No PHN",
    "Q_Life" = "QALY, Survival",
    "C_Hosp" = "Cost, Hospitalisation",     
    "Q_HZ" = "QALY, HZ"
  )
  
  g_sens_ce <- sens_ce %>%
    filter(Age0 %in% c(75, 80, 85)) %>% 
    mutate(
      Pars = factor(Pars, names(labs_lvs)),
      i = as.numeric(Pars),
      D = ifelse(Direction > 0, "+10%", "-10%"),
      a0 = paste0("Age of vaccination: ", Age0)
    ) %>% 
    ggplot() +
    geom_rect(aes(xmin = Thres0, xmax = Thres, ymin = i - 0.45, ymax = i + 0.45, 
                  y = Pars, fill = D)) +
    scale_fill_discrete("Changes") +
    scale_x_continuous("Threshold price, GBP per administration") +
    scale_y_discrete("Components", labels = labs_lvs) +
    facet_wrap(a0~.)
  
  ggsave(g_sens_ce, filename = here::here(root_fig, prefix_fig + "sens_ce1d" + ext), width = 8, height = 4)
  
  
  
}


exec_sens_uni <- function(pars, age0s) {
  # yss_uv <- tar_read(yss_uv, 3)
  # 
  # pars <- tar_read(pars_base, 3)
  # 
  # prefix <- gsub("Year0", "", names(pars)[[1]])
  # names(pars) <- gsub(prefix, "", names(pars))
  # age0s <- seq(60, 90, 5)
  
  pars_m <- pars
  
  to_investigate <- c("VE", "Inc", "Hosp", "PHN", "Mor", "Q_HZ", "C_GP_NonPHN", "C_GP_PHN", "C_Hosp")
  
  n_pars <- length(to_investigate)
  pars_m$N_Sims <- n_pars * 2 + 1
  
  pars_m$key_mapping <- tibble(
    Key = 1:(n_pars * 2 + 1),
    Pars = c(NA, rep(to_investigate, each = 2)),
    Direction = c(NA, rep(c("Q025", "Q975"), n_pars))
  )
  
  sts <- pars$Epidemiology %>% 
    group_by(Age) %>% 
    select(-Key) %>% 
    summarise_all(list(
      m = median,
      l = \(x) quantile(x, 0.025),
      u = \(x) quantile(x, 0.975)
    ))
  
  pars_m$Epidemiology <- pars$Epidemiology %>% 
    group_by(Age) %>% 
    select(-Key) %>% 
    summarise_all(median) %>% 
    crossing(Key = 1:pars_m$N_Sims)
  
  pars_m$Epidemiology[pars_m$Epidemiology$Key == 4, "r_hz"] <- sts$r_hz_l
  pars_m$Epidemiology[pars_m$Epidemiology$Key == 5, "r_hz"] <- sts$r_hz_u
  
  pars_m$Epidemiology[pars_m$Epidemiology$Key == 6, "p_gp"] <- sts$p_gp_u
  pars_m$Epidemiology[pars_m$Epidemiology$Key == 7, "p_gp"] <- sts$p_gp_l
  
  pars_m$Epidemiology[pars_m$Epidemiology$Key == 8, "p_phn"] <- sts$p_phn_l
  pars_m$Epidemiology[pars_m$Epidemiology$Key == 9, "p_phn"] <- sts$p_phn_u
  
  pars_m$Epidemiology[pars_m$Epidemiology$Key == 10, "r_mor_hz"] <- sts$r_mor_hz_l
  pars_m$Epidemiology[pars_m$Epidemiology$Key == 11, "r_mor_hz"] <- sts$r_mor_hz_u
  
  
  sts <- pars$CostEff %>% 
    group_by(Age) %>% 
    select(-Key) %>% 
    summarise_all(list(
      m = median,
      l = \(x) quantile(x, 0.025),
      u = \(x) quantile(x, 0.975)
    ))
  
  pars_m$CostEff <- pars$CostEff %>% 
    group_by(Age) %>% 
    select(-Key) %>% 
    summarise_all(median) %>% 
    crossing(Key = 1:pars_m$N_Sims)
  
  pars_m$CostEff[pars_m$CostEff$Key == 12, "QL_ph"] <- sts$QL_ph_l
  pars_m$CostEff[pars_m$CostEff$Key == 13, "QL_ph"] <- sts$QL_ph_u
  pars_m$CostEff[pars_m$CostEff$Key == 12, "QL_pn"] <- sts$QL_pn_l
  pars_m$CostEff[pars_m$CostEff$Key == 13, "QL_pn"] <- sts$QL_pn_u
  pars_m$CostEff[pars_m$CostEff$Key == 12, "QL_0"] <- sts$QL_0_l
  pars_m$CostEff[pars_m$CostEff$Key == 13, "QL_0"] <- sts$QL_0_u
  
  pars_m$CostEff[pars_m$CostEff$Key == 14, "cost_GP_pp_non_PHN_HZ_inf"] <- sts$cost_GP_pp_non_PHN_HZ_inf_l
  pars_m$CostEff[pars_m$CostEff$Key == 15, "cost_GP_pp_non_PHN_HZ_inf"] <- sts$cost_GP_pp_non_PHN_HZ_inf_u
  
  pars_m$CostEff[pars_m$CostEff$Key == 16, "cost_GP_pp_PHN_inf"] <- sts$cost_GP_pp_PHN_inf_l
  pars_m$CostEff[pars_m$CostEff$Key == 17, "cost_GP_pp_PHN_inf"] <- sts$cost_GP_pp_PHN_inf_u
  
  pars_m$CostEff[pars_m$CostEff$Key == 18, "cost_hosp_pp_inf"] <- sts$cost_hosp_pp_inf_l
  pars_m$CostEff[pars_m$CostEff$Key == 19, "cost_hosp_pp_inf"] <- sts$cost_hosp_pp_inf_u
  
  
  pars_m$VE_ZVL <- pars$VE_ZVL %>% 
    group_by(Vaccine, Age, TimeVac) %>% 
    select(-Key) %>% 
    summarise_all(median) %>% 
    crossing(Key = 1:pars_m$N_Sims)
  
  
  for (k in c("VE_RZV_2d", "VE_RZV_1d", "VE_ReRZV_2d", "VE_ReRZV_1d")) {
    sts <- pars[[k]] %>% 
      group_by(Vaccine, TimeVac) %>% 
      select(-Key) %>% 
      summarise_all(list(
        l = \(x) quantile(x, 0.025),
        u = \(x) quantile(x, 0.975)
      ))
    
    pars_m[[k]] <- pars[[k]] %>% 
      group_by(Vaccine, TimeVac) %>% 
      select(-Key) %>% 
      summarise_all(median) %>% 
      crossing(Key = 1:pars_m$N_Sims)
    
    pars_m[[k]][pars_m[[k]]$Key == 2, "Protection"] <- sts$l
    pars_m[[k]][pars_m[[k]]$Key == 3, "Protection"] <- sts$u
  }
  
  
  keys <- 1:pars_m$N_Sims
  
  pb <- txtProgressBar(min = min(age0s), max = max(age0s), style = 3,  width = 50, char = "=") 
  
  yss <- age0s %>% 
    lapply(\(age0) {
      setTxtProgressBar(pb, age0)
      keys %>% a_run_uv(pars_m, age0 = age0)
    }) %>% 
    bind_rows() %>% 
    relocate(Key, Scenario, Age0, Arm) %>% 
    left_join(pars_m$key_mapping)
  
  
  stats_uni <- local({
    temp <- yss %>% 
      pivot_longer(-c(Scenario, Age0, Arm, Key, N0, Year0, Pars, Direction), names_to = "Index")
    
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
        Thres20 = (dQ_All_d * 2e4 - dC_Med_d) / dN_VacRZV_d,
        Thres30 = (dQ_All_d * 3e4 - dC_Med_d) / dN_VacRZV_d,
      )
    
  }) %>% 
    filter(Arm == "RZV_2d") %>% 
    ungroup() %>% 
    select(Age0, Arm, Pars, Direction, Thres20, Thres30)
  
  
  
  sens_uni <- stats_uni %>% filter(!is.na(Pars))  %>% 
    left_join(stats_uni %>% filter(is.na(Pars)) %>% select(Age0, Thres20_0 = Thres20, Thres30_0 = Thres30))
  
  sens_uni
  
}



summarise_sens_uni <- function(sens_uni, prefix = "", folder = NA, ext = ".pdf") {
  require(tidyverse)
  
  if (!is.na(folder)) {
    root_tab <- here::here("docs", "tabs", folder)
    root_fig <- here::here("docs", "figs", folder)
    dir.create(root_tab, showWarnings = F)
    dir.create(root_fig, showWarnings = F)
  } else {
    root_tab <- here::here("docs", "tabs")
    root_fig <- here::here("docs", "figs")
  }
  
  if (prefix != "") {
    prefix_fig <- glue::as_glue("g_") + prefix + "_"
    prefix_tab <- glue::as_glue(prefix) + "_"
  } else {
    prefix_fig <- glue::as_glue("g_")
    prefix_tab <- glue::as_glue("")
  }
  
  write_csv(sens_uni, here::here(root_tab, prefix_tab + "sens_ce_uni.csv"))
  
  
  labs_lvs <- c(
    "VE" = "Vaccine effectiveness",
    "Inc" = "Incidence, HZ",
    "PHN" = "PHN %",
    "Hosp" = "Hospitalisation %",
    "Mor" = "Mortality, HZ",
    "Q_HZ" = "QALY, HZ",
    "C_GP_PHN" = "Cost, GP, PHN",
    "C_GP_NonPHN" = "Cost, GP, No PHN",
    "Q_Life" = "QALY, Survival",
    "C_Hosp" = "Cost, Hospitalisation"    
  )
  
  age0s <- sens_uni %>% pull(Age0) %>% unique() %>% sort()
  
  limit <- sens_uni %>% pull(Thres20) %>% range()
  
  gs_sens <- lapply(age0s, \(age0) {
    ord <- sens_uni %>%
      filter(Age0  == age0) %>% 
      group_by(Pars) %>% 
      summarise(size = abs(diff(range(Thres20)))) %>% 
      arrange(size) %>% 
      pull(Pars)
    
    sens_uni %>%
      filter(Age0  == age0) %>% 
      mutate(
        Pars = factor(Pars, ord),
        i = as.numeric(Pars),
        D = Direction,
        a0 = paste0("Age of vaccination: ", Age0)
      ) %>% 
      ggplot(aes(y = Pars)) +
      geom_rect(aes(xmin = Thres20_0, xmax = Thres20, ymin = i - 0.45, ymax = i + 0.45,
                    fill = D)) +
      scale_fill_discrete("Changes") +
      scale_x_continuous("Threshold price, GBP per administration", limits = limit) +
      scale_y_discrete("Components", labels = labs_lvs) +
      facet_wrap(a0~.)
  })
  
  g_sens_ce <- ggpubr::ggarrange(plotlist = gs_sens, 
                                 ncol = 2, nrow = ceiling(length(age0s) / 2), 
                                 common.legend = T, legend = "bottom")
  
  ggsave(g_sens_ce, filename = here::here(root_fig, prefix_fig + "sens_ce_uni" + ext), width = 6, height = 10)
  
  
  
}