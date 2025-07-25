

vis_thres <- function(stats_uv, stats_re) {
  require(ggplot2)

  gs <- list()

  d_main <- bind_rows(
    stats_uv$stats_uv_ce %>% 
      filter(Index %in% c("Thres20_50", "Thres30_90")) %>% 
      filter(Age0 >= 60) %>% 
      filter((Arm == "RZV_1d" & Age0 >= 80) | Arm == "RZV_2d") %>% 
      mutate(
        Arm = case_when(
          Age0 < 80 & Arm == "RZV_2d" ~ "Target",
          Arm == "RZV_2d" ~ Arm,
          Arm == "RZV_1d" ~ Arm,
          T ~ NA
        )
      ) %>% 
      filter(!is.na(Arm)) %>% 
      select(Age = Age0, Arm, Index, M),
    stats_re$stats_re_ce %>% 
      filter(Scenario != "Overall") %>% 
      filter(Index %in% c("Thres20_50", "Thres30_90")) %>% 
      filter(Age1 >= 80) %>% 
      filter(Age0 %in% c(70, 75)) %>%
      mutate(Arm = paste(Arm, Age0)) %>% 
      select(Age = Age1, Arm, Index, M)
  )
  
  bound <- d_main %>% 
      filter(Arm == "Target") %>% 
      pivot_wider(names_from = Index, values_from = M) %>% 
      mutate(thres = pmin(Thres20_50, Thres30_90)) %>% 
      pull(thres) %>% range()
    
  tags <- c(
    "Target" = "**A** Two-dose RZV<br>Existing programme",
    "RZV_2d" = "**B** Two-dose RZV<br>",
    "RZV_1d" = "**C** Single-dose RZV<br>",
    "ReRZV_2d 70" = "**D** Two-dose RZV<br>ZVL at age 70",
    "ReRZV_1d 70" = "**E** Single-dose RZV<br>ZVL at age 70"
  )
  
  gs$g_panel <- d_main %>% 
    filter(Age <= 95) %>% 
    filter(Arm %in% names(tags)) %>%
    mutate(Arm = factor(Arm, names(tags))) %>% 
    ggplot() +
    geom_line(aes(x = Age, y = M, colour = Index)) +
    geom_hline(yintercept = bound, linetype = 2) +
    scale_y_continuous("Threshold price per adminstration, GBP") +
    scale_color_discrete("Threshold", labels = c(Thres20_50 = "50% CE given 20,000 WTP", Thres30_90 = "90% CE given 30,000 WTP")) +
    facet_grid(.~Arm, scales = "free_x", labeller = labeller(Arm = tags)) +
    expand_limits(y = c(0, 200)) +
    theme(legend.position = "bottom", 
      strip.background = element_blank(), strip.text = ggtext::element_markdown(size = 11))
  
  
  ## Composition
  labs_comp <- c(
    "dQ_HZ_d" = "HZ prevention",
    "dQ_Life_d" = "Survival", 
    "dC_GP_d" = "GP care",
    "dC_Hosp_d" = "Hospitalization",
    "dC_VacRZV_d" = "Vaccination"
  )
  
  labs_arm <- c(RZV_2d = "Two doses", RZV_1d = "Single dose")
  
  waterfall <- stats_uv$stats_uv_ce %>% 
    select(Scenario, Age0, Arm, Index, M) %>% 
    filter(Index %in% c(names(labs_comp), "dN_VacRZV_d", "dQ_HZ_d")) %>% 
    pivot_wider(names_from = Index, values_from = M) %>%
    mutate(
      #dC_GP_d = dC_GP_NonPHN_d + dC_GP_PHN_d,
      Thres = (2e4 * (dQ_HZ_d + dQ_Life_d) - dC_Hosp_d - dC_GP_d) / dN_VacRZV_d,
      dC_VacRZV_d = Thres * dN_VacRZV_d
    ) %>% 
    select(Age0, Arm, dQ_HZ_d, dQ_Life_d, dC_Hosp_d, dC_GP_d, dC_VacRZV_d) %>% 
    pivot_longer(-c(Age0, Arm), names_to = "Index") %>% 
    mutate(
      Arm = factor(Arm, c("RZV_1d", "RZV_2d")),
      Type = ifelse(startsWith(Index, "dQ"), "QoL", "Cost"),
      Index = factor(Index, names(labs_comp)),
      id = as.numeric(Index),
      MB = ifelse(Type == "QoL", value * 2e4, -value)
    )
  
  peak <- waterfall %>% 
    filter(Index %in% c("dC_GP_d", "dC_Hosp_d")) %>% 
    group_by(Arm, Age0) %>% 
    mutate(cMB = c(0, cumsum(MB)[-2])) %>% 
    group_by(Arm, Index) %>% 
    filter(MB == max(MB))
  
  gs$g_wf <- waterfall %>% 
    group_by(Age0, Arm) %>% 
    arrange(Index) %>% 
    mutate(
      y1 = cumsum(MB),
      y0 = c(0, y1[-n()]),
      Direction = ifelse(MB > 0, "Beneficial", "Harmful")
    ) %>% 
    filter(Age0 %in% c(70, 80)) %>% 
    ggplot() +
    geom_rect(aes(x = Index, xmin = id - 0.4, xmax = id + 0.4, ymin = y0, ymax = y1, fill = Direction), alpha = 0.5) +
    geom_segment(aes(x=ifelse(id == last(id), last(id) - 0.5, id - 0.4), 
                     xend=ifelse(id == last(id), last(id) + 0.5, id + 1.4), 
                     y=y1, 
                     yend=y1)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 2.5) +
    annotate("text", x = 2.3, y = 10, label = "QoL", hjust = 1) +
    annotate("text", x = 2.7, y = 10, label = "Cost", hjust = 0) +
    scale_y_continuous("Net monetary benefit (WTP = 20,000 GBP)") + 
    scale_x_discrete("Source", labels = labs_comp) +
    facet_grid(Age0~Arm, labeller = labeller(
      Arm = c(Vac_2d = "Two doses", Vac_1d = "Single dose"),
      Age0 = c("70" = "Vaccinated at 70 YOA", "80" = "Vaccinated at 80 YOA")
    )) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  
  gs$g_comp_stack <- waterfall %>% 
    group_by(Age0, Arm) %>% 
    filter(Index != "dC_VacRZV_d") %>% 
    #mutate(Index = fct_rev(Index)) %>% 
    filter(Age0 <= 95) %>% 
    filter(Age0 >= 60) %>% 
    ggplot(aes(x = Age0)) +
    geom_bar(aes(y = MB, fill = Index), stat = "identity", position = "stack", colour = NA, width = 1) +
    scale_y_continuous("Monetary benefit in GBP") +
    scale_fill_discrete("Components", labels = labs_comp) +
    facet_grid(.~Arm, labeller = labeller(Arm = labs_arm)) +
    scale_x_continuous("Age at vaccination") +
    labs(caption = "1 QALY = 20,000 GBP") +
    theme(strip.background = element_blank(), strip.text = element_text(size = 11))

  
  gs$g_comp_stack_ann <- waterfall %>% 
    group_by(Age0, Arm) %>% 
    filter(Index != "dC_VacRZV_d") %>% 
    #mutate(Index = fct_rev(Index)) %>% 
    filter(Age0 <= 95) %>% 
    filter(Age0 >= 60) %>% 
    ggplot(aes(x = Age0)) +
    geom_bar(aes(y = MB, fill = Index), stat = "identity", position = "stack", colour = NA, width = 1, alpha = 0.7) +
    geom_linerange(data = peak, aes(x = Age0, ymin = cMB, ymax = cMB + MB), colour = "black") +
    geom_text(data = peak, aes(x = Age0, y = cMB + MB, 
      label = sprintf("£%s", round(MB, 1))), hjust = 0.5, vjust = -0.5, colour = "black") +
    scale_y_continuous("Monetary benefit in GBP") +
    scale_fill_discrete("Components", labels = labs_comp) +
    facet_grid(.~Arm, labeller = labeller(Arm = labs_arm)) +
    scale_x_continuous("Age at vaccination") +
    labs(caption = "1 QALY = 20,000 GBP") +
    theme(strip.background = element_blank(), strip.text = element_text(size = 11))
  
  
  gs$g_comp_stack_ann_al <- waterfall %>% 
    group_by(Age0, Arm) %>% 
    filter(Index != "dC_VacRZV_d") %>% 
    #mutate(Index = fct_rev(Index)) %>% 
    filter(Age0 <= 95) %>% 
    filter(Age0 >= 60) %>% 
    ggplot(aes(x = Age0)) +
    geom_bar(aes(y = MB, fill = Index, alpha = ifelse(Type == "Cost", 1, 0.9)), stat = "identity", position = "stack", colour = NA, width = 1) +
    geom_linerange(data = peak, aes(x = Age0, ymin = cMB, ymax = cMB + MB)) +
    geom_text(data = peak, aes(x = Age0, y = cMB + MB, label = sprintf("£%s", round(MB, 1))), hjust = 0.5, vjust = -0.5) +
    scale_y_continuous("Monetary benefit in GBP") +
    scale_fill_discrete("Components", labels = labs_comp) +
    scale_alpha(range = c(0.3, 1), guide = 'none') +
    facet_grid(.~Arm, labeller = labeller(Arm = labs_arm)) +
    scale_x_continuous("Age at vaccination") +
    labs(caption = "1 QALY = 20,000 GBP") +
    theme(strip.background = element_blank(), strip.text = element_text(size = 11))
  
  
  gs$g_comp_fill <- waterfall %>% 
    group_by(Age0, Arm) %>% 
    filter(Index != "dC_VacRZV_d") %>% 
    #mutate(Index = fct_rev(Index)) %>% 
    filter(Age0 <= 95) %>% 
    filter(Age0 >= 60) %>% 
    ggplot(aes(x = Age0)) +
    geom_bar(aes(y = MB, fill = Index), stat = "identity", position = "fill", colour = NA, width = 1) +
    scale_fill_discrete("Components", labels = labs_comp) +
    scale_y_continuous("Share of monetary benefit, %", labels = scales::percent) +
    facet_grid(.~Arm, labeller = labeller(Arm = labs_arm)) +
    scale_x_continuous("Age at vaccination") +
    labs(caption = "1 QALY = 20,000 GBP") +
    theme(strip.background = element_blank(), strip.text = element_text(size = 11))
  
  return(gs)
}


save_fig_thres <- function(gs, prefix = "", folder = NA, ext = ".pdf") {
  require(tidyverse)
  
  if (!is.na(folder)) {
    root <- here::here("docs", "figs", folder)
    dir.create(root, showWarnings = F)
  } else {
    root <- here::here("docs", "figs")
  }
  
  if (prefix != "") {
    prefix <- glue::as_glue("g_") + prefix + "_"
  } else {
    prefix <- glue::as_glue("g_")
  }
  
  ggsave(gs$g_panel, filename = here::here(root, prefix + "thres_panel" + ext), width = 8, height = 3.5)
  
  ggsave(gs$g_wf, filename = here::here(root, prefix + "rzv_waterfall" + ext), width = 8, height = 7.5)
  ggsave(gs$g_comp_stack, filename = here::here(root, prefix + "rzv_mb_stack" + ext), width = 8, height = 5.5)
  ggsave(gs$g_comp_stack_ann, filename = here::here(root, prefix + "rzv_mb_stack_ann" + ext), width = 8, height = 5.5)
  ggsave(gs$g_comp_stack_ann_al, filename = here::here(root, prefix + "rzv_mb_stack_ann_alp" + ext), width = 8, height = 5.5)
  ggsave(gs$g_comp_fill, filename = here::here(root, prefix + "rzv_mb_fill" + ext), width = 8, height = 5.5)
  
}



