
load_inputs_rzv <- function(pars_ce, f_ve_zvl, f_ve_rzv, seed = 11667) {
  vtype <- "rw"
  
  pars <- list()
  
  f_rzv <- here::here("pars", "pars_ve_rzv_uv2_rw_y10_zie.rdata")
  pars$y10_zie <- load_inputs(pars_ce, vtype, f_ve_zvl = f_ve_zvl, f_ve_rzv = f_rzv, seed = seed)
  
  f_rzv <- here::here("pars", "pars_ve_rzv_uv2_rw_y11_zie.rdata")
  pars$y11_zie <- load_inputs(pars_ce, vtype, f_ve_zvl = f_ve_zvl, f_ve_rzv = f_rzv, seed = seed)
  
  f_rzv <- here::here("pars", "pars_ve_rzv_uv2_rw_y11m_zie.rdata")
  pars$y11m_zie <- load_inputs(pars_ce, vtype, f_ve_zvl = f_ve_zvl, f_ve_rzv = f_rzv, seed = seed)
  
  
  f_rzv <- here::here("pars", "pars_ve_rzv_uv2_rw_y10_zig.rdata")
  pars$y10_zig <- load_inputs(pars_ce, vtype, f_ve_zvl = f_ve_zvl, f_ve_rzv = f_rzv, seed = seed)
  
  f_rzv <- here::here("pars", "pars_ve_rzv_uv2_rw_y11_zig.rdata")
  pars$y11_zig <- load_inputs(pars_ce, vtype, f_ve_zvl = f_ve_zvl, f_ve_rzv = f_rzv, seed = seed)
  
  f_rzv <- here::here("pars", "pars_ve_rzv_uv2_rw_y11m_zig.rdata")
  pars$y11m_zig <- load_inputs(pars_ce, vtype, f_ve_zvl = f_ve_zvl, f_ve_rzv = f_rzv, seed = seed)
  
  return(pars)
}


sens_rzv <- function(pars_rzv, age0s = 60:95) {
  
  stats <- lapply(pars_rzv, \(pars) {
    exec_cohort_rzv(pars, age0s = age0s) %>% 
      summarise_cohort()
  })
  
  fns <- c("y10_zie", "y11_zie", "y11m_zie", "y10_zig", "y11_zig", "y11m_zig")
  n_fns <- length(fns)
  
  
  ss <- list()
  
  ss$ves_2d <- bind_rows(lapply(1:n_fns, \(i) {
    pars_rzv[[i]]$VE_RZV_2d %>% mutate(fn = fns[i])
  })) %>% 
    filter(TimeVac <= 20) %>% 
    group_by(TimeVac, fn) %>% 
    summarise(
      M = median(Protection),
      L = quantile(Protection, 0.025),
      U = quantile(Protection, 0.975)
    )
  
  ss$ves_1d <- bind_rows(lapply(1:n_fns, \(i) {
    pars_rzv[[i]]$VE_RZV_1d %>% mutate(fn = fns[i])
  })) %>% 
    filter(TimeVac <= 20) %>% 
    group_by(TimeVac, fn) %>% 
    summarise(
      M = median(Protection),
      L = quantile(Protection, 0.025),
      U = quantile(Protection, 0.975)
    )
  
  ss$ces <- bind_rows(lapply(1:n_fns, \(i) {
    stats[[i]]$stats_uv_ce %>% mutate(fn = fns[i])
  }))
  
  ss$yss <- bind_rows(lapply(1:n_fns, \(i) {
    stats[[i]]$stats_uv_ys %>% mutate(fn = fns[i])
  }))
  
  return(ss)
}


summarise_sens_rzv <- function(sens_rzv, prefix = "", folder = NA, ext = ".pdf") {
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
    prefix_tab <- glue::as_glue(prefix) + "_"
    prefix_fig <- glue::as_glue("g_") + prefix + "_"
  } else {
    prefix_tab <- glue::as_glue("")
    prefix_fig <- glue::as_glue("g_")
  }
  
  
  labs_fn <- c(
    y10_zie = "ZI Exponential\nUp to 10th year",
    y11_zie = "ZI Exponential\nUp to 11th year",
    y11m_zie = "ZI Exponential\nUp to 11th year exc. 9th, 10th",
    y10_zig = "ZI Gamma\nUp to 10th year",
    y11_zig = "ZI Gamma\nUp to 11th year",
    y11m_zig = "ZI Gamma\nUp to 11th year exc. 9th, 10th"
  )
  
  ces <- sens_rzv$ces %>% 
    mutate(
      fn = factor(fn, levels = names(labs_fn))
    ) %>% 
    filter(Index %in% c("Thres20_50", "Thres30_90")) %>% 
    filter(Age0 >= 60) %>% 
    filter(Age0 <= 95) %>% 
    filter((Arm == "RZV_1d" & Age0 >= 80) | Arm == "RZV_2d") %>% 
    filter(!is.na(Arm)) %>% 
    select(Age = Age0, Arm, Index, M, fn)
  
  write_csv(ces, here::here(root_tab, prefix_tab + "sens_rzv.csv"))
  
  
  gs <- list()
  
  gs$g_sens_waning_ves <- sens_rzv$ves_2d%>% 
    mutate(
      fn = factor(fn, levels = names(labs_fn))
    ) %>% 
    ggplot() +
    geom_ribbon(aes(x = TimeVac, ymin = L, ymax = U), alpha = 0.3) +
    geom_line(aes(x = TimeVac, y = M)) +
    scale_y_continuous("Vaccine effectiveness, %, two doses", limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous("Year since vaccination") +
    guides(fill = guide_none(), colour = guide_none()) +
    facet_wrap(.~fn, labeller = labeller(fn = labs_fn), nrow = 1)
  
  gs$g_sens_waning_thres <- ces %>% 
    mutate(Arm = factor(Arm, c("RZV_2d", "RZV_1d"))) %>% 
    ggplot() +
    geom_line(aes(x = Age, y = M, colour = Index, linetype = Arm)) +
    scale_y_continuous("Threshold price per adminstration, GBP") +
    scale_color_discrete("Threshold", labels = c(Thres20_50 = "50% CE given 20,000 WTP", Thres30_90 = "90% CE given 30,000 WTP")) +
    scale_linetype_discrete("", labels = c(RZV_1d = "single dose", RZV_2d = "two doses")) +
    expand_limits(y = c(0, 200)) +
    facet_wrap(.~fn, labeller = labeller(fn = labs_fn), nrow = 1) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 60, hjust = 1))
  
  
  gs$g_sens_waning_bind <- ggpubr::ggarrange(
    gs$g_sens_waning_ves + labs(subtitle = "A"), 
    gs$g_sens_waning_thres + labs(subtitle = "B"), 
    nrow = 2)
  
  
  ggsave(gs$g_sens_waning_ves, filename = here::here(root_fig, prefix_fig + "sens_rzv_ves" + ext), width = 14, height = 4.5)
  ggsave(gs$g_sens_waning_thres, filename = here::here(root_fig, prefix_fig + "sens_rzv_thres" + ext), width = 14, height = 4.5)
  ggsave(gs$g_sens_waning_bind, filename = here::here(root_fig, prefix_fig + "sens_rzv_bind" + ext), width = 14, height = 8)
  
  return(gs)
}





