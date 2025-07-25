library(tidyverse)
library(readxl)
library(ggpubr)
library(tidybayes)

theme_set(theme_bw())


apply_lor <- function(p0, lor) 1 / (1 + exp(-log(p0 / (1 - p0)) - lor))


dat_ve <- read_xlsx(here::here("data", "processed_vaccine", "VE.xlsx"), sheet = 1) %>% 
  filter(Use) %>% 
  filter(Type == "HZ") %>% 
  mutate(
    M = M / 100,
    L = L / 100, 
    U = U / 100,
    Yr = gsub("yr", "", `Sub-group`),
    Yr = ifelse(is.na(Yr), 1.5, as.numeric(Yr))
  )



load(here::here("pars", "pars_ve_lor.rdata"))
load(here::here("pars", "pars_ve_zvl_rwa_zlg.rdata"))
load(here::here("pars", "pars_ve_rzv_rw_zlg.rdata"))



## ZVL VE -----
pars_ve_zvl %>% 
  mutate(AgeVac = Age - Yr) %>% 
  filter(AgeVac %in% seq(70, 75, 5)) %>% 
  filter(Yr < 20) %>% 
  #filter(Yr %in% c(5, 10)) %>% 
  group_by(IC, Vaccine, AgeVac, Yr) %>% 
  summarise(
    VE = mean(VE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = Yr, y = VE, colour = as.character(AgeVac))) +
  scale_y_continuous("Vaccine effectiveness, %", labels = scales::percent) +
  scale_x_continuous("Year since vaccination",) +
  expand_limits(y = 0:1, x = 0)



ve_zvl <- pars_ve_zvl %>% 
  mutate(
    AgeVac = Age - Yr,
    Tag = paste0("ZVL vaccinated at ", AgeVac)
  ) %>% 
  filter(AgeVac %in% seq(70, 75, 5)) %>% 
  filter(Yr <= 20) %>% 
  #filter(Yr %in% c(5, 10)) %>% 
  group_by(Tag, Yr) %>% 
  summarise(
    VE = mean(VE)
  )



g_zvl_gof <- pars_ve_zvl %>%
  filter(Age - Yr == 70) %>% 
  filter(Yr <= 20) %>% 
  group_by(Yr) %>% 
  summarise(
    M = mean(VE),
    L = quantile(VE, 0.025),
    U = quantile(VE, 0.975)
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Yr, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Yr, y = M)) +
  geom_pointrange(data = dat_ve %>% filter(Source == "Klein 2023"), aes(x = Yr, y = M, ymin = L, ymax = U)) +
  scale_y_continuous("Vaccine effectiveness, %", label = scales::percent) +
  scale_x_continuous("Year since vaccinated") +
  expand_limits(y = 0:1)


g_zvl_var <- ve_zvl %>% 
  ggplot() +
  geom_line(aes(x = Yr, y = VE, colour = Tag)) +
  geom_pointrange(data = dat_ve %>% filter(Source == "Klein 2023"), aes(x = Yr, y = M, ymin = L, ymax = U)) +
  scale_y_continuous("Vaccine effectiveness, %", labels = scales::percent) +
  scale_x_continuous("Year since vaccinated", breaks = c(1, seq(5, 20, 5))) +
  scale_colour_discrete("") +
  expand_limits(y = 0:1, x = 20) +
  theme(legend.position = c(0, 0), legend.justification = c(-0.1, -0.1))



g_zvl <- ggpubr::ggarrange(
  g_zvl_gof + labs(subtitle = "(A) Goodness of fit"), 
  g_zvl_var + labs(subtitle = "(B) Vaccine effectiveness by age of vaccination")
)


g_zvl


# Alternative plot
d_zvl2 <- pars_ve_zvl %>%
  filter(Age - Yr == 70) %>% 
  filter(Yr <= 20) %>% 
  group_by(Yr) %>% 
  summarise(
    M = mean(VE),
    L = quantile(VE, 0.025),
    U = quantile(VE, 0.975)
  )
g_zvl_gof_alt <- ggplot(data = d_zvl2) +
  geom_ribbon(aes(x = Yr, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Yr, y = M), linetype = "32") +
  geom_pointrange(data = dat_ve %>% filter(Source == "Klein 2023"), aes(x = Yr, y = M, ymin = L, ymax = U), fatten = 1) +
  scale_y_continuous("Vaccine effectiveness, %", label = scales::percent) +
  scale_x_continuous("Year since vaccinated") +
  expand_limits(y = 0:1)


g_zvl_var_alt <- ve_zvl %>% 
  ggplot() +
  geom_line(aes(x = Yr, y = VE, colour = Tag), linewidth = 1.0) +
  geom_line(data = d_zvl2, aes(x = Yr, y = M, colour = "Klein 2023"), linetype = "32") +
  scale_y_continuous("Vaccine effectiveness, %", labels = scales::percent) +
  scale_x_continuous("Year since vaccinated", breaks = c(1, seq(5, 20, 5))) +
  scale_colour_manual(NULL, values = c("black", "orange", "purple")) +
  expand_limits(y = 0:1, x = 20) +
  theme(legend.position = c(0.5, 0.7), legend.justification = c(-0.1, -0.1),
    legend.background = element_blank())


## RZV VE -----
d <- pars_ve_rzv%>% 
  filter(Yr <= 20) %>% 
  #filter(Yr %in% c(5, 10)) %>% 
  group_by(Vaccine, Yr) %>% 
  summarise(
    VE = median(VE)
  ) %>% 
  mutate(
    Tag = "Real-world, two doses (Baseline)"
  )


ve_rzv <- bind_rows(
  d, 
  d %>% 
    mutate(
      VE = apply_lor(VE, - lor_rw),
      Tag = "Trial, two doses"
    ),
  d %>% 
    mutate(
      VE = apply_lor(VE, lor_single),
      Tag = "Real-world, single dose"
    ),
  d %>%
    mutate(
      VE = apply_lor(VE, lor_re),
      Tag = "Real-world, two doses after ZVL"
    ),
  d %>%
    mutate(
      VE = apply_lor(VE, lor_re + lor_single),
      Tag = "Real-world, single dose after ZVL"
    )
) %>% 
  mutate(
    Tag = factor(Tag, c("Trial, two doses", 
                        "Real-world, two doses (Baseline)", 
                        "Real-world, two doses after ZVL", 
                        "Real-world, single dose",
                        "Real-world, single dose after ZVL"))
  )


g_rzv_gof <- pars_ve_rzv %>% 
  filter(Yr <= 20) %>% 
  group_by(Yr) %>% 
  summarise(
    M = mean(apply_lor(VE, - lor_rw)),
    L = quantile(apply_lor(VE, - lor_rw), 0.025),
    U = quantile(apply_lor(VE, - lor_rw), 0.975)
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Yr, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Yr, y = M)) +
  geom_pointrange(data = dat_ve %>% filter(!Realworld), aes(x = Yr, y = M, ymin = L, ymax = U)) +
  scale_y_continuous("Vaccine efficacy, %", label = scales::percent) +
  scale_x_continuous("Year since vaccinated") +
  expand_limits(y = 0)


g_rzv_var <- ve_rzv %>% 
  ggplot() +
  geom_line(aes(x = Yr, y = VE, colour = Tag)) +
  geom_pointrange(data = dat_ve %>% filter(!Realworld), aes(x = Yr, y = M, ymin = L, ymax = U)) +
  scale_y_continuous("Vaccine efficacy/effectiveness, %", labels = scales::percent) +
  scale_x_continuous("Year since vaccinated", breaks = c(1, seq(5, 20, 5))) +
  scale_colour_discrete("") +
  expand_limits(y = 0:1, x = 20) +
  theme(legend.position = c(0, 0), legend.justification = c(-0.1, -0.1))



g_rzv <- ggpubr::ggarrange(
  g_rzv_gof + labs(subtitle = "(A) Goodness of fit"), 
  g_rzv_var + labs(subtitle = "(B) Variants of implementation")
)

g_rzv

# Plot focusing on main VE estimate
d_rw <- pars_ve_rzv %>% 
  filter(Yr <= 20) %>% 
  group_by(Yr) %>% 
  summarise(
    M = mean(VE),
    L = quantile(VE, 0.025),
    U = quantile(VE, 0.975)
  )
d_rw$Tag <- "Real-world, two doses (Baseline)"

g_rzv_gof_alt <- pars_ve_rzv %>% 
  filter(Yr <= 20) %>% 
  group_by(Yr) %>% 
  summarise(
    M = mean(apply_lor(VE, - lor_rw)),
    L = quantile(apply_lor(VE, - lor_rw), 0.025),
    U = quantile(apply_lor(VE, - lor_rw), 0.975)
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Yr, ymin = L, ymax = U, fill = "Trial, two doses"), alpha = 0.2) +
  geom_line(aes(x = Yr, y = M, colour = "Trial, two doses", linetype = "Trial, two doses")) +
  geom_pointrange(data = dat_ve %>% filter(!Realworld), aes(x = Yr, y = M, ymin = L, ymax = U), fatten = 1) +
  scale_y_continuous("Vaccine efficacy, %", label = scales::percent) +
  scale_x_continuous("Year since vaccinated") +
  scale_colour_manual(NULL, values = c("black"), aesthetics = c("fill", "colour")) +
  scale_linetype_manual(NULL, values = c("32")) +
  expand_limits(y = 0) +
  theme(legend.position = c(0, 0), legend.justification = c(-0.1, -0.1))

g_rzv_var_alt <- ggplot(ve_rzv) +
  geom_ribbon(data = d_rw, aes(x = Yr, ymin = L, ymax = U, fill = Tag), alpha = 0.2) +
  geom_line(aes(x = Yr, y = VE, colour = Tag, linetype = Tag)) +
  scale_y_continuous("Vaccine efficacy/effectiveness, %", labels = scales::percent) +
  scale_x_continuous("Year since vaccinated", breaks = c(1, seq(5, 20, 5))) +
  scale_colour_manual(NULL, values = c("black", "#36f", "#36f", "#e3a", "#e3a"), aesthetics = c("fill", "colour")) +
  scale_linetype_manual(NULL, values = c("32", "solid", "32", "solid", "32")) +
  expand_limits(y = 0:1, x = 20) +
  theme(legend.position = c(0, 0), legend.justification = c(-0.01, -0.1), legend.background = element_blank())

g_alt <- ggpubr::ggarrange(
  g_rzv_gof_alt + labs(subtitle = "(A) RZV: Goodness of fit"), 
  g_rzv_var_alt + labs(subtitle = "(B) RZV: Variants of implementation"),
  g_zvl_gof_alt + labs(subtitle = "(C) ZVL: Goodness of fit"), 
  g_zvl_var_alt + labs(subtitle = "(D) ZVL: Effectiveness by age"),
  nrow = 2, ncol = 2
)

# zvl_prev <- read_csv(here::here("data", "previous", "2019_fig9_ZVL.csv"))
# 
# g_vax_comp_zvl <- ggplot(data = d_zvl2[d_zvl2$Yr <= 10, ]) +
#   geom_ribbon(aes(x = Yr, ymin = L, ymax = U), fill = "#36f", alpha = 0.2) +
#   geom_line(aes(x = Yr, y = M), colour = "#36f") +
#   geom_line(data = zvl_prev[zvl_prev$time >= 1, ], aes(x = time, y = Age60), col = "#db2", linewidth = 2) +
#   geom_line(data = zvl_prev[zvl_prev$time >= 1, ], aes(x = time, y = Age70), col = "#ec3", linewidth = 2) +
#   geom_line(data = zvl_prev[zvl_prev$time >= 1, ], aes(x = time, y = Age80), col = "#fd4", linewidth = 2) +
#   scale_y_continuous("VE for ZVL, %", label = scales::percent) +
#   scale_x_continuous("Year since vaccinated", breaks = c(2, 4, 6, 8, 10)) +
#   expand_limits(y = 0:1)
# 
# rzv_prev <- read_csv(here::here("data", "previous", "2019_fig10_RZV.csv"))
# 
# g_vax_comp_rzv <- ggplot(ve_rzv[ve_rzv$Tag == "Real-world, two doses (Baseline)" & ve_rzv$Yr <= 10, ]) +
#   geom_ribbon(data = d_rw[d_rw$Yr <= 10, ], aes(x = Yr, ymin = L, ymax = U), fill = "#36f", alpha = 0.2) +
#   geom_line(aes(x = Yr, y = VE), colour = "#36f") +
#   geom_line(data = rzv_prev[rzv_prev$time >= 1, ], aes(x = time, y = Age70), col = "#ec3", linewidth = 2) +
#   scale_y_continuous("VE for RZV, %", labels = scales::percent) +
#   scale_x_continuous("Year since vaccinated", breaks = c(2, 4, 6, 8, 10)) +
#   expand_limits(y = 0:1, x = 10) +
#   theme(legend.position = c(0, 0), legend.justification = c(-0.01, -0.1), legend.background = element_blank())
# 
# g_vax_comp <- ggpubr::ggarrange(
#   g_vax_comp_zvl + labs(subtitle = "(A) ZVL: Comparison"), 
#   g_vax_comp_rzv + labs(subtitle = "(B) RZV: Comparison"),
#   nrow = 1, ncol = 2
# )



## ZVL uptake


load(here::here("data", "processed_vaccine", "coverage.rdata"))
load(here::here("pars", "fitted_coverage.rdata"))


d <- coverage %>% 
  mutate(
    Age = Age - 1,
    Year = as.numeric(Year),
    Cohort = Year - Age + 70
  ) %>% 
  filter(Cohort >= 2014)



g_uptake_gof <- pred1$pred %>% 
  filter(Cohort == "2014") %>% 
  ggplot(aes(x = Age - 1)) +
  geom_ribbon(aes(ymin = Coverage_l, ymax = Coverage_u), alpha = 0.2) +
  geom_line(aes(y = Coverage)) +
  geom_point(data = d, aes(x = Age, y = value, colour = as.character(Cohort)), size = rel(2)) +
  geom_line(data = d, aes(x = Age, y = value, colour = as.character(Cohort))) +
  scale_y_continuous("Coverage, %", labels = scales::percent) +
  scale_colour_discrete("Cohort (year of 70 YOA)") +
  scale_x_continuous("Age", breaks = seq(70, 80, 2)) +
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(1, 0), legend.justification = c(1.1, -0.1))



g_uptake_gof


dir.create("docs/figs/inputs", showWarnings = T)
ggsave(g_zvl, filename = here::here("docs", "figs", "inputs", "g_vaccine_ve_zvl.png"), width = 12, height = 6)
ggsave(g_rzv, filename = here::here("docs", "figs", "inputs", "g_vaccine_ve_rzv.png"), width = 12, height = 6)
ggsave(g_alt, filename = here::here("docs", "figs", "inputs", "g_vaccine_ve_alt.png"), width = 10, height = 9)
ggsave(g_uptake_gof, filename = here::here("docs", "figs", "inputs", "g_vaccine_uptake.png"), width = 6, height = 4.5)
#ggsave(g_vax_comp, filename = here::here("docs", "figs", "g_vaccine_comparison.png"), width = 7, height = 3)


