### Not run without 'data/raw/'

library(targets)
library(tidyverse)
library(tidybayes)
library(readxl)
library(ggpubr)

theme_set(theme_bw(base_size = 14))


pars <- local({load(here::here("pars", "pars_base_35_rw.rdata")); pars})

prefix <- names(pars)[1]
prefix <- gsub("Year0", "", prefix)

names(pars) <- gsub(prefix, "", names(pars))

pars$Epidemiology
load(file = here::here("data", "raw", "Epi_HZ_CPRD_23Nov19.rdata"))
epi_hz <- epi_hz %>% filter(IC == 0)
epi_gp <- epi_gp %>% filter(IC == 0)
epi_phn <- epi_phn %>% filter(IC == 0)

load(here::here("data", "processed_demography", "Population_ONS.rdata"))


gs <- list()
isize <- 2.8 # adjust width of intervals for plotting
alpha_unused <- 0.4


#### Epidemiology -----

gs$g_r_hz <- pars$Epidemiology %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = r_hz), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  geom_pointinterval(data = epi_hz %>% filter(Age < 90), aes(x = Age, y = M, ymin = L, ymax = U)) +
  geom_pointinterval(data = epi_hz %>% filter(Age >= 90), aes(x = Age, y = M, ymin = L, ymax = U), alpha = alpha_unused, linetype = 2) +
  scale_y_continuous("Cases per 1,000 person-years", labels = scales::number_format(scale = 1e3),
    breaks = seq(0, 12, 2) * 0.001) +
  coord_cartesian(ylim = c(0, 12e-3), xlim = c(48, 102), expand = FALSE) +
  labs(subtitle = "Herpes zoster incidence, total", colour = "Level")

gs$g_r_hz

gs$g_r_hz_gp <- pars$Epidemiology %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = p_gp * r_hz), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  geom_pointinterval(data = epi_gp %>% filter(Age < 90), aes(x = Age, y = M, ymin = L, ymax = U)) +
  geom_pointinterval(data = epi_gp %>% filter(Age >= 90), aes(x = Age, y = M, ymin = L, ymax = U), alpha = alpha_unused, linetype = 2) +
  scale_y_continuous("Cases per 1,000 person-years", labels = scales::number_format(scale = 1e3),
    breaks = seq(0, 12, 2) * 0.001) +
  coord_cartesian(ylim = c(0, 12e-3), xlim = c(48, 102), expand = FALSE) +
  labs(subtitle = "Herpes zoster incidence, general practice only", colour = "Level")

gs$g_r_hz_gp


gs$g_p_hosp <- pars$Epidemiology %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = 1 - p_gp), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  scale_y_continuous("Proportion of cases", labels = scales::percent) +
  coord_cartesian(ylim = c(0, 0.2), xlim = c(48, 102), expand = FALSE) +
  labs(subtitle = "Herpes zoster hospitalisation", colour = "Level")

gs$g_p_hosp


gs$g_r_hz_hosp <- pars$Epidemiology %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = (1 - p_gp) * r_hz), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  scale_colour_brewer() +
  scale_y_continuous("per 1,000 person-year", labels = scales::number_format(scale = 1e3)) +
  labs(subtitle = "Epidemiology: Incidence HZ, hospitalised episode", colour = "Level") +
  expand_limits(y = 0)

gs$g_r_hz_hosp


gs$g_p_phn <- pars$Epidemiology %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = p_phn), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  geom_pointinterval(data = epi_phn %>% filter(Age < 90), aes(x = Age, y = M, ymin = L, ymax = U)) +
  geom_pointinterval(data = epi_phn %>% filter(Age >= 90), aes(x = Age, y = M, ymin = L, ymax = U), alpha = alpha_unused, linetype = 2) +
  scale_y_continuous("Proportion of cases", labels = scales::percent) +
  labs(subtitle = "Postherpetic neuralgia", colour = "Level") +
  coord_cartesian(ylim = c(0, 0.4), xlim = c(48, 102), expand = FALSE)

gs$g_p_phn


gs$g_p_mor <- pars$Epidemiology %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = r_mor_hz), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  scale_y_continuous("Deaths per 1,000,000 person-years", labels = scales::number_format(scale = 1e6)) +
  labs(subtitle = "Herpes zoster mortality", colour = "Level") +
  coord_cartesian(ylim = c(0, 3.3e-5), xlim = c(48, 102), expand = FALSE)

gs$g_p_mor


# gs$g_r_hz_mor <- pars$Epidemiology %>% 
#   ggplot() +
#   stat_interval(aes(x = Age, y = r_mor_hz * r_hz), linewidth = isize) +
#   scale_colour_brewer() +
#   scale_y_continuous("Mortality HZ, per 1,000", labels = scales::number_format(scale = 1e3)) +
#   labs(colour = "Level") +
#   expand_limits(y = 0)
# 
# gs$g_r_hz_mor


#### QoL -----

gs$g_ql_ph <- pars$CostEff %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = QL_ph), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  scale_y_continuous("QALYs per case") +
  labs(subtitle = "Health-related QoL loss due to herpes zoster", colour = "Level") +
  coord_cartesian(ylim = c(0, 0.22), xlim = c(48, 102), expand = FALSE)

gs$g_ql_ph


gs$g_ql_pn <- pars$CostEff %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = QL_pn), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  scale_y_continuous("QALYs per case") +
  labs(subtitle = "Health-related QoL loss due to HZ", colour = "Level") +
  coord_cartesian(ylim = c(0, 0.065), xlim = c(48, 102), expand = FALSE)

gs$g_ql_pn


gs$g_ql_0 <- pars$CostEff %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = QL_0), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  scale_y_continuous("QALYs per case") +
  labs(subtitle = "Health-related QoL loss due to HZ", colour = "Level") +
  coord_cartesian(ylim = c(0, 0.1501), xlim = c(48, 102), expand = FALSE)

gs$g_ql_0


#### Demography -----

demo <- demo_ons %>% filter(Location == "England") %>% 
  filter(Age < 100) %>% 
  mutate(
    AgeGrp = cut(Age, seq(0, 100, 10), right = F)
  ) %>% 
  group_by(AgeGrp, Year) %>% 
  summarise(
    Deaths = sum(N * mortality),
    N = sum(N),
    Dr = Deaths / N
  ) %>% 
  filter(Year %in% c(2023, 2028, 2033, 2045)) 

gs$g_demo_pop <- demo %>% 
  ggplot() +
  geom_bar(aes(x = N, y = AgeGrp), stat = "identity") +
  scale_x_continuous("Population size, million", labels = scales::number_format(scale = 1e-6)) +
  scale_y_discrete("Age group") +
  facet_grid(Year~.)


gs$g_demo_deaths <- demo %>% 
  ggplot() +
  geom_bar(aes(x = Deaths, y = AgeGrp), stat = "identity") +
  scale_x_continuous("Deaths, thousand", labels = scales::number_format(scale = 1e-3)) +
  scale_y_discrete("Age group") +
  facet_grid(Year~.)


gs$g_demo_dr <- demo %>% 
  ggplot() +
  geom_bar(aes(x = Dr, y = AgeGrp), stat = "identity") +
  scale_x_continuous("Background mortality, %", labels = scales::percent) +
  scale_y_discrete("Age group") +
  facet_grid(Year~.)


demo_ons %>% filter(Location == "England") %>% 
  filter(Age < 100 & Age >= 50) %>% 
  filter(Year == 2023) %>% 
  ggplot() +
  geom_line(aes(x = Age, y = mortality))


gs$g_demo <- ggarrange(
  gs$g_demo_pop + labs(subtitle = "(A)"), 
  gs$g_demo_dr + labs(subtitle = "(B)"), 
  nrow = 1, align = "v"
)

#### Cost -----

gs$g_cost_hosp <- pars$CostEff %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = cost_hosp_pp_inf), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  scale_y_continuous("GBP (£)") +
  coord_cartesian(xlim = c(48, 102), ylim = c(0, 5600), expand = FALSE) +
  labs(subtitle = "Direct medical cost, hospitalised episode", colour = "Level")


iqr <- pars$CostEff %>% 
  filter(Age == 70) %>% 
  pull(cost_GP_pp_non_PHN_HZ_inf) %>% 
  quantile(c(0.1, 0.5, 0.9))

gs$g_cost_gp <- pars$CostEff %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = cost_GP_pp_non_PHN_HZ_inf), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  scale_y_continuous("GBP (£)") +
  coord_cartesian(xlim = c(48, 102), iqr[2] + 8 * c(iqr[1] - iqr[2], iqr[3] - iqr[2]), expand = FALSE) +
  labs(subtitle = "Direct medical cost, GP episode without PHN", colour = "Level")


iqr <- pars$CostEff %>% 
  filter(Age == 70) %>% 
  pull(cost_GP_pp_PHN_inf) %>% 
  quantile(c(0.1, 0.5, 0.9))

gs$g_cost_gpphn <- pars$CostEff %>% 
  ggplot() +
  stat_lineribbon(aes(x = Age, y = cost_GP_pp_PHN_inf), .width = c(.95, .8, .5), color = "#08519C", alpha = 0.6) +
  scale_fill_brewer("Interval") +
  scale_y_continuous("GBP (£)") +
  coord_cartesian(xlim = c(48, 102), iqr[2] + 8 * c(iqr[1] - iqr[2], iqr[3] - iqr[2]), expand = FALSE) +
  labs(subtitle = "Direct medical cost, GP episode with PHN", colour = "Level")


#### Vaccine
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


## ZVL VE -----
# load(here::here("pars", "pars_ve_zvl_rwa.rdata"))
# pars_ve_zvl %>% 
#   mutate(AgeVac = Age - Yr) %>% 
#   filter(AgeVac %in% seq(70, 75, 5)) %>% 
#   filter(Yr < 20) %>% 
#   #filter(Yr %in% c(5, 10)) %>% 
#   group_by(IC, Vaccine, AgeVac, Yr) %>% 
#   summarise(
#     VE = mean(VE)
#   ) %>% 
#   ggplot() +
#   geom_line(aes(x = Yr, y = VE, colour = as.character(AgeVac))) +
#   scale_y_continuous("Vaccine effectiveness, %", labels = scales::percent) +
#   scale_x_continuous("Year since vaccination",) +
#   expand_limits(y = 0:1, x = 0)


# ve_zvl <- pars_ve_zvl %>% 
#   mutate(
#     AgeVac = Age - Yr,
#     Tag = paste0("ZVL vaccinated at ", AgeVac)
#   ) %>% 
#   filter(AgeVac %in% seq(70, 75, 5)) %>% 
#   filter(Yr <= 20) %>% 
#   #filter(Yr %in% c(5, 10)) %>% 
#   group_by(Tag, Yr) %>% 
#   summarise(
#     VE = mean(VE)
#   )


gs$g_vac_ve_zvl <- pars$VE_ZVL %>%  
  mutate(Yr = TimeVac) %>%
  filter(Age - Yr == 70) %>% 
  filter(Yr <= 20)%>% 
  group_by(Yr) %>% 
  summarise(
    M = mean(Protection),
    L = quantile(Protection, 0.025),
    U = quantile(Protection, 0.975)
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Yr, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Yr, y = M)) +
  geom_pointrange(data = dat_ve %>% filter(Source == "Klein 2023"), aes(x = Yr, y = M, ymin = L, ymax = U)) +
  scale_y_continuous("Protection", label = scales::percent) +
  coord_cartesian(xlim = c(0, 21), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous("Year of vaccination") +
  labs(subtitle = "ZVL vaccine effectiveness") + 
  expand_limits(y = 0:1)


## RZV VE -----
load(here::here("pars", "pars_ve_rzv_rw_zlg.rdata"))

d <- pars_ve_rzv%>% 
  filter(Yr <= 20) %>% 
  #filter(Yr %in% c(5, 10)) %>% 
  group_by(Vaccine, Yr) %>% 
  summarise(
    VE = mean(VE)
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


gs$g_vac_ve_rzv <- pars_ve_rzv %>% 
  filter(Yr <= 20) %>% 
  group_by(Yr) %>% 
  summarise(
    M = median(apply_lor(VE, -lor_rw)),
    L = quantile(apply_lor(VE, -lor_rw), 0.025),
    U = quantile(apply_lor(VE, -lor_rw), 0.975)
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Yr, ymin = L, ymax = U), alpha = 0.2) +
  #geom_line(aes(x = Yr, y = M)) +
  geom_pointrange(data = dat_ve %>% filter(!Realworld), aes(x = Yr, y = M, ymin = L, ymax = U)) +
  geom_line(data = ve_rzv, aes(x = Yr, y = VE, colour = Tag)) +
  scale_x_continuous("Year of vaccination") +
  scale_y_continuous("Protection", label = scales::percent) +
  coord_cartesian(xlim = c(0, 21), ylim = c(0, 1), expand = FALSE) +
  labs(subtitle = "RZV vaccine efficacy / effectiveness") +
  scale_color_discrete("") +
  expand_limits(y = 0, x = 0) +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
    legend.background = element_blank(), legend.text = element_text(size = 14))

gs$g_vac_ve_rzv 


## ZVL uptake

load(here::here("data", "processed_vaccine", "coverage.rdata"))
# load(here::here("pars", "fitted_coverage.rdata"))
load(here::here("pars", "fitted_vac_uptake.rdata"))


dat <- coverage %>% 
  mutate(
    VacT = Age - 71,
    Cohort = as.numeric(Year) - VacT
  ) %>% 
  filter(Cohort >= 2014) %>% 
  mutate(Cohort = as.character(Cohort))

# 
# gs$g_uptake_gof <- pred1$pred %>% 
#   filter(Cohort == "2014") %>% 
#   ggplot(aes(x = Age - 1)) +
#   geom_ribbon(aes(ymin = Coverage_l, ymax = Coverage_u), alpha = 0.2) +
#   geom_line(aes(y = Coverage)) +
#   geom_point(data = d, aes(x = Age, y = value, colour = as.character(Cohort)), size = rel(2)) +
#   geom_line(data = d, aes(x = Age, y = value, colour = as.character(Cohort))) +
#   scale_y_continuous("Coverage", labels = scales::percent) +
#   coord_cartesian(xlim = c(69.5, 79.5), ylim = c(0, 1), expand = FALSE) +
#   scale_colour_discrete("Cohort at age 70") +
#   scale_x_continuous("Age", breaks = seq(70, 80, 2)) +
#   expand_limits(y = c(0, 1)) +
#   labs(subtitle = "Vaccine uptake") +
#   theme(legend.position = c(1, -0.03), legend.justification = c(1.1, -0.1),
#     legend.background = element_blank(), legend.text = element_text(size = 14))
# 
# 
# gs$g_uptake_gof


gs$g_uptake_pred <- pred %>% 
  ggplot() +
  stat_lineribbon(aes(x = VacT, y = Coverage), .width = c(.99, .95, .8, .5), color = "#08519C", alpha = 0.3) +
  scale_fill_brewer("Interval") +
  #geom_pointrange(data = dat_ve %>% filter(!Realworld), aes(x = Yr, y = M, ymin = L, ymax = U)) +
  geom_point(data = dat, aes(x = VacT, y = value, colour = as.character(Cohort)), size = rel(2)) +
  geom_line(data = dat, aes(x = VacT, y = value, colour = as.character(Cohort))) +
  scale_y_continuous("Coverage, %", label = scales::percent) +
  scale_x_continuous("Years since vaccinated", breaks = 0:9) +
  scale_colour_discrete("Cohort (year of 70 YOA)") +
  expand_limits(y = 0:1) +
  theme(legend.position = c(1, -0.03), legend.justification = c(1.1, -0.1),
        legend.background = element_blank(), legend.text = element_text(size = 14))


gs$g_uptake_fitted <- fitted %>% 
  ggplot() +
  stat_lineribbon(aes(x = VacT, y = Coverage), .width = c(.99, .95, .8, .5), color = "#08519C", alpha = 0.3) +
  scale_fill_brewer("Interval") +
  #geom_pointrange(data = dat_ve %>% filter(!Realworld), aes(x = Yr, y = M, ymin = L, ymax = U)) +
  geom_point(data = dat, aes(x = VacT, y = value, colour = as.character(Cohort)), size = rel(2)) +
  geom_line(data = dat, aes(x = VacT, y = value, colour = as.character(Cohort))) +
  scale_y_continuous("Coverage, %", label = scales::percent) +
  scale_colour_discrete("Cohort (year of 70 YOA)") +
  scale_x_continuous("Years since vaccinated", breaks = 0:9) +
  expand_limits(y = 0:1) +
  theme(legend.position = c(1, -0.03), legend.justification = c(1.1, -0.1),
        legend.background = element_blank(), legend.text = element_text(size = 14))




#### Binding / output -----

gs$g_bind <- ggarrange(
  ggarrange(
    gs$g_r_hz,
    gs$g_r_hz_gp,
    gs$g_p_phn,
    gs$g_p_mor,
    gs$g_cost_gp,
    gs$g_cost_gpphn,
    gs$g_cost_hosp, 
    gs$g_ql_ph,
    nrow = 3,
    ncol = 3,
    common.legend = T, legend = "bottom"
  ),
  ggarrange(
    gs$g_vac_ve_zvl,
    gs$g_vac_ve_rzv,
    gs$g_uptake_gof,
    nrow = 1,
    ncol = 3
  ),
  nrow = 2,
  heights = c(4, 1) / 5
)
  

gs$g_base <- ggarrange(
  gs$g_r_hz,
  gs$g_r_hz_gp,
  gs$g_p_hosp,
  gs$g_p_phn,
  gs$g_p_mor,
  gs$g_ql_ph,
  gs$g_cost_gp,
  gs$g_cost_gpphn,
  gs$g_cost_hosp, 
  nrow = 3,
  ncol = 3,
  common.legend = T, legend = "bottom", align = "hv",
  labels = LETTERS, font.label = list(size = 16)
)

gs$g_vac <- ggarrange(
  gs$g_vac_ve_zvl,
  gs$g_vac_ve_rzv,
  gs$g_uptake_pred,
  nrow = 3,
  ncol = 1,
  align = "hv",
  labels = LETTERS, font.label = list(size = 16)
)


gs$g_epi <- ggarrange(
  gs$g_r_hz, 
  gs$g_r_hz_gp, 
  gs$g_p_phn, 
  nrow = 3,
  common.legend = T, legend = "bottom"
)


gs$g_cost <- ggarrange(
  gs$g_cost_gp,
  gs$g_cost_gpphn,
  gs$g_cost_hosp, 
  nrow = 3,
  common.legend = T, legend = "bottom"
)



ggsave(gs$g_epi, filename = here::here("docs", "figs", "inputs", "g_epi.png"), width = 5.5, height = 10)
ggsave(gs$g_vac, filename = here::here("docs", "figs", "inputs", "g_vac.png"), width = 7.5, height = 14)
ggsave(gs$g_base, filename = here::here("docs", "figs", "inputs", "g_base.png"), width = 15, height = 12)

ggsave(gs$g_demo, filename = here::here("docs", "figs", "inputs", "g_demo.png"), width = 7.5, height = 10)
ggsave(gs$g_cost, filename = here::here("docs", "figs", "inputs", "g_cost.png"), width = 5.5, height = 10)
ggsave(gs$g_bind, filename = here::here("docs", "figs", "inputs", "g_bind.png"), width = 15, height = 20)

