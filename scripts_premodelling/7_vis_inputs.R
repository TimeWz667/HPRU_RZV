library(tidyverse)
library(tidybayes)
library(ggpubr)

theme_set(theme_bw())




load(here::here("data", "processed_epi", "Epi_HZ_CPRD_23Nov19.rdata"))
load(here::here("data", "processed_epi", "Epi_HZ_CPRD_GPR.rdata"))


gs <- list()


gs$g_gpr_r_hz <- Epi_HZ %>% 
  ggplot() +
  stat_interval(aes(x = Age, y = r_hz)) +
  geom_pointinterval(data = epi_hz %>% filter(Age <= 90), aes(x = Age, y = M, ymin = L, ymax = U)) +
  scale_colour_brewer() +
  facet_wrap(.~ IC, scales = "free_y", labeller = label_both) +
  scale_y_continuous("Incidence HZ per 1,000", labels = scales::number_format(scale = 1e3)) +
  expand_limits(y = 0)


gs$g_gpr_r_hz_gp <- Epi_HZ %>% 
  ggplot() +
  stat_interval(aes(x = Age, y = p_gp * r_hz)) +
  geom_pointinterval(data = epi_gp %>% filter(Age <= 90), aes(x = Age, y = M, ymin = L, ymax = U)) +
  scale_colour_brewer() +
  facet_wrap(.~ IC, scales = "free_y", labeller = label_both) +
  scale_y_continuous("Incidence HZ, outpatients, per 1,000", labels = scales::number_format(scale = 1e3)) +
  expand_limits(y = 0)


gs$g_gpr_p_phn <- Epi_HZ %>% 
  ggplot() +
  stat_interval(aes(x = Age, y = p_phn)) +
  geom_pointinterval(data = epi_phn %>% filter(Age <= 90), aes(x = Age, y = M, ymin = L, ymax = U)) +
  scale_colour_brewer() +
  facet_wrap(.~ IC, scales = "free_y", labeller = label_both) +
  scale_y_continuous("Postherpetic neuralgia, %", labels = scales::percent) +
  expand_limits(y = 0)


gs$g_gpr_r_hz_nic <- Epi_HZ %>% 
  filter(IC == 0) %>% 
  ggplot() +
  stat_interval(aes(x = Age, y = r_hz)) +
  geom_pointinterval(data = epi_hz %>% filter(IC == 0), aes(x = Age, y = M, ymin = L, ymax = U)) +
  scale_colour_brewer() +
  scale_y_continuous("Incidence HZ per 1,000", labels = scales::number_format(scale = 1e3)) +
  expand_limits(y = 0)


gs$g_gpr_r_hz_gp_nic <- Epi_HZ %>% 
  filter(IC == 0) %>% 
  ggplot() +
  stat_interval(aes(x = Age, y = p_gp * r_hz)) +
  geom_pointinterval(data = epi_gp %>% filter(IC == 0), aes(x = Age, y = M, ymin = L, ymax = U)) +
  scale_colour_brewer() +
  scale_y_continuous("Incidence HZ, GP-only, per 1,000", labels = scales::number_format(scale = 1e3)) +
  expand_limits(y = 0)


gs$g_gpr_p_phn_nic <- Epi_HZ %>% 
  filter(IC == 0) %>% 
  ggplot() +
  stat_interval(aes(x = Age, y = p_phn)) +
  geom_pointinterval(data = epi_phn %>% filter(IC == 0), aes(x = Age, y = M, ymin = L, ymax = U)) +
  scale_colour_brewer() +
  scale_y_continuous("Postherpetic neuralgia, %", labels = scales::percent) +
  expand_limits(y = 0)


gs$g_nlm_r_hz_hosp_nic <- local({
  load(here::here("data", "processed_epi", "Epi_HZ_NIC_CPRD.rdata"))
  
  Hospitalisation_HZ %>% 
    filter(Key <= 300) %>% 
    filter(Age >= 18) %>% 
    ggplot(aes(x = Age, y = r_hospitalisation_hz)) +
    stat_lineribbon(.width = c(.99, .95, .8, .5), color = "#08519C") +
    scale_y_continuous("Incidence HZ, Hospitalised, per 100 000", labels = scales::number_format(scale = 1e5)) +
    expand_limits(y = 0) +
    scale_fill_brewer()
})



gs$g_nlm_r_mor_hz_nic <- local({
  load(here::here("data", "processed_epi", "Epi_HZ_NIC_CPRD_bind.rdata"))
  Epi_HZ %>% 
    select(Key, Age, r_mor_hz) %>% 
    filter(Age >= 50)
}) %>% 
  ggplot(aes(x = Age, y = r_mor_hz)) +
  stat_lineribbon(.width = c(.99, .95, .8, .5), color = "#08519C") +
  scale_y_continuous("Mortality HZ, per 100 000", labels = scales::number_format(scale = 1e5)) +
  expand_limits(y = 0) +
  scale_fill_brewer()


gs$g_nlm_c_hosp_nic <- local({
  load(here::here("data", "processed_ce", "Cost_Hospitalisation_NIC.rdata"))
  
  Cost_Hospitalisation_HZ %>% 
    filter(Age >= 18) %>% 
    ggplot(aes(x = Age, y = Hospitalisation_costs_pp_HZ )) +
    stat_lineribbon(.width = c(.99, .95, .8, .5), color = "#08519C") +
    scale_y_continuous("HZ-related hospitalisation cost \nper patient, GBP") +
    expand_limits(y = 0) +
    scale_fill_brewer()

})



gs$g_gpr_r_hz_gp
gs$g_gpr_r_hz
gs$g_gpr_p_phn

gs$g_gpr_r_hz_nic
gs$g_gpr_r_hz_gp_nic
gs$g_gpr_p_phn_nic


gs$g_nlm_r_mor_hz_nic
gs$g_nlm_c_hosp_nic


gs$g_epi <- ggarrange(
  gs$g_gpr_r_hz_nic + theme(legend.position = "none") + labs(subtitle = "(A)"),
  gs$g_nlm_r_mor_hz_nic + theme(legend.position = "none") + labs(subtitle = "(D)"),
  gs$g_gpr_r_hz_gp_nic + theme(legend.position = "none") + labs(subtitle = "(B)"),
  gs$g_nlm_c_hosp_nic + theme(legend.position = "none") + labs(subtitle = "(E)"),
  gs$g_gpr_p_phn_nic + theme(legend.position = "none") + labs(subtitle = "(C)"),
  nrow = 3, ncol = 2, align = "v"
)

gs$g_epi


# Comparison plots with Rosello & Jit report to JCVI, 2018
# Incidence of cases in primary care
dp1 <- read_csv(here::here("data", "previous", "2019_fig1.csv"))
g_comp_gp <- Epi_HZ %>% 
  filter(IC == 0) %>% 
  ggplot() +
  stat_interval(aes(x = Age, y = r_hz)) +
  geom_ribbon(data = dp1, aes(x = Age, ymin = lower, ymax = upper), fill = "#ec3", alpha = 0.5) +
  geom_line(data = dp1, aes(x = Age, y = middle), colour = "#440") +
  geom_pointinterval(data = epi_hz %>% filter(IC == 0), aes(x = Age, y = M, ymin = L, ymax = U)) +
  scale_colour_brewer() +
  scale_y_continuous("Incidence HZ per 1,000", labels = scales::number_format(scale = 1e3), breaks = c(0, 2, 4, 6, 8, 10, 12) * 0.001) +
  expand_limits(y = 0) +
  theme(legend.position = "none")

# Proportion with PHN
dp4 <-read_csv(here::here("data", "previous", "2019_fig4.csv"))
g_comp_phn <- Epi_HZ %>% 
  filter(IC == 0) %>% 
  ggplot() +
  stat_interval(aes(x = Age, y = p_phn)) +
  geom_ribbon(data = dp4, aes(x = Age, ymin = lower, ymax = upper), fill = "#ec3", alpha = 0.5) +
  geom_line(data = dp4, aes(x = Age, y = middle), colour = "#440") +
  geom_pointinterval(data = epi_phn %>% filter(IC == 0), aes(x = Age, y = M, ymin = L, ymax = U)) +
  scale_colour_brewer() +
  scale_y_continuous("Postherpetic neuralgia, %", labels = scales::percent) +
  expand_limits(y = 0) +
  theme(legend.position = "none")

# Hospitalisation
dp2 <- read_csv(here::here("data", "previous", "2019_fig2.csv"))
dp2 <- dp2[!is.na(dp2$middle), ]
g_comp_hosp <- Epi_HZ %>% 
  filter(IC == 0) %>% 
  ggplot() +
  stat_interval(aes(x = Age, y = r_hz * (1 - p_gp))) +
  geom_ribbon(data = dp2, aes(x = Age, ymin = lower, ymax = upper), fill = "#ec3", alpha = 0.5) +
  geom_line(data = dp2, aes(x = Age, y = middle), colour = "#440") +
  scale_colour_brewer() +
  scale_y_continuous("Hospitalization HZ per 100,000", labels = scales::number_format(scale = 1e5)) +
  expand_limits(y = 0) +
  theme(legend.position = "none")


# Mortality: not simulated in previous report
g_comp_mort <- local({
  load(here::here("data", "processed_epi", "Epi_HZ_NIC_CPRD_bind.rdata"))
  Epi_HZ %>% 
    select(Key, Age, r_mor_hz) %>% 
    filter(Age >= 50)
}) %>% 
  ggplot(aes(x = Age, y = r_mor_hz)) +
  stat_lineribbon(.width = c(.99, .95, .8, .5), color = "#08519C") +
  scale_y_continuous("Mortality HZ per 100,000", labels = scales::number_format(scale = 1e5)) +
  expand_limits(y = 0) +
  scale_fill_brewer() +
  theme(legend.position = "none")

g_comp <- ggarrange(
  g_comp_gp   + labs(subtitle = "(A)"),
  g_comp_phn  + labs(subtitle = "(B)"),
  g_comp_hosp + labs(subtitle = "(C)"),
  g_comp_mort + labs(subtitle = "(D)"),
  nrow = 2, ncol = 2, align = "v"
)

g_comp

ggsave(g_comp, filename = here::here("outputs", "figs", "g_epi_nic_comp.png"), width = 7.5, height = 6.8)




load(here::here("data", "processed_demography", "Population_ONS.rdata"))

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



gs$g_demo <- ggarrange(
    gs$g_demo_pop + labs(subtitle = "(A)"), 
    gs$g_demo_dr + labs(subtitle = "(B)"), 
    nrow = 1, align = "v"
  )



ggsave(gs$g_epi, filename = here::here("outputs", "figs", "g_epi_nic.png"), width = 7.5, height = 10)
ggsave(gs$g_demo, filename = here::here("outputs", "figs", "g_demo_nic.png"), width = 7.5, height = 10)



gs$g_nlm_c_hosp_nic_50plus <- local({
  load(here::here("data", "processed_ce", "Cost_Hospitalisation_NIC.rdata"))
  
  Cost_Hospitalisation_HZ %>% 
    filter(Age >= 50) %>% 
    ggplot(aes(x = Age, y = Hospitalisation_costs_pp_HZ )) +
    stat_lineribbon(.width = c(.99, .95, .8, .5), color = "#08519C") +
    scale_y_continuous("HZ-related hospitalisation cost \nper patient, GBP") +
    expand_limits(y = 0) +
    scale_fill_brewer()

})

ggsave(gs$g_nlm_c_hosp_nic_50plus + theme(legend.position = "none"), filename = here::here("outputs", "figs", "g_hosp_cost.png"), width = 4, height = 3)



