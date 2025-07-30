library(tidyverse)
library(tidybayes)
library(ggpubr)


theme_seggpubrtheme_set(theme_bw())


load(here::here("pars", "pars_base_35_rw.rdata"))

names(pars)



g_ip <- crossing(Key = 1:100, Age = 60:95, AgeVac = c(60, 65, 70)) %>% 
  mutate(TimeVac = Age - AgeVac + 1) %>% 
  left_join(pars$VE_RZV_2d, relationship = "many-to-many") %>% 
  mutate(
    Protection = ifelse(Age >= AgeVac, Protection, 0),
    AgeVac = as.character(AgeVac)) %>% 
  group_by(Age, AgeVac) %>% 
  summarise(IP = mean(Protection)) %>% 
  ggplot() +
  geom_bar(aes(x = Age, y = IP, fill = AgeVac), stat = "identity", position = "dodge", width = 1) +
  scale_y_continuous("Immunity to HZ (vaccine effectiveness), %", labels = scales::percent) +
  scale_fill_discrete("Age of RZV vaccination") +
  expand_limits(y = 0:1) +
  theme(legend.position = c(1, 1), legend.justification = c(1.1, 1.1))

g_ip
  

g_inc <- crossing(Key = 1:100, Age = 60:95) %>% 
  left_join(pars$Epidemiology, relationship = "many-to-many") %>% 
  group_by(Age) %>% 
  summarise(r_hz = mean(r_hz)) %>% 
  ggplot() +
  geom_bar(aes(x = Age, y = r_hz), stat = "identity", position = "dodge", width = 1, alpha = 0.3) +
  scale_y_continuous("HZ incidence, all, per 100,000", labels = scales::number_format(scale = 1e5)) +
  expand_limits(y = 0.01)

g_inc



d_inc <- crossing(Key = 1:100, Age = 60:95, AgeVac = c(60, 65, 70)) %>% 
  mutate(TimeVac = Age - AgeVac + 1) %>% 
  left_join(pars$VE_RZV_2d, relationship = "many-to-many") %>% 
  left_join(pars$Epidemiology, relationship = "many-to-many") %>% 
  mutate(
    Protection = ifelse(Age >= AgeVac, Protection, 0),
    AgeVac = as.character(AgeVac)
  ) %>% 
  group_by(Age, AgeVac) %>% 
  summarise(
    r_hz_re = mean(r_hz * (1 - Protection)),
    r_hz_im = mean(r_hz * Protection),
    r_hz = mean(r_hz)
  )


g_inc <- ggplot(d_inc) +
  geom_bar(aes(x = Age, y = r_hz, fill = AgeVac), stat = "identity", position = "dodge", width = 1) +
  scale_y_continuous("HZ incidence, all, per 100,000", labels = scales::number_format(scale = 1e5)) +
  scale_fill_discrete("Age of RZV vaccination") +
  expand_limits(y = 0.01) +
  theme(legend.position = "bottom")


g_inc_re <- ggplot(d_inc) +
  geom_bar(aes(x = Age, y = r_hz_re, fill = AgeVac), stat = "identity", position = "identity", width = 1, alpha = 0.3) +
  scale_y_continuous("HZ incidence, all, per 100,000", labels = scales::number_format(scale = 1e5)) +
  scale_fill_discrete("Age of RZV vaccination") +
  expand_limits(y = 0.01) +
  theme(legend.position = "none", legend.justification = c(1.1, 1.1))


g_inc_im <- ggplot(d_inc) +
  geom_bar(aes(x = Age, y = r_hz_im, fill = AgeVac), stat = "identity", position = "identity", width = 1, alpha = 0.3) +
  scale_y_continuous("HZ incidence, all, per 100,000", labels = scales::number_format(scale = 1e5)) +
  scale_fill_discrete("Age of RZV vaccination") +
  expand_limits(y = 0.01) +
  theme(legend.position = "none", legend.justification = c(1.1, 1.1))




g1 <- ggarrange(g_ip, g_inc, nrow = 2)

g2 <- ggarrange(
  g_inc_re + labs(subtitle = "Breakthrough cases"), 
  g_inc_im + labs(subtitle = "Averted cases"),  
  nrow = 2)

g1
g2


g_31 <- d_inc %>% 
  select(Age, AgeVac, r_hz_im) %>% 
  filter(AgeVac %in% c(60, 70)) %>% 
  ggplot() +
  geom_bar(aes(x = Age, y = r_hz_im, fill = AgeVac), stat = "identity", position = "identity", width = 1, alpha = 0.2) +
  scale_y_continuous("HZ incidence, all, per 100,000", labels = scales::number_format(scale = 1e5)) +
  scale_fill_discrete("Age of RZV vaccination") +
  expand_limits(y = 0.01) +
  theme(legend.position = "none", legend.justification = c(1.1, 1.1)) +
  labs(subtitle = "Averted cases from those vaccinated at their 60 and 70 YOA")


g_32 <- d_inc %>% 
  select(Age, AgeVac, r_hz_im) %>% 
  filter(AgeVac %in% c(60, 70)) %>% 
  ggplot() +
  geom_bar(aes(x = Age, y = r_hz_im, fill = AgeVac), stat = "identity", position = "identity", width = 1, alpha = 0.2) +
  geom_bar(data = d_inc %>% filter(AgeVac == 60) %>% filter(Age >= 70), aes(x = Age, y = r_hz_im), stat = "identity", position = "identity", width = 1, fill = "white") +
  scale_y_continuous("HZ incidence, all, per 100,000", labels = scales::number_format(scale = 1e5)) +
  scale_fill_discrete("Age of RZV vaccination") +
  expand_limits(y = 0.01) +
  theme(legend.position = "none", legend.justification = c(1.1, 1.1)) +
  labs(subtitle = " ")


g_33 <- d_inc %>% 
  select(Age, AgeVac, r_hz_im) %>% 
  filter(AgeVac %in% c(60, 70)) %>% 
  pivot_wider(names_from = AgeVac, values_from = r_hz_im, names_prefix = "V") %>% 
  mutate(
    V60u = pmax(V60 - V70, 0),
    V70u = pmax(V70 - V60, 0)
  ) %>% 
  pivot_longer(c(V60u, V70u), names_to = "AgeVac", names_pattern = "V(\\d+)u", values_to = "diff_inc") %>% 
  ggplot() +
  geom_bar(aes(x = Age, y = diff_inc, fill = AgeVac), stat = "identity", position = "identity", width = 1, alpha = 0.2) +
  scale_y_continuous("HZ incidence, all, per 100,000", labels = scales::number_format(scale = 1e5)) +
  scale_fill_discrete("Age of RZV vaccination") +
  expand_limits(y = 0.01) +
  theme(legend.position = "none", legend.justification = c(1.1, 1.1)) +
  labs(subtitle = "Difference between those vaccinated at their 60 and 70 YOA")



g3 <- ggarrange(g_31, g_32, g_33, ncol = 1)

g3


ggsave(g1, filename = here::here("docs", "immuprofile", "g1.png"), width = 5, height = 8)
ggsave(g2, filename = here::here("docs", "immuprofile", "g2.png"), width = 5, height = 8)
ggsave(g3, filename = here::here("docs", "immuprofile", "g3.png"), width = 5, height = 8)






stats_uv_ys <- read_csv(here::here("docs", "tabs", "rw_35", "stats_uv_ys.csv"))
stats_uv_ce <- read_csv(here::here("docs", "tabs", "rw_35", "stats_uv_ce.csv"))




stats_uv_ys %>% 
  filter(Age0 <= 70) %>% 
  filter(Index == "Year_Immunised") %>% 
  filter(Arm %in% c("RZV_2d")) %>% 
  data.frame()

 

g_41 <- stats_uv_ys %>% 
  filter(Age0 %in% c(60, 65, 70)) %>% 
  filter(Index == "Year_Immunised") %>% 
  filter(Arm == "RZV_2d") %>% 
  ggplot() +
  geom_pointrange(aes(x = as.character(Age0), y = M, ymin = L, ymax = U)) +
  scale_x_discrete("Age of vaccination") +
  scale_y_continuous("Years with full immunisation") +
  expand_limits(y = 0)


ggsave(g_41, filename = here::here("docs", "immuprofile", "g41.png"), width = 5, height = 4)



g_42 <- stats_uv_ce %>% 
  filter(Age0 %in% c(60, 65, 70)) %>% 
  filter(Index %in% c("dC_GP_d", "dC_Hosp_d")) %>% 
  filter(Arm == "RZV_2d") %>% 
  ggplot() +
  geom_pointrange(aes(x = as.character(Age0), y = -M, ymin = -L, ymax = -U)) +
  facet_wrap(.~Index, scale = "free_y", labeller = labeller(Index = c(dC_GP_d = "Cost, GP", dC_Hosp_d = "Cost, Hospitalisation"))) +
  scale_x_discrete("Age of vaccination") +
  scale_y_continuous("GBP per vaccinatee") +
  expand_limits(y = 0)

ggsave(g_42, filename = here::here("docs", "immuprofile", "g42.png"), width = 8, height = 4)



g_43 <- stats_uv_ce %>% 
  filter(Age0 %in% c(60, 65, 70)) %>% 
  filter(Index %in% c("dQ_HZ_d")) %>% 
  filter(Arm == "RZV_2d") %>% 
  ggplot() +
  geom_pointrange(aes(x = Age0, y = M, ymin = L, ymax = U)) +
  scale_x_discrete("Age of vaccination") +
  scale_y_continuous("QALY loss") +
  expand_limits(y = 0)

ggsave(g_43, filename = here::here("docs", "immuprofile", "g43.png"), width = 5, height = 4)




