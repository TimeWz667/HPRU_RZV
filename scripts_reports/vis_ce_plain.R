library(targets)
library(tidyverse)


theme_set(theme_bw())



yss <- tar_read(yss_uv, 3)



sel <- yss %>% 
  filter(Age0 %in% c(80, 85, 90)) %>% 
  select(Key, Age0, Arm, C_Med_d, Q_All_d, N_VacRZV_d)


sel %>% 
  filter(Arm != "SOC") %>% 
  left_join(sel %>% filter(Arm == "SOC") %>% select(-Arm), by = c("Key", "Age0"), suffix = c("_1", "_0")) %>% 
  crossing(price = c(75, 100, 125)) %>% 
  mutate(
    dC_All_d = price * N_VacRZV_d_1 + C_Med_d_1 - C_Med_d_0,
    dQ_All_d = Q_All_d_1 - Q_All_d_0
  ) %>% 
  filter(Key <= 500) %>% 
  ggplot() +
  geom_point(aes(x = dQ_All_d, y = dC_All_d, colour = Arm), alpha = 0.1) +
  geom_abline(slope = 2e4, linetype = 2) +
  geom_abline(slope = 3e4, linetype = 3) +
  scale_y_continuous("Incremental cost per vaccinator, £") + 
  scale_colour_discrete("Scenario", label = c(RZV_1d = "One-dose", RZV_2d = "Two-dose")) +
  facet_grid(price~Age0) + 
  expand_limits(x = 0, y = 0)


