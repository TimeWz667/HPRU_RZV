library(rstan)
library(tidyverse)
library(tidybayes)
library(loo)

theme_set(theme_bw())


## Load model
options(mc.cores = 6)
rstan_options(auto_write = TRUE)

model <- stan_model(here::here("models", "Coverage.stan"))


## Load data
load(here::here("data", "processed_vaccine", "coverage.rdata"))

dat <- coverage %>% 
  mutate(
    VacT = Age - 71,
    Cohort = as.numeric(Year) - VacT
  ) %>% 
  filter(Cohort >= 2014) %>% 
  mutate(Cohort = as.character(Cohort))



ds <- dat %>% select(Coverage = value, VacT) %>% as.list()

ds$N <- nrow(dat)

ds


## Model fitting
n_collect <- 2000
n_warmup <- 1500
n_chains = 4
n_iter <- n_warmup + round(n_collect / n_chains)

post <- sampling(model, data = ds, 
                 chains = n_chains, iter = n_iter, warmup = n_warmup)


pars <- data.frame(rstan::extract(post, pars = c("p_initial", "p_catch", "sigma"))) %>%
  as_tibble() %>% 
  mutate(Key = 1:n()) %>% 
  filter(Key <= n_collect)


pred <- data.frame(rstan::extract(post, pars = c("Coverage_pred"))) %>%
  as_tibble() %>% 
  mutate(Key = 1:n()) %>% 
  pivot_longer(- "Key", values_to = "Coverage") %>% 
  tidyr::extract(name, "VacT", "Coverage_pred.(\\d+)", convert = T) %>% 
  mutate(VacT = VacT - 1)


fitted <- data.frame(rstan::extract(post, pars = c("Coverage_fitted"))) %>%
  as_tibble() %>% 
  mutate(Key = 1:n()) %>% 
  pivot_longer(- "Key", values_to = "Coverage") %>% 
  tidyr::extract(name, "VacT", "Coverage_fitted.(\\d+)", convert = T) %>% 
  mutate(VacT = VacT - 1)


save(pars, pred, fitted, file = here::here("pars", "fitted_vac_uptake.rdata"))



g_uptake_pred <- pred %>% 
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
  theme(legend.position = c(1, 0), legend.justification = c(1.1, -0.1))


g_uptake_fitted <- fitted %>% 
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
  theme(legend.position = c(1, 0), legend.justification = c(1.1, -0.1))





