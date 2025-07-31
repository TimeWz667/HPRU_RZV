library(rstan)
library(tidyverse)
library(readxl)

theme_set(theme_bw())


options(mc.cores = 6)
rstan_options(auto_write = TRUE)


## Data loading
dat <- read_xlsx(here::here("data", "processed_vaccine", "VE_Strezova.xlsx")) %>% 
  mutate(
    Year = gsub("Year ", "", Year),
    Year = as.numeric(Year),
    n = n0 + n1,
    across(c(M, LL, UL), \(d) gsub("·", ".", d) %>% as.double() / 100)
  ) %>% 
  select(Yr = Year, Age_group, n, M = M, L = LL, U = UL)


dat_ve <- dat %>% filter(Age_group == "Overall")


ds <- dat_ve %>% 
  mutate(
    sd = (U - L) / 2 / 1.96,,
    n = round(n / n()),
    y = round(n * M)
  ) %>% 
  select(n, yr = Yr, y) %>% 
  #filter(yr > 1) %>% 
  as.list()

ds$N <- length(ds$yr)


## No y9 y10

dsm <- dat_ve %>% 
  mutate(
    sd = (U - L) / 2 / 1.96,,
    n = round(n / n()),
    y = round(n * M)
  ) %>% 
  filter(Yr != 9 & Yr != 10) %>% 
  select(n, yr = Yr, y) %>% 
  #filter(yr > 1) %>% 
  as.list()

dsm$N <- length(dsm$yr)



for (key_model in c("zl_exp", "zl_gamma")) {
  model <-  stan_model(here::here("models", key_model + glue::as_glue(".stan")))
  
  for (skip910 in c(F, T)) {
    if (skip910) {
      post <- sampling(model, data = dsm, chains = 3, iter = 2000, warmup = floor(2000 - 1000))
    } else {
      post <- sampling(model, data = ds, chains = 3, iter = 2000, warmup = floor(2000 - 1000))
    }
    
    if (key_model == "zl_gamma") {
      sel <- data.frame(rstan::extract(post, pars = c("p0", "alpha", "beta"))) %>%
        as_tibble()
    } else if (key_model == "zl_exp") {
      sel <- data.frame(rstan::extract(post, pars = c("p0", "beta"))) %>%
        as_tibble() %>%
        mutate(alpha = 1)
    }
    sel <- sel %>% mutate(Key = 1:n())
    
    sims <- sel %>%
      filter(Key <= 1000) %>%
      full_join(crossing(Key = 1:1000, Yr = 1:30)) %>%
      mutate(
        VE = p0 * (1 - pgamma(Yr, alpha, beta))
      )  %>%
      select(Key, Yr, VE)
    
    
    g_gof <- sims %>% 
      group_by(Yr) %>% 
      summarise(
        M = mean(VE),
        L = quantile(VE, 0.025),
        U = quantile(VE, 0.975)
      ) %>% 
      ggplot() +
      geom_ribbon(aes(x = Yr, ymin = L, ymax = U), alpha = 0.2) +
      geom_line(aes(x = Yr, y = M)) +
      geom_pointrange(data = dat_ve, aes(x = Yr, y = M, ymin = L, ymax = U)) +
      scale_y_continuous("Vaccine efficacy, %", label = scales::percent) +
      scale_x_continuous("Year since vaccinated") +
      expand_limits(y = 0)
    
    
    if (skip910) {
      tag <- glue::as_glue("y11m")
    } else {
      tag <- glue::as_glue("y11")
    }
    
    
    save(sel, file = here::here("pars", "pars_ve_rzv_" + tag + "_" + glue::as_glue(key_model) + ".rdata"))
    write_csv(sims, here::here("pars", "sims_ve_rzv_" + tag + "_" + glue::as_glue(key_model) + ".csv"))
    ggsave(g_gof, filename = here::here("docs", "figs", "inputs", "g_pars_ve_rzv_" + tag + "_" + glue::as_glue(key_model) + ".png"), width = 7, height = 5.5)
  }
}






