library(rstan)
library(tidyverse)
library(readxl)

theme_set(theme_bw())


options(mc.cores = 6)
rstan_options(auto_write = TRUE)


dat_ve <- read_xlsx(here::here("data", "processed_vaccine", "VE.xlsx"), sheet = 1) %>% 
  filter(Use) %>% 
  filter(Type == "HZ") %>% 
  filter(Vaccine == "Zostavax") %>% 
  filter(Source == "Klein 2023") %>% 
  mutate(
    M = M / 100,
    L = L / 100, 
    U = U / 100,
    Yr = gsub("yr", "", `Sub-group`),
    Yr = ifelse(is.na(Yr), 1.5, as.numeric(Yr))
  )


dat_ve %>% 
  mutate(
    sd = (U - L) / 2 / 1.96,
    n = M * (1 - M) / sd ** 2,
    n = round(n),
    l = qbinom(0.025, size =n, M) / n,
    u = qbinom(0.975, size =n, M) / n
  )


ds <- dat_ve %>% 
  filter(Follow_up == 1) %>% 
  mutate(
    sd = (U - L) / 2 / 1.96,
    n = M * (1 - M) / sd ** 2,
    n = round(n / n()),
    y = round(n * M)
  ) %>% 
  select(n, yr = Yr, y) %>% 
  as.list()

ds$N <- length(ds$yr)


n_collect <- 2000
log_lik <- list()

for (key_model in c("zie", "zig")) {
  if (key_model == "zie") {
    model <-  stan_model(here::here("models", "zi_exp.stan"))
  } else {
    model <-  stan_model(here::here("models", "zi_gamma.stan"))
  }
  
  post <- sampling(model, data = ds, chains = 3, iter = 2000, warmup = floor(2000 - 1000))
  
  
  if (key_model == "zig") {
    sel <- data.frame(rstan::extract(post, pars = c("p0", "alpha", "beta"))) %>%
      as_tibble()
  } else {
    sel <- data.frame(rstan::extract(post, pars = c("p0", "beta"))) %>%
      as_tibble() %>%
      mutate(alpha = 1)
  }
  sel <- sel %>% mutate(Key = 1:n()) %>% filter(Key <= n_collect)
  
  sims <- sel %>%
    full_join(crossing(Key = (sel %>% pull(Key) %>% unique()), Yr = 1:30)) %>%
    mutate(
      VE = p0 * (1 - pgamma(Yr, alpha, beta))
    )  %>%
    select(Key, Yr, VE)
  
  
  g_gof <- sims %>% 
    ggplot() +
    stat_lineribbon(aes(x = Yr, y = VE), .width = c(.99, .95, .8, .5), color = "#08519C", alpha = 0.6) +
    scale_fill_brewer("Interval") +
    geom_pointrange(data = dat_ve, aes(x = Yr, y = M, ymin = L, ymax = U)) +
    scale_y_continuous("Vaccine efficacy, %", label = scales::percent) +
    scale_x_continuous("Year since vaccinated") +
    expand_limits(y = 0)
  
  
  
  tag <- glue::as_glue(key_model)
  
  log_lik[[tag]] <- extract_log_lik(post)
  
  save(sel, file = here::here("pars", "fitted_ve_zvl_" + tag + ".rdata"))
  # write_csv(sims, here::here("pars", "sims_ve_zvl_" + tag + ".csv"))
  ggsave(g_gof, filename = here::here("docs", "figs", "inputs", "g_pars_ve_zvl_" + tag + ".png"), width = 7, height = 5.5)
  
}


loo_compare(loo::waic(log_lik[[1]]), loo::waic(log_lik[[2]]))
loo_compare(loo::loo(log_lik[[1]]), loo::loo(log_lik[[2]]))

