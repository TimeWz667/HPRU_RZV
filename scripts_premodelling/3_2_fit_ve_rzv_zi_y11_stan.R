library(rstan)
library(tidyverse)
library(readxl)
library(loo)

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



n_collect <- 2000
log_lik <- list()

for (key_model in c("zie", "zig")) {
  if (key_model == "zie") {
    model <-  stan_model(here::here("models", "zi_exp.stan"))
  } else {
    model <-  stan_model(here::here("models", "zi_gamma.stan"))
  }
  
  for (skip910 in c(F, T)) {
    if (skip910) {
      post <- sampling(model, data = dsm, chains = 3, iter = 2000, warmup = floor(2000 - 1000))
    } else {
      post <- sampling(model, data = ds, chains = 3, iter = 2000, warmup = floor(2000 - 1000))
    }
    
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
      #geom_pointrange(data = dat_ve %>% filter(!Realworld), aes(x = Yr, y = M, ymin = L, ymax = U)) +
      scale_y_continuous("Vaccine efficacy, %", label = scales::percent) +
      scale_x_continuous("Year since vaccinated") +
      expand_limits(y = 0)
    
    
    if (skip910) {
      tag <- glue::as_glue("y11m") + "_" + key_model
      g_gof <- g_gof + 
        geom_pointrange(data = dat_ve %>% filter(!(Yr %in% c(9, 10))), aes(x = Yr, y = M, ymin = L, ymax = U)) + 
        geom_pointrange(data = dat_ve %>% filter((Yr %in% c(9, 10))), aes(x = Yr, y = M, ymin = L, ymax = U), linetype = 2, shape = 1)
      
    } else {
      tag <- glue::as_glue("y11") + "_" + key_model
      g_gof <- g_gof + 
        geom_pointrange(data = dat_ve, aes(x = Yr, y = M, ymin = L, ymax = U))
    }
    
    
    log_lik[[tag]] <- extract_log_lik(post)
    
    save(sel, file = here::here("pars", "fitted_ve_rzv_" + tag + ".rdata"))
    # write_csv(sims, here::here("pars", "sims_ve_rzv_" + tag + ".csv"))
    ggsave(g_gof, filename = here::here("docs", "figs", "inputs", "g_pars_ve_rzv_" + tag + ".png"), width = 7, height = 5.5)
  }
}


print(names(log_lik)[c(1, 3)])
loo_compare(loo::waic(log_lik[[1]]), loo::waic(log_lik[[3]]))
loo_compare(loo::loo(log_lik[[1]]), loo::loo(log_lik[[3]]))

print(names(log_lik)[c(2, 4)])
loo_compare(loo::waic(log_lik[[2]]), loo::waic(log_lik[[4]]))
loo_compare(loo::loo(log_lik[[2]]), loo::loo(log_lik[[4]]))


## Zero inflated exponential distributions are better than gammas


