library(tidyverse)

scenario <- "rw_35"

tab_ce_re <- read_csv(here::here("docs", "tabs", scenario, "stats_re_ce.csv"))
tab_ce_uv <- read_csv(here::here("docs", "tabs", scenario, "stats_uv_ce.csv"))


bind_rows(
  tab_ce_uv %>% 
    filter(Arm == "RZV_2d") %>% 
    filter(Index == "Thres") %>% 
    select(Age = Age0, Arm, Thres = M) %>% 
    filter(Age %in% seq(60, 100, 5)),
  tab_ce_uv %>% 
    filter(Arm == "RZV_1d") %>% 
    filter(Index == "Thres") %>% 
    select(Age = Age0, Arm, Thres = M) %>% 
    filter(Age %in% seq(60, 100, 5)),
  tab_ce_re %>% 
    filter(Scenario != "Overall") %>% 
    #filter(Arm == "ReRZV_1d") %>% 
    filter(Index == "Thres") %>% 
    select(Age0, Age = Age1, Arm, Thres = M) %>% 
    filter(Age %in% seq(60, 100, 5))
) %>% 
  relocate(Age0, Age, Arm, Thres) %>%
  arrange(Arm, Age) %>% 
  write_csv(here::here("docs", "manuscript", "threshold_price.csv"))

