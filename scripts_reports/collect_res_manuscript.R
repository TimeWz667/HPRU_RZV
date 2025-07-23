library(tidyverse)


dir.create(here::here("docs", "manuscript"), showWarnings = T)


scenario <- "rw_35"

root <- here::here("docs", "manuscript")

# Figures
file.copy(here::here("docs", "figs", "inputs", "g_base.png"), here::here(root, "Fig 1.png"), overwrite = T)
file.copy(here::here("docs", "figs", "inputs", "g_vac.png"), here::here(root, "Fig 2.png"), overwrite = T)

file.copy(here::here("docs", "figs", scenario, "g_thres_panel.png"), here::here(root, "Fig 3.png"), overwrite = T)
file.copy(here::here("docs", "figs", scenario, "g_rzv_mb_stack_ann.png"), here::here(root, "Fig 4.png"), overwrite = T)

file.copy(here::here("docs", "figs", "proj", "g_rw_inc_bind.png"), here::here(root, "Fig 5.png"), overwrite = T)

# Tables
file.copy(here::here("docs", "tabs", "Tab_programme_realworld.xlsx"), here::here(root, "Table 1.xlsx"), overwrite = T)

