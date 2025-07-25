
## Parameters: Demography

source(here::here("scripts_premodelling", "2_1_prepare_pars_demo.R"))


## Parameters: Epidemiology
source(here::here("scripts_premodelling", "2_2_model_epi_gpr.R"))


## Parameters: Costing
source(here::here("scripts_premodelling", "2_3_prepare_pars_cost.R"))


## Parameters: Vaccine efficacy/effectiveness
source(here::here("scripts_premodelling", "3_1_prepare_ve.R"))
source(here::here("scripts_premodelling", "3_2_fit_ve_rzv_zi_stan.R"))
source(here::here("scripts_premodelling", "3_3_fit_ve_zvl_zi_stan.R"))
source(here::here("scripts_premodelling", "3_4_compile_ve.R"))
source(here::here("scripts_premodelling", "3_5_fit_ve_rzv_sens.R"))


## Parameters: Vaccine coverage
source(here::here("scripts_premodelling", "4_1_vis_coverage.R"))
source(here::here("scripts_premodelling", "4_2_fit_coverage.R"))
source(here::here("scripts_premodelling", "4_3_sim_immunity.R"))


## Legacy simulations
# source(here::here("scripts_premodelling", "1_vis_inputs.R"))
# source(here::here("scripts_premodelling", "5_1_run_shingrix_nic.R"))
# source(here::here("scripts_premodelling", "5_2_run_shingrix_ic.R"))
# source(here::here("scripts_premodelling", "5_3_run_zostavax_nic_aj.R"))
# source(here::here("scripts_premodelling", "5_4_run_shingrix_nic_updated.R"))
# source(here::here("scripts_premodelling", "5_5_run_shingrix_ic_updated.R"))
# 
# source(here::here("scripts_premodelling", "6_compare_versions.R"))
# source(here::here("scripts_premodelling", "6_validate_with_previous_version.R"))
