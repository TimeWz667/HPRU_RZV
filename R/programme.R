
strategy_null <- function(df, year) {
  require(tidyverse)
  return(df %>% mutate(eli = NA, tp_uptake = 0))
}


strategy_zvl <- function(df, year) {
  require(tidyverse)

  df <- df %>% 
    mutate(
      eli = case_when(
        Vaccine != "None" ~ NA,
        Age < 70 ~ NA,
        Age < 80 ~ "ZVL",
        T ~ NA
      ),
      tp_uptake = case_when(
        is.na(eli) ~ NA,
        Age == 70 ~ "Ini",
        T ~ "Cat"
      )
    )
  
  return(df)
}


strategy_changeonly <- function(df, year) {
  require(tidyverse)

  df <- df %>% 
    mutate(
      eli = case_when(
        Vaccine != "None" ~ NA,
        Age < 70 ~ NA,
        Age < 80 ~ ifelse(year < 2023, "ZVL", "RZV_2d"),
        T ~ NA
      ),
      tp_uptake = case_when(
        is.na(eli) ~ NA,
        Age == 70 ~ "Ini",
        T ~ "Cat"
      )
    )
  
  return(df)
}


strategy_scheduled <- function(df, year) {
  require(tidyverse)
  
  if (year < 2023) {
    df <- df %>% 
      mutate(
        eli = case_when(
          Vaccine != "None" ~ NA,
          Age < 70 ~ NA,
          Age < 80 ~ "ZVL",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age == 70 ~ "Ini",
          T ~ "Cat"
        )
      )
  } else if (year < 2028) {
    df <- df %>% 
      mutate(
        eli = case_when(
          Vaccine != "None" ~ NA,
          Age < 65 ~ NA,
          Age >= 80 ~ NA,
          Age >= 70 ~ "RZV_2d",
          Age < (65 + year - 2023) ~ "RZV_2d",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age %in% c(65, 70) ~ "Ini",
          T ~ "Cat"
        )
      )
  } else if (year < 2033) {
    df <- df %>% 
      mutate(
        eli = case_when(
          Vaccine != "None" ~ NA,
          Age < 60 ~ NA,
          Age >= 80 ~ NA,
          Age >= 65 ~ "RZV_2d",
          Age < (60 + year - 2028) ~ "RZV_2d",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age %in% c(60, 65) ~ "Ini",
          T ~ "Cat"
        )
      )
  } else {
    df <- df %>% 
      mutate(
        eli = case_when(
          Vaccine != "None" ~ NA,
          Age < 60 ~ NA,
          Age < 80 ~ "RZV_2d",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age == 60 ~ "Ini",
          T ~ "Cat"
        )
      )
  }
}


strategy_scheduled65 <- function(df, year) strategy_scheduled(df, min(year, 2027))


strategy_scheduled_old <- function(df, year, year0 = 2024, vaccine_old = "1d", cap_age = 85) {
  require(tidyverse)
  
  if (vaccine_old == "1d") {
    vo_pri <- "RZV_1d"
    vo_re <- "ReRZV_1d"
  } else if (vaccine_old == "2d") {
    vo_pri <- "RZV_2d"
    vo_re <- "ReRZV_2d"
  } else {
    vo_pri <- NA
    vo_re <- NA
  }
  
  
  if (year < 2023) {
    df <- df %>% 
      mutate(
        eli = case_when(
          Vaccine != "None" ~ NA,
          Age < 70 ~ NA,
          Age < 80 ~ "ZVL",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age == 70 ~ "Ini",
          T ~ "Cat"
        )
      )
  } else if (year < 2028) {
    df <- df %>% 
      mutate(
        eli = case_when(
          Age >= cap_age ~ NA,
          (Age >= 80) & (Vaccine == "None") ~ vo_pri,
          (Age >= 80) & (Vaccine == "ZVL") ~ vo_re,
          Age < 65 ~ NA,
          Vaccine != "None" ~ NA,
          Age >= 70 ~ "RZV_2d",
          Age < (65 + year - 2023) ~ "RZV_2d",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age %in% c(65, 70, 80) ~ "Ini",
          T ~ "Cat"
        )
      )
  } else if (year < 2033) {
    df <- df %>% 
      mutate(
        eli = case_when(
          Age >= cap_age ~ NA,
          (Age >= 80) & (Vaccine == "None") ~ vo_pri,
          (Age >= 80) & (Vaccine == "ZVL") ~ vo_re,
          Age < 60 ~ NA,         
          Vaccine != "None" ~ NA,
          Age >= 65 ~ "RZV_2d",
          Age < (60 + year - 2028) ~ "RZV_2d",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age %in% c(60, 65, 80) ~ "Ini",
          T ~ "Cat"
        )
      )
  } else {
    df <- df %>% 
      mutate(
        eli = case_when(
          Age >= cap_age ~ NA,
          (Age >= 80) & (Vaccine == "None") ~ vo_pri,
          (Age >= 80) & (Vaccine == "ZVL") ~ vo_re,
          Age < 60 ~ NA,         
          Vaccine != "None" ~ NA,
          Age < 80 ~ "RZV_2d",
          T ~ NA
        ),
        tp_uptake = case_when(
          is.na(eli) ~ NA,
          Age %in% c(60, 80) ~ "Ini",
          T ~ "Cat"
        )
      )
  }
}


strategy_scheduled_1d85 <- function(df, year) strategy_scheduled_old(df, year, year0 = 2024, vaccine_old = "1d", cap_age = 85) 
strategy_scheduled_1d90 <- function(df, year) strategy_scheduled_old(df, year, year0 = 2024, vaccine_old = "1d", cap_age = 90) 
strategy_scheduled_1d95 <- function(df, year) strategy_scheduled_old(df, year, year0 = 2024, vaccine_old = "1d", cap_age = 95)

strategy_scheduled_2d85 <- function(df, year) strategy_scheduled_old(df, year, year0 = 2024, vaccine_old = "2d", cap_age = 85) 
strategy_scheduled_2d90 <- function(df, year) strategy_scheduled_old(df, year, year0 = 2024, vaccine_old = "2d", cap_age = 90) 
strategy_scheduled_2d95 <- function(df, year) strategy_scheduled_old(df, year, year0 = 2024, vaccine_old = "2d", cap_age = 95) 

