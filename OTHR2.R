#table generation 

#**************************************************************
# TABLE: GESTATIONAL AGE-SPECIFIC REFERENCE VALUES
# All Biomarkers - Box-Cox FPR
# Weeks 5, 10, 15, 20, 25, 30, 35
#**************************************************************

library(dplyr)
library(gt)
library(gamlss)
library(gamlss.dist)

selected_weeks <- seq(10, 35, by = 5)
#**************************************************************
# EXTRACT EXISTING FITTED REFERENCE VALUES
#**************************************************************

library(dplyr)

selected_weeks <- seq(10, 35, by = 5)

reference_values <- bind_rows(
  
  centile_agp_bc %>%
    mutate(Biomarker = "AGP"),
  
  centile_nlr_bc %>%
    mutate(Biomarker = "NLR"),
  
  centile_lmr_bc %>%
    mutate(Biomarker = "LMR"),
  
  centile_plr_bc %>%
    mutate(Biomarker = "PLR"),
  
  centile_sii_bc %>%
    mutate(Biomarker = "SII")
  
) %>%
  filter(
    dplyr::near(ga_wks, round(ga_wks)),
    round(ga_wks) %in% selected_weeks
  ) %>%
  mutate(
    ga_wks = round(ga_wks)
  ) %>%
  select(
    Biomarker,
    ga_wks,
    centile_5,
    centile_50,
    centile_95
  ) %>%
  arrange(Biomarker, ga_wks) %>%
  mutate(
    across(
      c(centile_5, centile_50, centile_95),
      ~ round(.x, 2)
    )
  )

# Display table
reference_values
