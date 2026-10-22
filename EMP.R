#**************************************************************
# Empirical CRP percentiles and trimester participation summary
#**************************************************************

library(dplyr)
library(tidyr)
library(gt)


#--------------------------------------------------------------
# 1. Prepare trimester variable
#--------------------------------------------------------------

crp_trimester <- crp.model %>%
  filter(
    !is.na(imputed_crp),
    !is.na(ga_wks)
  ) %>%
  mutate(
    trimester = case_when(
      ga_wks < 14 ~ "T1",
      ga_wks < 28 ~ "T2",
      ga_wks <= 40 ~ "T3",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(trimester)) %>%
  mutate(
    trimester = factor(
      trimester,
      levels = c("T1", "T2", "T3")
    )
  )


#--------------------------------------------------------------
# 2. Total number of unique pregnancies
#--------------------------------------------------------------

total_participants <- crp_trimester %>%
  summarise(
    total = n_distinct(pregid)
  ) %>%
  pull(total)


#--------------------------------------------------------------
# 3. Number and percentage of pregnancies represented
#    in each trimester
#
# Percentages may sum to more than 100% because the same
# pregnancy may contribute measurements to multiple trimesters.
#--------------------------------------------------------------

participants_by_trimester <- crp_trimester %>%
  group_by(trimester) %>%
  summarise(
    participants = n_distinct(pregid),
    .groups = "drop"
  ) %>%
  mutate(
    percent = round(
      100 * participants / total_participants,
      1
    )
  )

participants_by_trimester


#--------------------------------------------------------------
# 4. Pregnancies represented in all three trimesters
#--------------------------------------------------------------

participants_all_trimesters <- crp_trimester %>%
  distinct(
    pregid,
    trimester
  ) %>%
  count(
    pregid,
    name = "n_trimesters"
  ) %>%
  summarise(
    participants_all_three = sum(n_trimesters == 3),
    percent_all_three = round(
      100 * participants_all_three / total_participants,
      1
    )
  )

participants_all_trimesters

#--------------------------------------------------------------
# 4b. Pregnancies with CRP measurements in both T1 and T2
#--------------------------------------------------------------

participants_t1_t2 <- crp_trimester %>%
  filter(trimester %in% c("T1", "T2")) %>%
  distinct(
    pregid,
    trimester
  ) %>%
  count(
    pregid,
    name = "n_trimesters"
  ) %>%
  summarise(
    participants_t1_t2 = sum(n_trimesters == 2),
    
    percent_t1_t2 = round(
      100 * participants_t1_t2 / total_participants,
      1
    )
  )

participants_t1_t2


#--------------------------------------------------------------
# 5. Overall empirical CRP percentiles
#
# These are observation-level percentiles. Every CRP
# measurement contributes one value.
#--------------------------------------------------------------

empirical_centiles_overall <- crp_trimester %>%
  summarise(
    observations = n(),
    participants = n_distinct(pregid),
    
    P5 = unname(
      quantile(
        imputed_crp,
        probs = 0.05,
        na.rm = TRUE
      )
    ),
    
    P50 = unname(
      quantile(
        imputed_crp,
        probs = 0.50,
        na.rm = TRUE
      )
    ),
    
    P95 = unname(
      quantile(
        imputed_crp,
        probs = 0.95,
        na.rm = TRUE
      )
    )
  ) %>%
  mutate(
    across(
      c(P5, P50, P95),
      \(x) round(x, 2)
    )
  )

empirical_centiles_overall


#--------------------------------------------------------------
# 6. Empirical CRP percentiles by trimester
#
# These are calculated directly from all imputed CRP
# observations within each trimester.
#--------------------------------------------------------------

empirical_centiles_trimester <- crp_trimester %>%
  group_by(trimester) %>%
  summarise(
    observations = n(),
    
    participants = n_distinct(pregid),
    
    P5 = unname(
      quantile(
        imputed_crp,
        probs = 0.05,
        na.rm = TRUE
      )
    ),
    
    P50 = unname(
      quantile(
        imputed_crp,
        probs = 0.50,
        na.rm = TRUE
      )
    ),
    
    P95 = unname(
      quantile(
        imputed_crp,
        probs = 0.95,
        na.rm = TRUE
      )
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(P5, P50, P95),
      \(x) round(x, 2)
    )
  ) %>%
  arrange(trimester)

empirical_centiles_trimester


#--------------------------------------------------------------
# 7. Combine trimester and overall percentile results
#--------------------------------------------------------------

empirical_centiles_combined <- bind_rows(
  empirical_centiles_trimester %>%
    mutate(
      Group = as.character(trimester)
    ) %>%
    dplyr::select(
      Group,
      observations,
      P5,
      P50,
      P95
    ),
  
  empirical_centiles_overall %>%
    mutate(
      Group = "Overall"
    ) %>%
    dplyr::select(
      Group,
      observations,
      P5,
      P50,
      P95
    )
)

empirical_centiles_combined


#--------------------------------------------------------------
# 8. Publication-style gt table
#--------------------------------------------------------------

empirical_centiles_gt <- empirical_centiles_combined %>%
  gt() %>%
  cols_label(
    Group = "Group",
    observations = "CRP observations, n",
    P5 = "5th percentile",
    P50 = "50th percentile",
    P95 = "95th percentile"
  ) %>%
  fmt_number(
    columns = c(P5, P50, P95),
    decimals = 2
  ) %>%
  cols_align(
    align = "left",
    columns = Group
  ) %>%
  cols_align(
    align = "center",
    columns = -Group
  ) %>%
  tab_header(
    title = "Empirical CRP percentiles overall and by trimester"
  ) %>%
  tab_source_note(
    source_note = paste0(
      "Percentiles were calculated directly from imputed CRP ",
      "observations. The same pregnancy may contribute observations ",
      "in more than one trimester."
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(
      columns = everything()
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = Group == "Overall"
    )
  ) %>%
  tab_options(
    table.font.size = px(14),
    data_row.padding = px(8),
    column_labels.padding = px(10),
    table.width = pct(100)
  )

empirical_centiles_gt



#--------------------------------------------------------------
# 9. Optional gt table for trimester participation
#--------------------------------------------------------------

participants_by_trimester_gt <- participants_by_trimester %>%
  gt() %>%
  cols_label(
    trimester = "Trimester",
    participants = "Pregnancies, n",
    percent = "Pregnancies represented (%)"
  ) %>%
  fmt_number(
    columns = percent,
    decimals = 1
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_header(
    title = "Pregnancies represented by trimester"
  ) %>%
  tab_source_note(
    source_note = paste0(
      "Percentages may sum to more than 100% because a pregnancy ",
      "may have CRP measurements in multiple trimesters."
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(
      columns = everything()
    )
  ) %>%
  tab_options(
    table.font.size = px(14),
    data_row.padding = px(8),
    column_labels.padding = px(10),
    table.width = pct(100)
  )

participants_by_trimester_gt
