library(dplyr)
library(tibble)

# Representative GA for each trimester
trimester_ga <- tibble(
  Trimester = c("1st trimester", "2nd trimester", "3rd trimester"),
  ga_wks = c(
    10,    # midpoint of 6-14
    21,    # midpoint of 14-28
    31.5   # midpoint of 28-35
  )
)

get_trimester_fpr <- function(centile_df, marker_name) {
  
  trimester_ga %>%
    rowwise() %>%
    mutate(
      Marker = marker_name,
      
      P5 = approx(
        x = centile_df$ga_wks,
        y = centile_df$centile_5,
        xout = ga_wks
      )$y,
      
      P50 = approx(
        x = centile_df$ga_wks,
        y = centile_df$centile_50,
        xout = ga_wks
      )$y,
      
      P95 = approx(
        x = centile_df$ga_wks,
        y = centile_df$centile_95,
        xout = ga_wks
      )$y
    ) %>%
    ungroup() %>%
    select(
      Marker,
      Trimester,
      ga_wks,
      P5,
      P50,
      P95
    )
}
trimester_centiles <- bind_rows(
  get_trimester_fpr(centile_agp_bc, "AGP"),
  get_trimester_fpr(centile_nlr_bc, "NLR"),
  get_trimester_fpr(centile_lmr_bc, "LMR"),
  get_trimester_fpr(centile_plr_bc, "PLR"),
  get_trimester_fpr(centile_sii_bc, "SII")
) %>%
  mutate(
    P5 = round(P5, 2),
    P95 = round(P95, 2)
  )

View(trimester_centiles)

trapz_mean <- function(df, centile, start, end) {
  
  tmp <- df %>%
    filter(ga_wks >= start, ga_wks <= end) %>%
    arrange(ga_wks)
  
  x <- tmp$ga_wks
  y <- tmp[[centile]]
  
  auc <- sum(
    diff(x) * (head(y, -1) + tail(y, -1)) / 2
  )
  
  auc / (end - start)
}
trimester_fpr <- function(df, marker) {
  
  tibble(
    Marker = marker,
    Trimester = c("1st", "2nd", "3rd"),
    
    P5 = c(
      trapz_mean(df, "centile_5", 6, 14),
      trapz_mean(df, "centile_5", 14, 28),
      trapz_mean(df, "centile_5", 28, 35)
    ),
    
    P50 = c(
      trapz_mean(df, "centile_50", 6, 14),
      trapz_mean(df, "centile_50", 14, 28),
      trapz_mean(df, "centile_50", 28, 35)
    ),
    
    P95 = c(
      trapz_mean(df, "centile_95", 6, 14),
      trapz_mean(df, "centile_95", 14, 28),
      trapz_mean(df, "centile_95", 28, 35)
    )
  )
}
trimester_centiles <- bind_rows(
  trimester_fpr(centile_agp_bc, "AGP"),
  trimester_fpr(centile_nlr_bc, "NLR"),
  trimester_fpr(centile_lmr_bc, "LMR"),
  trimester_fpr(centile_plr_bc, "PLR"),
  trimester_fpr(centile_sii_bc, "SII")
) %>%
  mutate(
    across(c(P5, P50, P95), ~ round(.x, 2))
  )

library(gtsummary)
library(gt)

tri_cent <- trimester_centiles %>%
  select(Marker, Trimester, P5, P50, P95) %>%
  gt(groupname_col = "Marker") %>%
  cols_label(
    Trimester = "Trimester",
    P5 = "5th",
    P50 = "50th ",
    P95 = "95th"
  ) %>%
  fmt_number(
    columns = c(P5, P50, P95),
    decimals = 2
  ) %>%
  tab_header(
    title = "Centiles by Trimester")

tri_cent

trimester_centiles

#**************************************************************
# CALCULATE MEAN FITTED CENTILES BY TRIMESTER
#**************************************************************

library(dplyr)
library(tibble)

trapz_mean <- function(df, centile, start, end) {
  
  # Create prediction points including exact boundaries
  x <- sort(unique(c(
    start,
    df$ga_wks[df$ga_wks > start & df$ga_wks < end],
    end
  )))
  
  # Interpolate fitted centile values
  y <- approx(
    x = df$ga_wks,
    y = df[[centile]],
    xout = x
  )$y
  
  # Calculate area under the fitted centile curve
  auc <- sum(
    diff(x) *
      (head(y, -1) + tail(y, -1)) / 2
  )
  
  # Divide by trimester duration
  auc / (end - start)
}

trimester_fpr <- function(df, marker) {
  
  intervals <- tibble(
    Trimester = c(
      "1st trimester",
      "2nd trimester",
      "3rd trimester"
    ),
    start = c(6, 14, 28),
    end = c(14, 28, 35)
  )
  
  intervals %>%
    rowwise() %>%
    mutate(
      Marker = marker,
      
      P5 = trapz_mean(
        df, "centile_5", start, end
      ),
      
      P50 = trapz_mean(
        df, "centile_50", start, end
      ),
      
      P95 = trapz_mean(
        df, "centile_95", start, end
      )
    ) %>%
    ungroup() %>%
    select(Marker, Trimester, P5, P50, P95)
}

#**************************************************************
# COMBINE ALL BIOMARKERS
#**************************************************************

trimester_centiles <- bind_rows(
  
  trimester_fpr(centile_agp_bc, "AGP"),
  trimester_fpr(centile_nlr_bc, "NLR"),
  trimester_fpr(centile_lmr_bc, "LMR"),
  trimester_fpr(centile_plr_bc, "PLR"),
  trimester_fpr(centile_sii_bc, "SII")
  
) %>%
  mutate(
    across(
      c(P5, P50, P95),
      ~ round(.x, 2)
    )
  )

#**************************************************************
# DISPLAY TABLE
#**************************************************************

trimester_centiles %>%
  gt(groupname_col = "Marker") %>%
  cols_label(
    Trimester = "Trimester",
    P5 = "5th centile",
    P50 = "50th centile",
    P95 = "95th centile"
  ) %>%
  fmt_number(
    columns = c(P5, P50, P95),
    decimals = 2
  ) %>%
  tab_header(
    title = "Mean Gestational Age-Specific Centiles by Trimester",
    subtitle = "Box-Cox Fractional Polynomial Regression"
  )










10/22/26
library(dplyr)
library(tibble)

# Representative GA for each trimester
trimester_ga <- tibble(
  Trimester = c("1st trimester", "2nd trimester", "3rd trimester"),
  ga_wks = c(
    10,    # midpoint of 6-14
    21,    # midpoint of 14-28
    31.5   # midpoint of 28-35
  )
)

get_trimester_fpr <- function(centile_df, marker_name) {
  
  trimester_ga %>%
    rowwise() %>%
    mutate(
      Marker = marker_name,
      
      P5 = approx(
        x = centile_df$ga_wks,
        y = centile_df$centile_5,
        xout = ga_wks
      )$y,
      
      P50 = approx(
        x = centile_df$ga_wks,
        y = centile_df$centile_50,
        xout = ga_wks
      )$y,
      
      P95 = approx(
        x = centile_df$ga_wks,
        y = centile_df$centile_95,
        xout = ga_wks
      )$y
    ) %>%
    ungroup() %>%
    select(
      Marker,
      Trimester,
      ga_wks,
      P5,
      P50,
      P95
    )
}
trimester_centiles <- bind_rows(
  get_trimester_fpr(centile_agp_bc, "AGP"),
  get_trimester_fpr(centile_nlr_bc, "NLR"),
  get_trimester_fpr(centile_lmr_bc, "LMR"),
  get_trimester_fpr(centile_plr_bc, "PLR"),
  get_trimester_fpr(centile_sii_bc, "SII")
) %>%
  mutate(
    P5 = round(P5, 2),
    P95 = round(P95, 2)
  )

View(trimester_centiles)

trapz_mean <- function(df, centile, start, end) {
  
  tmp <- df %>%
    filter(ga_wks >= start, ga_wks <= end) %>%
    arrange(ga_wks)
  
  x <- tmp$ga_wks
  y <- tmp[[centile]]
  
  auc <- sum(
    diff(x) * (head(y, -1) + tail(y, -1)) / 2
  )
  
  auc / (end - start)
}
trimester_fpr <- function(df, marker) {
  
  tibble(
    Marker = marker,
    Trimester = c("1st", "2nd", "3rd"),
    
    P5 = c(
      trapz_mean(df, "centile_5", 6, 14),
      trapz_mean(df, "centile_5", 14, 28),
      trapz_mean(df, "centile_5", 28, 35)
    ),
    
    P50 = c(
      trapz_mean(df, "centile_50", 6, 14),
      trapz_mean(df, "centile_50", 14, 28),
      trapz_mean(df, "centile_50", 28, 35)
    ),
    
    P95 = c(
      trapz_mean(df, "centile_95", 6, 14),
      trapz_mean(df, "centile_95", 14, 28),
      trapz_mean(df, "centile_95", 28, 35)
    )
  )
}
trimester_centiles <- bind_rows(
  trimester_fpr(centile_agp_bc, "AGP"),
  trimester_fpr(centile_nlr_bc, "NLR"),
  trimester_fpr(centile_lmr_bc, "LMR"),
  trimester_fpr(centile_plr_bc, "PLR"),
  trimester_fpr(centile_sii_bc, "SII")
) %>%
  mutate(
    across(c(P5, P50, P95), ~ round(.x, 2))
  )

library(gtsummary)
library(gt)

tri_cent <- trimester_centiles %>%
  select(Marker, Trimester, P5, P50, P95) %>%
  gt(groupname_col = "Marker") %>%
  cols_label(
    Trimester = "Trimester",
    P5 = "5th",
    P50 = "50th ",
    P95 = "95th"
  ) %>%
  fmt_number(
    columns = c(P5, P50, P95),
    decimals = 2
  ) %>%
  tab_header(
    title = "Centiles by Trimester")

tri_cent

trimester_centiles

#**************************************************************
# CALCULATE MEAN FITTED CENTILES BY TRIMESTER
#**************************************************************

library(dplyr)
library(tibble)

trapz_mean <- function(df, centile, start, end) {
  
  # Create prediction points including exact boundaries
  x <- sort(unique(c(
    start,
    df$ga_wks[df$ga_wks > start & df$ga_wks < end],
    end
  )))
  
  # Interpolate fitted centile values
  y <- approx(
    x = df$ga_wks,
    y = df[[centile]],
    xout = x
  )$y
  
  # Calculate area under the fitted centile curve
  auc <- sum(
    diff(x) *
      (head(y, -1) + tail(y, -1)) / 2
  )
  
  # Divide by trimester duration
  auc / (end - start)
}

trimester_fpr <- function(df, marker) {
  
  intervals <- tibble(
    Trimester = c(
      "1st trimester",
      "2nd trimester",
      "3rd trimester"
    ),
    start = c(6, 14, 28),
    end = c(14, 28, 35)
  )
  
  intervals %>%
    rowwise() %>%
    mutate(
      Marker = marker,
      
      P5 = trapz_mean(
        df, "centile_5", start, end
      ),
      
      P50 = trapz_mean(
        df, "centile_50", start, end
      ),
      
      P95 = trapz_mean(
        df, "centile_95", start, end
      )
    ) %>%
    ungroup() %>%
    select(Marker, Trimester, P5, P50, P95)
}

#**************************************************************
# COMBINE ALL BIOMARKERS
#**************************************************************

trimester_centiles <- bind_rows(
  
  trimester_fpr(centile_agp_bc, "AGP"),
  trimester_fpr(centile_nlr_bc, "NLR"),
  trimester_fpr(centile_lmr_bc, "LMR"),
  trimester_fpr(centile_plr_bc, "PLR"),
  trimester_fpr(centile_sii_bc, "SII")
  
) %>%
  mutate(
    across(
      c(P5, P50, P95),
      ~ round(.x, 2)
    )
  )

#**************************************************************
# DISPLAY TABLE
#**************************************************************

trimester_centiles %>%
  gt(groupname_col = "Marker") %>%
  cols_label(
    Trimester = "Trimester",
    P5 = "5th centile",
    P50 = "50th centile",
    P95 = "95th centile"
  ) %>%
  fmt_number(
    columns = c(P5, P50, P95),
    decimals = 2
  ) %>%
  tab_header(
    title = "Mean Gestational Age-Specific Centiles by Trimester",
    subtitle = "Box-Cox Fractional Polynomial Regression"
  )
