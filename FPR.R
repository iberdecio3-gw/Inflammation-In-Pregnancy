library(dplyr)
library(tidyr)
library(ggplot2)
library(gamlss)
library(gamlss.dist)
library(MASS)
library(patchwork)

# Example variables -- replace/add your actual CBC-derived ratios
markers <- c(
  "agp",
  "nlr",
  "plr",
  "lmr",
  "sii"
)

missing_codes <- c(-5, -6, -7, -9)

df_bio_long_hc <- df_bio_long_hc %>%
  mutate(
    across(
      c(
        wbc,
        plt,
        m08_cbc_neu_fcc_lborres,
        m08_cbc_lymph_fcc_lborres,
        m08_cbc_mono_fcc_lborres
      ),
      ~ if_else(.x %in% missing_codes, NA_real_, .x)
    )
  )

df_bio_long_hc <- df_bio_long_hc %>%
  mutate(
    nlr = m08_cbc_neu_fcc_lborres /
      m08_cbc_lymph_fcc_lborres,
    
    lmr = m08_cbc_lymph_fcc_lborres /
      m08_cbc_mono_fcc_lborres,
    
    plr = plt /
      m08_cbc_lymph_fcc_lborres,
    
    sii = plt *
      m08_cbc_neu_fcc_lborres /
      m08_cbc_lymph_fcc_lborres
  )



#assign Inf values as non-usuable (denominator 0)
df_bio_long_hc <- df_bio_long_hc %>%
  mutate(
    nlr_analysis = if_else(is.finite(nlr) & nlr > 0, nlr, NA_real_),
    lmr_analysis = if_else(is.finite(lmr) & lmr > 0, lmr, NA_real_),
    plr_analysis = if_else(is.finite(plr) & plr > 0, plr, NA_real_),
    sii_analysis = if_else(is.finite(sii) & sii > 0, sii, NA_real_)
  )

#------------------------------------------------------------
# Basic distribution summary
#------------------------------------------------------------

distribution_summary <- function(data, variable) {
  
  x <- data[[variable]]
  
  tibble(
    Marker = variable,
    N = sum(!is.na(x)),
    Missing = sum(is.na(x)),
    Zero = sum(x == 0, na.rm = TRUE),
    Negative = sum(x < 0, na.rm = TRUE),
    Min = min(x, na.rm = TRUE),
    P05 = quantile(x, 0.05, na.rm = TRUE),
    P25 = quantile(x, 0.25, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Mean = mean(x, na.rm = TRUE),
    P75 = quantile(x, 0.75, na.rm = TRUE),
    P95 = quantile(x, 0.95, na.rm = TRUE),
    Max = max(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE)
  )
}

markers_analysis <- c(
  "agp",
  "nlr_analysis",
  "lmr_analysis",
  "plr_analysis",
  "sii_analysis"
)

marker_summary <- bind_rows(
  lapply(
    markers_analysis,
    \(v) distribution_summary(df_bio_long_hc, v)
  )
)

View(marker_summary)

##############################

# Plot distributions

############################


plot_marker_distribution <- function(data, variable) {
  
  df <- data %>%
    filter(
      !is.na(ga_wks),
      !is.na(.data[[variable]])
    )
  
  p1 <- ggplot(df, aes(x = .data[[variable]])) +
    geom_histogram(bins = 40) +
    labs(
      title = paste(variable, "- Distribution"),
      x = variable,
      y = "Count"
    ) +
    theme_bw()
  
  
  p2 <- ggplot(df, aes(y = .data[[variable]])) +
    geom_boxplot() +
    labs(
      title = paste(variable, "- Boxplot"),
      y = variable
    ) +
    theme_bw()
  
  
  p3 <- ggplot(df, aes(sample = .data[[variable]])) +
    stat_qq() +
    stat_qq_line() +
    labs(
      title = paste(variable, "- Normal Q-Q")
    ) +
    theme_bw()
  
  
  p4 <- ggplot(df, aes(x = ga_wks, y = .data[[variable]])) +
    geom_point(alpha = 0.15, size = 0.5) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      title = paste(variable, "vs Gestational Age"),
      x = "Gestational age (weeks)",
      y = variable
    ) +
    theme_bw()
  
  (p1 | p2) / (p3 | p4)
}

plot_marker_distribution(df_bio_long_hc, "agp")

plot_marker_distribution(df_bio_long_hc, "nlr_analysis")

plot_marker_distribution(df_bio_long_hc, "lmr_analysis")

plot_marker_distribution(df_bio_long_hc, "plr_analysis")

plot_marker_distribution(df_bio_long_hc, "sii_analysis")

normality_results <- bind_rows(
  
  shapiro_by_site(
    df_fpr_agp, "agp", "AGP"
  ),
  
  shapiro_by_site(
    df_fpr_nlr, "nlr_analysis", "NLR"
  ),
  
  shapiro_by_site(
    df_fpr_lmr, "lmr_analysis", "LMR"
  ),
  
  shapiro_by_site(
    df_fpr_plr, "plr_analysis", "PLR"
  ),
  
  shapiro_by_site(
    df_fpr_sii, "sii_analysis", "SII"
  )
)

write.csv(normality_results, "WS_normality_Secondary.csv")


###### ----- Test Different Distributions - AGP ----- #####

#**************************************************************
# AGP - EXCLUDE 128 mg/L OUTLIER
#**************************************************************

df_fpr_agp <- df_bio_long_hc %>%
  dplyr::select(
    site,
    ga_wks,
    agp
  ) %>%
  filter(
    !is.na(ga_wks),
    !is.na(agp),
    is.finite(ga_wks),
    is.finite(agp),
    agp > 0,
    agp <= 120
  )

summary(df_fpr_agp$agp)
nrow(df_fpr_agp)

#agp_gamma <- gamlss(
#  agp ~ fp(ga_wks, 2),
#  sigma.fo = ~ fp(ga_wks, 2),
#  family = GA(),
#  data = df_fpr_agp,
#  control = con
#)

#agp_lognormal <- gamlss(
#  agp ~ fp(ga_wks, 2),
#  sigma.fo = ~ fp(ga_wks, 2),
#  family = LOGNO(),
#  data = df_fpr_agp,
#  control = con
#)

#agp_weibull <- gamlss(
#  agp ~ fp(ga_wks, 2),
#  sigma.fo = ~ fp(ga_wks, 2),
#  family = WEI(),
#  data = df_fpr_agp,
#  control = con
#)

#rm(agp_weibull)

#agp_weibull <- gamlss(
#  agp ~ fp(ga_wks, 2),
#  sigma.fo = ~ fp(ga_wks, 2),
#  family = WEI(),
#  data = df_fpr_agp,
#  control = gamlss.control(
#    trace = TRUE,
#    n.cyc = 100
#  )
#)

#agp_invgauss <- gamlss(
#  agp ~ fp(ga_wks, 2),
#  sigma.fo = ~ fp(ga_wks, 2),
#  family = IG(),
#  data = df_fpr_agp,
#  control = con
#)

###### ----- Test Different Distributions - NLR ----- #####


df_fpr_nlr <- df_bio_long_hc %>%
  dplyr::select(
    site,
    ga_wks,
    nlr_analysis
  ) %>%
  filter(
    !is.na(ga_wks),
    !is.na(nlr_analysis),
    is.finite(ga_wks),
    is.finite(nlr_analysis),
    nlr_analysis > 0
  )

nlr_gamma <- gamlss(
  nlr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = GA(),
  data = df_fpr_nlr,
  control = con
)

nlr_lognormal <- gamlss(
  nlr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = LOGNO(),
  data = df_fpr_nlr,
  control = con
)

nlr_weibull <- gamlss(
  nlr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = WEI(),
  data = df_fpr_agp,
  control = con
)

nlr_invgauss <- gamlss(
  nlr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = IG(),
  data = df_fpr_nlr,
  control = con
)


###### ----- Test Different Distributions - LMR ----- #####

df_fpr_lmr <- df_bio_long_hc %>%
  dplyr::select(
    site,
    ga_wks,
    lmr_analysis
  ) %>%
  filter(
    !is.na(ga_wks),
    !is.na(lmr_analysis),
    is.finite(ga_wks),
    is.finite(lmr_analysis),
    lmr_analysis > 0
  )

lmr_gamma <- gamlss(
  lmr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = GA(),
  data = df_fpr_lmr,
  control = con
)

lmr_lognormal <- gamlss(
  lmr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = LOGNO(),
  data = df_fpr_lmr,
  control = con
)

lmr_weibull <- gamlss(
  lmr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = WEI(),
  data = df_fpr_lmr,
  control = con
)

lmr_invgauss <- gamlss(
  lmr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = IG(),
  data = df_fpr_lmr,
  control = con
)

###### ----- Test Different Distributions - PLR ----- #####
df_fpr_plr <- df_bio_long_hc %>%
  dplyr::select(
    site,
    ga_wks,
    plr_analysis
  ) %>%
  filter(
    !is.na(ga_wks),
    !is.na(plr_analysis),
    is.finite(ga_wks),
    is.finite(plr_analysis),
    plr_analysis > 0
  )

#plr_gamma <- gamlss(
  plr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = GA(),
  data = df_fpr_plr,
  control = con
)

plr_lognormal <- gamlss(
  plr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = LOGNO(),
  data = df_fpr_plr,
  control = con
)

plr_weibull <- gamlss(
  plr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = WEI(),
  data = df_fpr_plr,
  control = con
)

plr_invgauss <- gamlss(
  plr_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = IG(),
  data = df_fpr_plr,
  control = con
)





###### ----- Test Different Distributions - SII ----- #####

df_fpr_sii <- df_bio_long_hc %>%
  dplyr::select(
    site,
    ga_wks,
    sii_analysis
  ) %>%
  filter(
    !is.na(ga_wks),
    !is.na(sii_analysis),
    is.finite(ga_wks),
    is.finite(sii_analysis),
    sii_analysis > 0
  )

sii_gamma <- gamlss(
  sii_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = GA(),
  data = df_fpr_sii,
  control = con
)

sii_lognormal <- gamlss(
  sii_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = LOGNO(),
  data = df_fpr_sii,
  control = con
)

sii_weibull <- gamlss(
  sii_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = WEI(),
  data = df_fpr_sii,
  control = con
)

sii_invgauss <- gamlss(
  sii_analysis ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = IG(),
  data = df_fpr_sii,
  control = con
)



################## Box cox Transformations #############

inverse_boxcox <- function(y, lambda) {
  
  if (abs(lambda) < 1e-6) {
    
    exp(y)
    
  } else {
    
    (lambda * y + 1)^(1 / lambda)
    
  }
}

#**************************************************************
# AGP BOX-COX TRANSFORMATION
#**************************************************************

bc_agp <- MASS::boxcox(
  lm(agp ~ ga_wks, data = df_fpr_agp),
  plotit = FALSE
)

lambda_agp <- bc_agp$x[which.max(bc_agp$y)]

lambda_agp


if (abs(lambda_agp) < 1e-6) {
  
  df_fpr_agp <- df_fpr_agp %>%
    mutate(
      agp_bc = log(agp)
    )
  
} else {
  
  df_fpr_agp <- df_fpr_agp %>%
    mutate(
      agp_bc = (agp^lambda_agp - 1) / lambda_agp
    )
}


agp_boxcox <- gamlss(
  agp_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_agp,
  control = con
)

agp_boxcox_pred_bc <- centiles.pred(
  agp_boxcox,
  data = df_fpr_agp,
  xname = "ga_wks",
  xvalues = agp_ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()

names(agp_boxcox_pred_bc) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)

inverse_boxcox <- function(y, lambda) {
  
  if (abs(lambda) < 1e-6) {
    
    exp(y)
    
  } else {
    
    (lambda * y + 1)^(1 / lambda)
    
  }
}

agp_boxcox_pred <- agp_boxcox_pred_bc %>%
  mutate(
    centile_5 = inverse_boxcox(
      centile_5,
      lambda_agp
    ),
    
    centile_50 = inverse_boxcox(
      centile_50,
      lambda_agp
    ),
    
    centile_95 = inverse_boxcox(
      centile_95,
      lambda_agp
    )
  )

lambda_agp

head(agp_boxcox_pred)

summary(agp_boxcox_pred)

#**************************************************************
# NLR BOX-COX TRANSFORMATION
#**************************************************************

bc_nlr <- MASS::boxcox(
  lm(nlr_analysis ~ ga_wks, data = df_fpr_nlr),
  plotit = FALSE
)

lambda_nlr <- bc_nlr$x[which.max(bc_nlr$y)]

lambda_nlr


if (abs(lambda_nlr) < 1e-6) {
  
  df_fpr_nlr <- df_fpr_nlr %>%
    mutate(
      nlr_bc = log(nlr_analysis)
    )
  
} else {
  
  df_fpr_nlr <- df_fpr_nlr %>%
    mutate(
      nlr_bc = (nlr_analysis^lambda_nlr - 1) / lambda_nlr
    )
}


nlr_boxcox <- gamlss(
  nlr_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_nlr,
  control = con
)

#**************************************************************
# LMR BOX-COX TRANSFORMATION
#**************************************************************

bc_lmr <- MASS::boxcox(
  lm(lmr_analysis ~ ga_wks, data = df_fpr_lmr),
  plotit = FALSE
)

lambda_lmr <- bc_lmr$x[which.max(bc_lmr$y)]

lambda_lmr


if (abs(lambda_lmr) < 1e-6) {
  
  df_fpr_lmr <- df_fpr_lmr %>%
    mutate(
      lmr_bc = log(lmr_analysis)
    )
  
} else {
  
  df_fpr_lmr <- df_fpr_lmr %>%
    mutate(
      lmr_bc = (lmr_analysis^lambda_lmr - 1) / lambda_lmr
    )
}


lmr_boxcox <- gamlss(
  lmr_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_lmr,
  control = con
)

#**************************************************************
# PLR BOX-COX TRANSFORMATION
#**************************************************************

bc_plr <- MASS::boxcox(
  lm(plr_analysis ~ ga_wks, data = df_fpr_plr),
  plotit = FALSE
)

lambda_plr <- bc_plr$x[which.max(bc_plr$y)]

lambda_plr


if (abs(lambda_plr) < 1e-6) {
  
  df_fpr_plr <- df_fpr_plr %>%
    mutate(
      plr_bc = log(plr_analysis)
    )
  
} else {
  
  df_fpr_plr <- df_fpr_plr %>%
    mutate(
      plr_bc = (plr_analysis^lambda_plr - 1) / lambda_plr
    )
}


plr_boxcox <- gamlss(
  plr_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_plr,
  control = con
)

#**************************************************************
# SII BOX-COX TRANSFORMATION
#**************************************************************

bc_sii <- MASS::boxcox(
  lm(sii_analysis ~ ga_wks, data = df_fpr_sii),
  plotit = FALSE
)

lambda_sii <- bc_sii$x[which.max(bc_sii$y)]

lambda_sii


if (abs(lambda_sii) < 1e-6) {
  
  df_fpr_sii <- df_fpr_sii %>%
    mutate(
      sii_bc = log(sii_analysis)
    )
  
} else {
  
  df_fpr_sii <- df_fpr_sii %>%
    mutate(
      sii_bc = (sii_analysis^lambda_sii - 1) / lambda_sii
    )
}


sii_boxcox <- gamlss(
  sii_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_sii,
  control = con
)

tibble(
  Marker = c("AGP", "NLR", "LMR", "PLR", "SII"),
  Lambda = c(
    lambda_agp,
    lambda_nlr,
    lambda_lmr,
    lambda_plr,
    lambda_sii
  )
)

######### Compare

GAIC(
  agp_gamma,
  agp_lognormal,
  #agp_weibull,
  agp_invgauss
)

GAIC(
  nlr_gamma,
  nlr_lognormal,
  nlr_weibull,
  nlr_invgauss
)

GAIC(
  lmr_gamma,
  lmr_lognormal,
  lmr_weibull,
  lmr_invgauss
)

GAIC(
  plr_gamma,
  plr_lognormal,
  #plr_weibull,
  plr_invgauss
)

GAIC(
  sii_gamma,
  sii_lognormal,
  sii_weibull,
  sii_invgauss
)

################ plot best fits #########

#**************************************************************
# AGP LOG-NORMAL CENTILE CURVES
#**************************************************************
{
  ga_values <- seq(6, 40, by = 0.1)
  
  centile_agp_lognormal <- centiles.pred(
    agp_lognormal,
    data = df_fpr_agp,
    xname = "ga_wks",
    xvalues = ga_values,
    type = "centiles",
    cent = c(5, 50, 95)
  ) %>%
    as.data.frame()
  
  names(centile_agp_lognormal) <- c(
    "ga_wks",
    "centile_5",
    "centile_50",
    "centile_95"
  )
  
  plot(
    df_fpr_agp$ga_wks,
    df_fpr_agp$agp,
    xlab = "Gestational Age (weeks)",
    ylab = "AGP",
    main = "AGP Log-Normal Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
    bg = "transparent",
    xlim = c(6, 40),
    ylim = c(0, 2),
    pch = 15,
    cex = 0.2,
    col = "grey86"
  )
  
  lines(
    centile_agp_lognormal$ga_wks,
    centile_agp_lognormal$centile_5,
    col = "red",
    lwd = 1.5
  )
  
  lines(
    centile_agp_lognormal$ga_wks,
    centile_agp_lognormal$centile_50,
    col = "green",
    lwd = 1.5
  )
  
  lines(
    centile_agp_lognormal$ga_wks,
    centile_agp_lognormal$centile_95,
    col = "blue",
    lwd = 1.5
  )
  
  legend(
    "topright",
    legend = c("5th", "50th", "95th"),
    cex = 0.7,
    col = c("red", "green", "blue"),
    lty = 1,
    lwd = 1.5
  )
  
  centile_agp_gamma <- centiles.pred(
    agp_gamma,
    data = df_fpr_agp,
    xname = "ga_wks",
    xvalues = seq(6, 40, 0.1),
    type = "centiles",
    cent = c(5, 50, 95)
  ) %>%
    as.data.frame()
  
}

#**************************************************************
# NLR - LOG-NORMAL
#**************************************************************
{
  centile_nlr_lognormal <- centiles.pred(
    nlr_lognormal,
    data = df_fpr_nlr,
    xname = "ga_wks",
    xvalues = seq(6, 40, 0.1),
    type = "centiles",
    cent = c(5, 50, 95)
  ) %>%
    as.data.frame()
  
  names(centile_nlr_lognormal) <- c(
    "ga_wks",
    "centile_5",
    "centile_50",
    "centile_95"
  )
  
  plot(
    df_fpr_nlr$ga_wks,
    df_fpr_nlr$nlr_analysis,
    xlab = "Gestational Age (weeks)",
    ylab = "NLR",
    main = "NLR Log-Normal Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
    xlim = c(6, 40),
    ylim = c(0, 10),
    pch = 15,
    cex = 0.2,
    col = "grey86"
  )
  
  lines(
    centile_nlr_lognormal$ga_wks,
    centile_nlr_lognormal$centile_5,
    col = "red",
    lwd = 1.5
  )
  
  lines(
    centile_nlr_lognormal$ga_wks,
    centile_nlr_lognormal$centile_50,
    col = "green",
    lwd = 1.5
  )
  
  lines(
    centile_nlr_lognormal$ga_wks,
    centile_nlr_lognormal$centile_95,
    col = "blue",
    lwd = 1.5
  )
  
  legend(
    "topright",
    legend = c("5th", "50th", "95th"),
    col = c("red", "green", "blue"),
    lty = 1,
    lwd = 1.5,
    cex = 0.7
  )
}

#**************************************************************
# LMR - LOG-NORMAL
#**************************************************************
{
  centile_lmr_lognormal <- centiles.pred(
    lmr_lognormal,
    data = df_fpr_lmr,
    xname = "ga_wks",
    xvalues = seq(6, 40, 0.1),
    type = "centiles",
    cent = c(5, 50, 95)
  ) %>%
    as.data.frame()
  
  names(centile_lmr_lognormal) <- c(
    "ga_wks",
    "centile_5",
    "centile_50",
    "centile_95"
  )
  
  plot(
    df_fpr_lmr$ga_wks,
    df_fpr_lmr$lmr_analysis,
    xlab = "Gestational Age (weeks)",
    ylab = "LMR",
    main = "LMR Log-Normal Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
    xlim = c(6, 40),
    ylim = c(0, 10),
    pch = 15,
    cex = 0.2,
    col = "grey86"
  )
  
  lines(
    centile_lmr_lognormal$ga_wks,
    centile_lmr_lognormal$centile_5,
    col = "red",
    lwd = 1.5
  )
  
  lines(
    centile_lmr_lognormal$ga_wks,
    centile_lmr_lognormal$centile_50,
    col = "green",
    lwd = 1.5
  )
  
  lines(
    centile_lmr_lognormal$ga_wks,
    centile_lmr_lognormal$centile_95,
    col = "blue",
    lwd = 1.5
  )
  
  legend(
    "topright",
    legend = c("5th", "50th", "95th"),
    col = c("red", "green", "blue"),
    lty = 1,
    lwd = 1.5,
    cex = 0.7
  )
}

#**************************************************************
# PLR - GAMMA
#**************************************************************
{
  centile_plr_gamma <- centiles.pred(
    plr_gamma,
    data = df_fpr_plr,
    xname = "ga_wks",
    xvalues = seq(6, 40, 0.1),
    type = "centiles",
    cent = c(5, 50, 95)
  ) %>%
    as.data.frame()
  
  names(centile_plr_gamma) <- c(
    "ga_wks",
    "centile_5",
    "centile_50",
    "centile_95"
  )
  
  plot(
    df_fpr_plr$ga_wks,
    df_fpr_plr$plr_analysis,
    xlab = "Gestational Age (weeks)",
    ylab = "PLR",
    main = "PLR Gamma Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
    xlim = c(6, 40),
    ylim = c(0, 400),
    pch = 15,
    cex = 0.2,
    col = "grey86"
  )
  
  lines(
    centile_plr_gamma$ga_wks,
    centile_plr_gamma$centile_5,
    col = "red",
    lwd = 1.5
  )
  
  lines(
    centile_plr_gamma$ga_wks,
    centile_plr_gamma$centile_50,
    col = "green",
    lwd = 1.5
  )
  
  lines(
    centile_plr_gamma$ga_wks,
    centile_plr_gamma$centile_95,
    col = "blue",
    lwd = 1.5
  )
  
  legend(
    "topright",
    legend = c("5th", "50th", "95th"),
    col = c("red", "green", "blue"),
    lty = 1,
    lwd = 1.5,
    cex = 0.7
  )
}

#**************************************************************
# SII - GAMMA
#**************************************************************

centile_sii_gamma <- centiles.pred(
  sii_gamma,
  data = df_fpr_sii,
  xname = "ga_wks",
  xvalues = seq(6, 40, 0.1),
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()

names(centile_sii_gamma) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)

plot(
  df_fpr_sii$ga_wks,
  df_fpr_sii$sii_analysis,
  xlab = "Gestational Age (weeks)",
  ylab = "SII",
  main = "SII Gamma Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  xlim = c(6, 40),
  ylim = c(0, 2500),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(
  centile_sii_gamma$ga_wks,
  centile_sii_gamma$centile_5,
  col = "red",
  lwd = 1.5
)

lines(
  centile_sii_gamma$ga_wks,
  centile_sii_gamma$centile_50,
  col = "green",
  lwd = 1.5
)

lines(
  centile_sii_gamma$ga_wks,
  centile_sii_gamma$centile_95,
  col = "blue",
  lwd = 1.5
)

legend(
  "topright",
  legend = c("5th", "50th", "95th"),
  col = c("red", "green", "blue"),
  lty = 1,
  lwd = 1.5,
  cex = 0.7
)


######### ----- "Winner" vs Box-Cox Transformation ---- ######

library(dplyr)
library(tidyr)
library(gamlss)

#**************************************************************
# Helper: inverse Box-Cox
#**************************************************************

inverse_boxcox <- function(y, lambda) {
  
  if (abs(lambda) < 1e-6) {
    
    exp(y)
    
  } else {
    
    (lambda * y + 1)^(1 / lambda)
    
  }
}


#**************************************************************
# Helper: empirical centile coverage
#**************************************************************

coverage_native <- function(model, data, outcome, label) {
  
  ga_values <- sort(unique(data$ga_wks))
  
  pred <- centiles.pred(
    model,
    data = data,
    xname = "ga_wks",
    xvalues = ga_values,
    type = "centiles",
    cent = c(5, 50, 95)
  ) %>%
    as.data.frame()
  
  names(pred) <- c(
    "ga_wks",
    "centile_5",
    "centile_50",
    "centile_95"
  )
  
  check <- data %>%
    left_join(
      pred,
      by = "ga_wks"
    )
  
  x <- check[[outcome]]
  
  tibble(
    Distribution = label,
    
    P5 = mean(
      x <= check$centile_5,
      na.rm = TRUE
    ) * 100,
    
    P50 = mean(
      x <= check$centile_50,
      na.rm = TRUE
    ) * 100,
    
    P95 = mean(
      x <= check$centile_95,
      na.rm = TRUE
    ) * 100
  ) %>%
    mutate(
      MACE = (
        abs(P5 - 5) +
          abs(P50 - 50) +
          abs(P95 - 95)
      ) / 3
    )
}


#**************************************************************
# Helper: Box-Cox empirical centile coverage
#**************************************************************

coverage_boxcox <- function(
    model,
    data,
    outcome,
    lambda,
    label = "Box-Cox"
) {
  
  ga_values <- sort(unique(data$ga_wks))
  
  pred_bc <- centiles.pred(
    model,
    data = data,
    xname = "ga_wks",
    xvalues = ga_values,
    type = "centiles",
    cent = c(5, 50, 95)
  ) %>%
    as.data.frame()
  
  names(pred_bc) <- c(
    "ga_wks",
    "centile_5",
    "centile_50",
    "centile_95"
  )
  
  pred <- pred_bc %>%
    mutate(
      centile_5 = inverse_boxcox(
        centile_5,
        lambda
      ),
      
      centile_50 = inverse_boxcox(
        centile_50,
        lambda
      ),
      
      centile_95 = inverse_boxcox(
        centile_95,
        lambda
      )
    )
  
  check <- data %>%
    left_join(
      pred,
      by = "ga_wks"
    )
  
  x <- check[[outcome]]
  
  tibble(
    Distribution = label,
    
    P5 = mean(
      x <= check$centile_5,
      na.rm = TRUE
    ) * 100,
    
    P50 = mean(
      x <= check$centile_50,
      na.rm = TRUE
    ) * 100,
    
    P95 = mean(
      x <= check$centile_95,
      na.rm = TRUE
    ) * 100
  ) %>%
    mutate(
      MACE = (
        abs(P5 - 5) +
          abs(P50 - 50) +
          abs(P95 - 95)
      ) / 3
    )
}


#**************************************************************
# AGP
# Best native = Log-normal
#**************************************************************

coverage_agp_native <- coverage_native(
  model = agp_lognormal,
  data = df_fpr_agp,
  outcome = "agp",
  label = "Log-normal"
)

coverage_agp_bc <- coverage_boxcox(
  model = agp_boxcox,
  data = df_fpr_agp,
  outcome = "agp",
  lambda = lambda_agp
)

coverage_agp <- bind_rows(
  coverage_agp_native,
  coverage_agp_bc
) %>%
  mutate(
    Marker = "AGP",
    .before = 1
  )

coverage_agp


#**************************************************************
# NLR
# Best native = Log-normal
#**************************************************************

coverage_nlr_native <- coverage_native(
  model = nlr_lognormal,
  data = df_fpr_nlr,
  outcome = "nlr_analysis",
  label = "Log-normal"
)

coverage_nlr_bc <- coverage_boxcox(
  model = nlr_boxcox,
  data = df_fpr_nlr,
  outcome = "nlr_analysis",
  lambda = lambda_nlr
)

coverage_nlr <- bind_rows(
  coverage_nlr_native,
  coverage_nlr_bc
) %>%
  mutate(
    Marker = "NLR",
    .before = 1
  )

coverage_nlr


#**************************************************************
# LMR
# Best native = Log-normal
#**************************************************************

coverage_lmr_native <- coverage_native(
  model = lmr_lognormal,
  data = df_fpr_lmr,
  outcome = "lmr_analysis",
  label = "Log-normal"
)

coverage_lmr_bc <- coverage_boxcox(
  model = lmr_boxcox,
  data = df_fpr_lmr,
  outcome = "lmr_analysis",
  lambda = lambda_lmr
)

coverage_lmr <- bind_rows(
  coverage_lmr_native,
  coverage_lmr_bc
) %>%
  mutate(
    Marker = "LMR",
    .before = 1
  )

coverage_lmr


#**************************************************************
# PLR
# Best native = Gamma
#**************************************************************

coverage_plr_native <- coverage_native(
  model = plr_gamma,
  data = df_fpr_plr,
  outcome = "plr_analysis",
  label = "Gamma"
)

coverage_plr_bc <- coverage_boxcox(
  model = plr_boxcox,
  data = df_fpr_plr,
  outcome = "plr_analysis",
  lambda = lambda_plr
)

coverage_plr <- bind_rows(
  coverage_plr_native,
  coverage_plr_bc
) %>%
  mutate(
    Marker = "PLR",
    .before = 1
  )

coverage_plr


#**************************************************************
# SII
# Best native = Gamma
#**************************************************************

coverage_sii_native <- coverage_native(
  model = sii_gamma,
  data = df_fpr_sii,
  outcome = "sii_analysis",
  label = "Gamma"
)

coverage_sii_bc <- coverage_boxcox(
  model = sii_boxcox,
  data = df_fpr_sii,
  outcome = "sii_analysis",
  lambda = lambda_sii
)

coverage_sii <- bind_rows(
  coverage_sii_native,
  coverage_sii_bc
) %>%
  mutate(
    Marker = "SII",
    .before = 1
  )

coverage_sii


#**************************************************************
# COMBINE ALL RESULTS
#**************************************************************

coverage_all <- bind_rows(
  coverage_agp,
  coverage_nlr,
  coverage_lmr,
  coverage_plr,
  coverage_sii
) %>%
  arrange(
    Marker,
    MACE
  )

coverage_all


#**************************************************************
# OPTIONAL: nicely rounded table
#**************************************************************

coverage_all_rounded <- coverage_all %>%
  mutate(
    across(
      c(P5, P50, P95, MACE),
      ~ round(.x, 2)
    )
  )

coverage_all_rounded


#**************************************************************
# OPTIONAL: GT TABLE
#**************************************************************

library(gt)

coverage_all_rounded %>%
  gt(
    groupname_col = "Marker"
  ) %>%
  cols_label(
    Distribution = "Model",
    P5 = "P5 (%)",
    P50 = "P50 (%)",
    P95 = "P95 (%)",
    MACE = "MACE"
  ) %>%
  tab_spanner(
    label = "Empirical centile coverage",
    columns = c(P5, P50, P95)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  )


############ ----- Plot Box Cox --------- #####

#**************************************************************
# AGP BOX-COX CENTILE CURVES
#**************************************************************

ga_values <- seq(6, 40, by = 0.1)


# Predict centiles on Box-Cox transformed scale
centile_agp_bc <- centiles.pred(
  agp_boxcox,
  data = df_fpr_agp,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()


# Correct object name here
names(centile_agp_bc) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)


#**************************************************************
# Back-transform centiles to original AGP scale
#**************************************************************

centile_agp_bc <- centile_agp_bc %>%
  mutate(
    centile_5 = inverse_boxcox(
      centile_5,
      lambda_agp
    ),
    
    centile_50 = inverse_boxcox(
      centile_50,
      lambda_agp
    ),
    
    centile_95 = inverse_boxcox(
      centile_95,
      lambda_agp
    )
  )


# Check the back-transformed values
summary(centile_agp_bc)


#**************************************************************
# Plot
#**************************************************************

plot(
  df_fpr_agp$ga_wks,
  df_fpr_agp$agp,
  xlab = "Gestational Age (weeks)",
  ylab = "AGP",
  main = "AGP Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  bg = "transparent",
  xlim = c(6, 40),
  ylim = c(0, 2),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)


# 5th centile
lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_5,
  col = "red",
  lwd = 1.5
)


# 50th centile
lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_50,
  col = "green",
  lwd = 1.5
)


# 95th centile
lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_95,
  col = "blue",
  lwd = 1.5
)


legend(
  "topright",
  legend = c("5th", "50th", "95th"),
  cex = 0.7,
  col = c("red", "green", "blue"),
  lty = 1,
  lwd = 1.5
)

# Fit smoothed empirical centiles
loess_5 <- loess(
  empirical_5 ~ ga_wks,
  data = emp_agp,
  span = 0.3
)

loess_50 <- loess(
  empirical_50 ~ ga_wks,
  data = emp_agp,
  span = 0.3
)

loess_95 <- loess(
  empirical_95 ~ ga_wks,
  data = emp_agp,
  span = 0.3
)

# Predict onto the full GA grid
emp_agp <- emp_agp %>%
  mutate(
    smooth_5 = predict(loess_5, newdata = data.frame(ga_wks = ga_wks)),
    smooth_50 = predict(loess_50, newdata = data.frame(ga_wks = ga_wks)),
    smooth_95 = predict(loess_95, newdata = data.frame(ga_wks = ga_wks))
  )
lines(
  emp_agp$ga_wks,
  emp_agp$smooth_5,
  col = "red",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_agp$ga_wks,
  emp_agp$smooth_50,
  col = "green",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_agp$ga_wks,
  emp_agp$smooth_95,
  col = "blue",
  lty = 2,
  lwd = 1.5
)

#**************************************************************
# NLR BOX-COX CENTILE CURVES
#**************************************************************

ga_values <- seq(6, 40, by = 0.1)

centile_nlr_bc <- centiles.pred(
  nlr_boxcox,
  data = df_fpr_nlr,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()

names(centile_nlr_bc) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)

centile_nlr_bc <- centile_nlr_bc %>%
  mutate(
    centile_5  = inverse_boxcox(centile_5, lambda_nlr),
    centile_50 = inverse_boxcox(centile_50, lambda_nlr),
    centile_95 = inverse_boxcox(centile_95, lambda_nlr)
  )

plot(
  df_fpr_nlr$ga_wks,
  df_fpr_nlr$nlr_analysis,
  xlab = "Gestational Age (weeks)",
  ylab = "NLR",
  main = "NLR Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  xlim = c(6, 40),
  ylim = c(0, 10),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(centile_nlr_bc$ga_wks, centile_nlr_bc$centile_5,
      col = "red", lwd = 1.5)

lines(centile_nlr_bc$ga_wks, centile_nlr_bc$centile_50,
      col = "green", lwd = 1.5)

lines(centile_nlr_bc$ga_wks, centile_nlr_bc$centile_95,
      col = "blue", lwd = 1.5)

legend(
  "topright",
  legend = c("5th", "50th", "95th"),
  col = c("red", "green", "blue"),
  lty = 1,
  lwd = 1.5,
  cex = 0.7
)

# Fit smoothed empirical centiles
loess_5 <- loess(
  empirical_5 ~ ga_wks,
  data = emp_nlr,
  span = 0.3
)

loess_50 <- loess(
  empirical_50 ~ ga_wks,
  data = emp_nlr,
  span = 0.3
)

loess_95 <- loess(
  empirical_95 ~ ga_wks,
  data = emp_nlr,
  span = 0.3
)

# Predict onto the full GA grid
emp_agp <- emp_agp %>%
  mutate(
    smooth_5 = predict(loess_5, newdata = data.frame(ga_wks = ga_wks)),
    smooth_50 = predict(loess_50, newdata = data.frame(ga_wks = ga_wks)),
    smooth_95 = predict(loess_95, newdata = data.frame(ga_wks = ga_wks))
  )
lines(
  emp_agp$ga_wks,
  emp_agp$smooth_5,
  col = "red",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_agp$ga_wks,
  emp_agp$smooth_50,
  col = "green",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_agp$ga_wks,
  emp_agp$smooth_95,
  col = "blue",
  lty = 2,
  lwd = 1.5
)


#**************************************************************
# LMR BOX-COX CENTILE CURVES
#**************************************************************

centile_lmr_bc <- centiles.pred(
  lmr_boxcox,
  data = df_fpr_lmr,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()

names(centile_lmr_bc) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)

centile_lmr_bc <- centile_lmr_bc %>%
  mutate(
    centile_5  = inverse_boxcox(centile_5, lambda_lmr),
    centile_50 = inverse_boxcox(centile_50, lambda_lmr),
    centile_95 = inverse_boxcox(centile_95, lambda_lmr)
  )

plot(
  df_fpr_lmr$ga_wks,
  df_fpr_lmr$lmr_analysis,
  xlab = "Gestational Age (weeks)",
  ylab = "LMR",
  main = "LMR Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  xlim = c(6, 40),
  ylim = c(0, 10),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(centile_lmr_bc$ga_wks, centile_lmr_bc$centile_5,
      col = "red", lwd = 1.5)

lines(centile_lmr_bc$ga_wks, centile_lmr_bc$centile_50,
      col = "green", lwd = 1.5)

lines(centile_lmr_bc$ga_wks, centile_lmr_bc$centile_95,
      col = "blue", lwd = 1.5)

legend(
  "topright",
  legend = c("5th", "50th", "95th"),
  col = c("red", "green", "blue"),
  lty = 1,
  lwd = 1.5,
  cex = 0.7
)
loess_5 <- loess(
  empirical_5 ~ ga_wks,
  data = emp_lmr,
  span = 0.3
)

loess_50 <- loess(
  empirical_50 ~ ga_wks,
  data = emp_lmr,
  span = 0.3
)

loess_95 <- loess(
  empirical_95 ~ ga_wks,
  data = emp_lmr,
  span = 0.3
)

# Predict onto the full GA grid
emp_lmr <- emp_lmr %>%
  mutate(
    smooth_5 = predict(loess_5, newdata = data.frame(ga_wks = ga_wks)),
    smooth_50 = predict(loess_50, newdata = data.frame(ga_wks = ga_wks)),
    smooth_95 = predict(loess_95, newdata = data.frame(ga_wks = ga_wks))
  )
lines(
  emp_lmr$ga_wks,
  emp_lmr$smooth_5,
  col = "red",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_lmr$ga_wks,
  emp_lmr$smooth_50,
  col = "green",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_lmr$ga_wks,
  emp_lmr$smooth_95,
  col = "blue",
  lty = 2,
  lwd = 1.5
)



#**************************************************************
# PLR BOX-COX CENTILE CURVES
#**************************************************************

centile_plr_bc <- centiles.pred(
  plr_boxcox,
  data = df_fpr_plr,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()

names(centile_plr_bc) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)

centile_plr_bc <- centile_plr_bc %>%
  mutate(
    centile_5  = inverse_boxcox(centile_5, lambda_plr),
    centile_50 = inverse_boxcox(centile_50, lambda_plr),
    centile_95 = inverse_boxcox(centile_95, lambda_plr)
  )

plot(
  df_fpr_plr$ga_wks,
  df_fpr_plr$plr_analysis,
  xlab = "Gestational Age (weeks)",
  ylab = "PLR",
  main = "PLR Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  xlim = c(6, 40),
  ylim = c(0, 400),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(centile_plr_bc$ga_wks, centile_plr_bc$centile_5,
      col = "red", lwd = 1.5)

lines(centile_plr_bc$ga_wks, centile_plr_bc$centile_50,
      col = "green", lwd = 1.5)

lines(centile_plr_bc$ga_wks, centile_plr_bc$centile_95,
      col = "blue", lwd = 1.5)

legend(
  "topright",
  legend = c("5th", "50th", "95th"),
  col = c("red", "green", "blue"),
  lty = 1,
  lwd = 1.5,
  cex = 0.7
)

loess_5 <- loess(
  empirical_5 ~ ga_wks,
  data = emp_plr,
  span = 0.3
)

loess_50 <- loess(
  empirical_50 ~ ga_wks,
  data = emp_plr,
  span = 0.3
)

loess_95 <- loess(
  empirical_95 ~ ga_wks,
  data = emp_plr,
  span = 0.3
)

# Predict onto the full GA grid
emp_plr <- emp_plr %>%
  mutate(
    smooth_5 = predict(loess_5, newdata = data.frame(ga_wks = ga_wks)),
    smooth_50 = predict(loess_50, newdata = data.frame(ga_wks = ga_wks)),
    smooth_95 = predict(loess_95, newdata = data.frame(ga_wks = ga_wks))
  )
lines(
  emp_plr$ga_wks,
  emp_plr$smooth_5,
  col = "red",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_plr$ga_wks,
  emp_plr$smooth_50,
  col = "green",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_plr$ga_wks,
  emp_plr$smooth_95,
  col = "blue",
  lty = 2,
  lwd = 1.5
)
#**************************************************************
# SII BOX-COX CENTILE CURVES
#**************************************************************

centile_sii_bc <- centiles.pred(
  sii_boxcox,
  data = df_fpr_sii,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()

names(centile_sii_bc) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)

centile_sii_bc <- centile_sii_bc %>%
  mutate(
    centile_5  = inverse_boxcox(centile_5, lambda_sii),
    centile_50 = inverse_boxcox(centile_50, lambda_sii),
    centile_95 = inverse_boxcox(centile_95, lambda_sii)
  )

plot(
  df_fpr_sii$ga_wks,
  df_fpr_sii$sii_analysis,
  xlab = "Gestational Age (weeks)",
  ylab = "SII",
  main = "SII Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  xlim = c(6, 40),
  ylim = c(0, 2500),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(centile_sii_bc$ga_wks, centile_sii_bc$centile_5,
      col = "red", lwd = 1.5)

lines(centile_sii_bc$ga_wks, centile_sii_bc$centile_50,
      col = "green", lwd = 1.5)

lines(centile_sii_bc$ga_wks, centile_sii_bc$centile_95,
      col = "blue", lwd = 1.5)

legend(
  "topright",
  legend = c("5th", "50th", "95th"),
  col = c("red", "green", "blue"),
  lty = 1,
  lwd = 1.5,
  cex = 0.7
)

loess_5 <- loess(
  empirical_5 ~ ga_wks,
  data = emp_sii,
  span = 0.3
)

loess_50 <- loess(
  empirical_50 ~ ga_wks,
  data = emp_sii,
  span = 0.3
)

loess_95 <- loess(
  empirical_95 ~ ga_wks,
  data = emp_sii,
  span = 0.3
)

# Predict onto the full GA grid
emp_sii <- emp_sii %>%
  mutate(
    smooth_5 = predict(loess_5, newdata = data.frame(ga_wks = ga_wks)),
    smooth_50 = predict(loess_50, newdata = data.frame(ga_wks = ga_wks)),
    smooth_95 = predict(loess_95, newdata = data.frame(ga_wks = ga_wks))
  )
lines(
  emp_sii$ga_wks,
  emp_sii$smooth_5,
  col = "red",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_sii$ga_wks,
  emp_sii$smooth_50,
  col = "green",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_sii$ga_wks,
  emp_sii$smooth_95,
  col = "blue",
  lty = 2,
  lwd = 1.5
)

#------------- Add in Empirical Curves ------------

empirical_centiles <- function(data, outcome) {
  
  data %>%
    filter(
      !is.na(ga_wks),
      !is.na(.data[[outcome]])
    ) %>%
    group_by(ga_wks) %>%
    summarise(
      empirical_5  = quantile(.data[[outcome]], 0.05, na.rm = TRUE),
      empirical_50 = quantile(.data[[outcome]], 0.50, na.rm = TRUE),
      empirical_95 = quantile(.data[[outcome]], 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}
emp_agp <- empirical_centiles(df_fpr_agp, "agp")
emp_nlr <- empirical_centiles(df_fpr_nlr, "nlr_analysis")
emp_lmr <- empirical_centiles(df_fpr_lmr, "lmr_analysis")
emp_plr <- empirical_centiles(df_fpr_plr, "plr_analysis")
emp_sii <- empirical_centiles(df_fpr_sii, "sii_analysis")

#**************************************************************
# NORMAL APPROXIMATION - ALL BIOMARKERS
#**************************************************************

library(dplyr)

get_normal_centiles <- function(model, data, lambda, marker) {
  
  # Mean and SD from fitted FPR model
  mu <- fitted(model, what = "mu")
  sigma <- fitted(model, what = "sigma")
  
  # Calculate centiles and back-transform
  data.frame(
    Biomarker = marker,
    ga_wks = data$ga_wks,
    P5 = inverse_boxcox(mu - 1.645 * sigma, lambda),
    P50 = inverse_boxcox(mu, lambda),
    P95 = inverse_boxcox(mu + 1.645 * sigma, lambda)
  )
}

#**************************************************************
# COMBINE ALL FIVE BIOMARKERS
#**************************************************************

normal_centiles <- bind_rows(
  
  get_normal_centiles(
    agp_boxcox, df_fpr_agp_35, lambda_agp, "AGP"
  ),
  
  get_normal_centiles(
    nlr_boxcox, df_fpr_nlr_35, lambda_nlr, "NLR"
  ),
  
  get_normal_centiles(
    lmr_boxcox, df_fpr_lmr_35, lambda_lmr, "LMR"
  ),
  
  get_normal_centiles(
    plr_boxcox, df_fpr_plr_35, lambda_plr, "PLR"
  ),
  
  get_normal_centiles(
    sii_boxcox, df_fpr_sii_35, lambda_sii, "SII"
  )
)

# View results
View(normal_centiles)

