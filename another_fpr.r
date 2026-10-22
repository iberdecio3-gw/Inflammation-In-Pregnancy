#*********************************
# Analysis 1:
# FPR Reference Curves for CRP by Gestational Age
# Distribution assumptions:
# Box-Cox, Gamma, Inverse Gaussian, Weibull
#*********************************

########## HEALTHY COHORT ################
library(dplyr)
library(gamlss)
library(gamlss.dist)
library(ggplot2)

con <- gamlss.control(trace = FALSE)

df_fpr <- crp.model %>%
  dplyr::select(site, ga_wks, crp, imputed_crp, delta) %>%
  filter(
    !is.na(ga_wks),
    !is.na(imputed_crp),
    imputed_crp > 0
  )

# Gamma
fpr_gamma <- gamlss(
  imputed_crp ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = GA(),
  data = df_fpr,
  control = con
)

# Inverse Gaussian
fpr_invgauss <- gamlss(
  imputed_crp ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = IG(),
  data = df_fpr,
  control = con
)

# Weibull
fpr_weibull <- gamlss(
  imputed_crp ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = WEI(),
  data = df_fpr,
  control = con
)

# Box-Cox transformation
# Estimate Box-Cox lambda using the cleaned positive CRP data

bc_fit <- MASS::boxcox(
  lm(imputed_crp ~ ga_wks, data = df_fpr),
  plotit = FALSE
)

lambda <- bc_fit$x[which.max(bc_fit$y)]

lambda


if (abs(lambda) < 1e-6) {
  
  df_crp_imp_bc <- df_fpr %>%
    mutate(
      crp_bc = log(imputed_crp)
    )
  
} else {
  
  df_crp_imp_bc <- df_fpr %>%
    mutate(
      crp_bc = (imputed_crp^lambda - 1) / lambda
    )
}

df_fpr_bc <- df_crp_imp_bc %>%
  dplyr::select(
    site,
    ga_wks,
    imputed_crp,
    crp_bc
  )

fpr_boxcox <- gamlss(
  crp_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_bc,
  control = con
)


#Compare Model Fit
GAIC(
  fpr_boxcox,
  fpr_gamma,
  fpr_invgauss,
  fpr_weibull
)


### pLOTS: Centile curves using FPR model (5th, 50th, 95th centiles)

centiles_fun <- function(model, data, title) {
  
  centiles(
    model,
    xvar = data$ga_wks,
    cent = c(5, 50, 95),
    legend = FALSE,
    xlab = "Gestational Age (weeks)",
    ylab = "CRP (mg/L)",
    main = title,
    bg = "transparent",
    xlim = range(6, 40),
    ylim = c(0, 20),
    pch = 15,
    cex = 0.2,
    col = "grey86",
    col.cent = c("red", "green", "blue")
  )
  
  legend(
    "topright",
    inset = c(0, 0),
    legend = c(5, 50, 95),
    cex = 0.6,
    col = c("red", "green", "blue"),
    lty = c(1, 1, 1)
  )
}


ga_values <- 6:40

#GAMMA CURVES + VALUES
centiles_gamma <- centiles_fun(
  fpr_gamma,
  df_fpr,
  "Gamma Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)"
)

centile_gamma <- centiles.pred(
  fpr_gamma,
  data = df_fpr,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()

names(centile_gamma) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)



centiles_invgauss <- centiles_fun(
  fpr_invgauss,
  df_fpr,
  "Inverse Gaussian Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)"
)

centile_invgauss <- centiles.pred(
  fpr_invgauss,
  data = df_fpr,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()

names(centile_invgauss) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)

centiles_weibull <- centiles_fun(
  fpr_weibull,
  df_fpr,
  "Weibull Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)"
)
centile_weibull <- centiles.pred(
  fpr_weibull,
  data = df_fpr,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()

names(centile_weibull) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)

#**************************************************************
# Box-Cox transformation
#**************************************************************

# Estimate Box-Cox lambda
bc_fit <- MASS::boxcox(
  lm(imputed_crp ~ ga_wks, data = df_fpr),
  plotit = FALSE
)

lambda <- bc_fit$x[which.max(bc_fit$y)]
lambda


# Transform CRP
if (abs(lambda) < 1e-6) {
  
  df_crp_imp_bc <- df_fpr %>%
    mutate(
      crp_bc = log(imputed_crp)
    )
  
} else {
  
  df_crp_imp_bc <- df_fpr %>%
    mutate(
      crp_bc = (imputed_crp^lambda - 1) / lambda
    )
  
}


df_fpr_bc <- df_crp_imp_bc %>%
  dplyr::select(
    site,
    ga_wks,
    imputed_crp,
    crp_bc
  )


# Fit Box-Cox model
fpr_boxcox <- gamlss(
  crp_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_bc,
  control = con
)


#**************************************************************
# Inverse Box-Cox function
#**************************************************************

inverse_boxcox <- function(y, lambda){
  
  if(abs(lambda) < 1e-6){
    
    exp(y)
    
  } else{
    
    (lambda * y + 1)^(1/lambda)
    
  }
  
}


#**************************************************************
# Plot Box-Cox centiles
#**************************************************************

centiles_bc_pred <- centiles.pred(
  fpr_boxcox,
  data = df_fpr_bc,
  xname = "ga_wks",
  xvalues = seq(6,40,0.1),
  type = "centiles",
  cent = c(5,50,95)
) %>%
  as.data.frame()

names(centiles_bc_pred) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)


# Back-transform
centiles_bc_pred <- centiles_bc_pred %>%
  mutate(
    centile_5  = inverse_boxcox(centile_5,lambda),
    centile_50 = inverse_boxcox(centile_50,lambda),
    centile_95 = inverse_boxcox(centile_95,lambda)
  )


# Plot
plot(
  df_fpr_bc$ga_wks,
  df_fpr_bc$imputed_crp,
  xlab = "Gestational Age (weeks)",
  ylab = "CRP (mg/L)",
  main = "Box-Cox Transformation Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  bg = "transparent",
  xlim = c(6,40),
  ylim = c(0,20),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(
  centiles_bc_pred$ga_wks,
  centiles_bc_pred$centile_5,
  col = "red"
)

lines(
  centiles_bc_pred$ga_wks,
  centiles_bc_pred$centile_50,
  col = "green"
)

lines(
  centiles_bc_pred$ga_wks,
  centiles_bc_pred$centile_95,
  col = "blue"
)

legend(
  "topright",
  inset = c(0,0),
  legend = c(5,50,95),
  cex = 0.6,
  col = c("red","green","blue"),
  lty = 1
)


#**************************************************************
# Numeric centiles (weeks 6-40)
#**************************************************************

centile_boxcox <- centiles.pred(
  fpr_boxcox,
  data = df_fpr_bc,
  xname = "ga_wks",
  xvalues = 6:40,
  type = "centiles",
  cent = c(5,50,95)
) %>%
  as.data.frame()

names(centile_boxcox) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)


# Back-transform to CRP (mg/L)
centile_boxcox <- centile_boxcox %>%
  mutate(
    centile_5  = inverse_boxcox(centile_5,lambda),
    centile_50 = inverse_boxcox(centile_50,lambda),
    centile_95 = inverse_boxcox(centile_95,lambda)
  )

centile_boxcox

########### Centile Values by Distribution
selected_weeks <- c(5, 10, 15, 20, 25, 30, 35, 40)

centile_values <- bind_rows(
  centile_boxcox %>%
    mutate(Distribution = "Box-Cox"),
  
  centile_gamma %>%
    mutate(Distribution = "Gamma"),
  
  centile_invgauss %>%
    mutate(Distribution = "Inverse Gaussian"),
  
  centile_weibull %>%
    mutate(Distribution = "Weibull")
) %>%
  filter(ga_wks %in% selected_weeks) %>%
  dplyr::select(
    Distribution,
    ga_wks,
    centile_5,
    centile_50,
    centile_95
  ) %>%
  arrange(Distribution, ga_wks) %>%
  mutate(
    across(
      c(centile_5, centile_50, centile_95),
      round,
      digits = 2
    )
  )
gt(centile_values)

#box cox only
centile_values <- centile_boxcox %>%
    mutate(Distribution = "Box-Cox") %>%
  filter(ga_wks %in% selected_weeks) %>%
  dplyr::select(
    Distribution,
    ga_wks,
    centile_5,
    centile_50,
    centile_95
  ) %>%
  arrange(Distribution, ga_wks) %>%
  mutate(
    across(
      c(centile_5, centile_50, centile_95),
      round,
      digits = 2
    )
  )
gt(centile_values)

centile_table <- centile_boxcox %>% mutate(Distribution = "Box-Cox") %>%
  filter(ga_wks %in% selected_weeks) %>%
  tidyr::pivot_longer(
    cols = c(centile_5, centile_50, centile_95),
    names_to = "Centile",
    values_to = "Value"
  ) %>%
  mutate(
    Centile = recode(
      Centile,
      centile_5 = "5th",
      centile_50 = "50th",
      centile_95 = "95th"
    ),
    ga_wks = paste0("Week ", ga_wks)
  ) %>%
  tidyr::pivot_wider(
    names_from = ga_wks,
    values_from = Value
  ) %>%
  mutate(
    across(starts_with("Week"), \(x) round(x, 2))
  )

library(gt)

centile_table %>%
  gt(
    groupname_col = "Distribution",
    rowname_col = "Centile"
  ) %>%
  cols_label(
    Centile = "Centile",
    `Week 10` = "10",
    `Week 15` = "15",
    `Week 20` = "20",
    `Week 25` = "25",
    `Week 30` = "30",
    `Week 35` = "35",
    `Week 40` = "40"
  ) %>%
  tab_stubhead(label = "Centile") %>%
  tab_spanner(
    label = "Gestational age (weeks)",
    columns = starts_with("Week")
  )



#MACE
coverage <- data.frame(
  Distribution = c("Box-Cox", "Gamma", "Inverse Gaussian", "Weibull"),
  P5  = c(4.4722, 6.678082, 6.144239, 7.534247),
  P50 = c(48.6805, 45.3362, 35.19339, 47.33078),
  P95 = c(99.58703, 100, 100.0000, 99.76833)
)

coverage$MACE <- with(
  coverage,
  (abs(P5 - 5) + abs(P50 - 50) + abs(P95 - 95)) / 3
)

coverage


library(dplyr)

comparison_table <- centile_boxcox %>%
  rename(
    Healthy_5 = centile_5,
    Healthy_50 = centile_50,
    Healthy_95 = centile_95
  ) %>%
  inner_join(
    centile_boxcox_nh %>%
      rename(
        NonHealthy_5 = centile_5,
        NonHealthy_50 = centile_50,
        NonHealthy_95 = centile_95
      ),
    by = "ga_wks"
  ) %>%
  mutate(
    Difference_5 = NonHealthy_5 - Healthy_5,
    Difference_50 = NonHealthy_50 - Healthy_50,
    Difference_95 = NonHealthy_95 - Healthy_95
  ) %>%
  mutate(
    across(
      -ga_wks,
      \(x) round(x, 2)
    )
  )

comparison_table
gt(comparison_table)




plot1 <- df_fpr_agp %>% plot(agp ~ ga_wks)

hist(df_fpr_agp$agp, breaks = 30)

boxplot(df_fpr_agp$agp)

summary(df_fpr_agp$agp)

qqnorm(df$agp)
qqline(df$agp)

