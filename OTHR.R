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
# Restrict analysis to <= 35 weeks
#**************************************************************
{
df_fpr_agp_35 <- df_fpr_agp %>%
  filter(ga_wks <= 35)

df_fpr_nlr_35 <- df_fpr_nlr %>%
  filter(ga_wks <= 35)

df_fpr_lmr_35 <- df_fpr_lmr %>%
  filter(ga_wks <= 35)

df_fpr_plr_35 <- df_fpr_plr %>%
  filter(ga_wks <= 35)

df_fpr_sii_35 <- df_fpr_sii %>%
  filter(ga_wks <= 35)


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


emp_agp <- empirical_centiles(
  df_fpr_agp_35,
  "agp"
)

emp_nlr <- empirical_centiles(
  df_fpr_nlr_35,
  "nlr_analysis"
)

emp_lmr <- empirical_centiles(
  df_fpr_lmr_35,
  "lmr_analysis"
)

emp_plr <- empirical_centiles(
  df_fpr_plr_35,
  "plr_analysis"
)

emp_sii <- empirical_centiles(
  df_fpr_sii_35,
  "sii_analysis"
)

#**************************************************************
# AGP BOX-COX TRANSFORMATION
# Analysis restricted to <= 35 weeks
#**************************************************************

bc_agp <- MASS::boxcox(
  lm(
    agp ~ ga_wks,
    data = df_fpr_agp_35
  ),
  plotit = FALSE
)

lambda_agp <- bc_agp$x[
  which.max(bc_agp$y)
]


if (abs(lambda_agp) < 1e-6) {
  
  df_fpr_agp_35 <- df_fpr_agp_35 %>%
    mutate(
      agp_bc = log(agp)
    )
  
} else {
  
  df_fpr_agp_35 <- df_fpr_agp_35 %>%
    mutate(
      agp_bc = (agp^lambda_agp - 1) / lambda_agp
    )
}


agp_boxcox <- gamlss(
  agp_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_agp_35,
  control = con
)

plot(
  df_fpr_agp_35$ga_wks,
  df_fpr_agp_35$agp,
  xlab = "Gestational Age (weeks)",
  ylab = "AGP",
  main = "AGP Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  bg = "transparent",
  xlim = c(6, 35),
  ylim = c(0, 2),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

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

############ ----- Plot Box Cox --------- #####

#**************************************************************
# AGP BOX-COX CENTILE CURVES
# Analysis and plot restricted to <= 35 weeks
#**************************************************************

ga_values <- seq(6, 35, by = 0.1)


# Predict centiles on Box-Cox transformed scale
centile_agp_bc <- centiles.pred(
  agp_boxcox,
  data = df_fpr_agp_35,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()


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


summary(centile_agp_bc)


#**************************************************************
# Plot
#**************************************************************

plot(
  df_fpr_agp_35$ga_wks,
  df_fpr_agp_35$agp,
  xlab = "Gestational Age (weeks)",
  ylab = "AGP",
  main = "AGP Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  bg = "transparent",
  xlim = c(6, 35),
  ylim = c(0, 2),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)


#**************************************************************
# Fitted centile curves
#**************************************************************

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


#**************************************************************
# Smoothed empirical centile curves
# emp_agp should already be based on df_fpr_agp_35
#**************************************************************

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


# Predict smoothed empirical centiles
emp_agp <- emp_agp %>%
  mutate(
    smooth_5 = predict(
      loess_5,
      newdata = data.frame(ga_wks = ga_wks)
    ),
    
    smooth_50 = predict(
      loess_50,
      newdata = data.frame(ga_wks = ga_wks)
    ),
    
    smooth_95 = predict(
      loess_95,
      newdata = data.frame(ga_wks = ga_wks)
    )
  )


#**************************************************************
# Add empirical centiles as dashed lines
#**************************************************************

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
}

emp_agp <- df_fpr_agp_35 %>%
  filter(
    !is.na(ga_wks),
    !is.na(agp)
  ) %>%
  mutate(
    ga_bin = round(ga_wks * 2) / 2
  ) %>%
  group_by(ga_bin) %>%
  summarise(
    ga_wks = mean(ga_wks),
    empirical_5  = quantile(agp, 0.05, na.rm = TRUE),
    empirical_50 = quantile(agp, 0.50, na.rm = TRUE),
    empirical_95 = quantile(agp, 0.95, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
loess_5 <- loess(
  empirical_5 ~ ga_wks,
  data = emp_agp,
  weights = n,
  span = 0.5
)

loess_50 <- loess(
  empirical_50 ~ ga_wks,
  data = emp_agp,
  weights = n,
  span = 0.5
)

loess_95 <- loess(
  empirical_95 ~ ga_wks,
  data = emp_agp,
  weights = n,
  span = 0.5
)

emp_agp <- emp_agp %>%
  mutate(
    smooth_5 = predict(
      loess_5,
      newdata = data.frame(ga_wks = ga_wks)
    ),
    smooth_50 = predict(
      loess_50,
      newdata = data.frame(ga_wks = ga_wks)
    ),
    smooth_95 = predict(
      loess_95,
      newdata = data.frame(ga_wks = ga_wks)
    )
  )
df_fpr_agp_35 %>%
  count(ga_wks) %>%
  arrange(n)

#**************************************************************
# Predict smoothed empirical centiles
#**************************************************************

emp_agp <- emp_agp %>%
  mutate(
    smooth_5 = predict(
      loess_5,
      newdata = data.frame(ga_wks = ga_wks)
    ),
    
    smooth_50 = predict(
      loess_50,
      newdata = data.frame(ga_wks = ga_wks)
    ),
    
    smooth_95 = predict(
      loess_95,
      newdata = data.frame(ga_wks = ga_wks)
    )
  )


#**************************************************************
# Plot AGP
#**************************************************************

plot(
  df_fpr_agp_35$ga_wks,
  df_fpr_agp_35$agp,
  xlab = "Gestational Age (weeks)",
  ylab = "AGP",
  main = "AGP Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  bg = "transparent",
  xlim = c(6, 35),
  ylim = c(0, 2),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)


#**************************************************************
# Fitted Box-Cox centile curves - solid
#**************************************************************

lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_5,
  col = "red",
  lwd = 1.5
)

lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_50,
  col = "green",
  lwd = 1.5
)

lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_95,
  col = "blue",
  lwd = 1.5
)


#**************************************************************
# Smoothed empirical centiles - dashed
#**************************************************************

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
# Legend
#**************************************************************

legend(
  "topright",
  legend = c(
    "5th fitted",
    "50th fitted",
    "95th fitted",
    "5th empirical",
    "50th empirical",
    "95th empirical"
  ),
  col = c(
    "red",
    "green",
    "blue",
    "red",
    "green",
    "blue"
  ),
  lty = c(
    1, 1, 1,
    2, 2, 2
  ),
  lwd = 1.5,
  cex = 0.65
)

library(dplyr)
library(gamlss)
library(gamlss.dist)
library(MASS)


###############################################################
# BOX-COX ANALYSIS
# ANALYSIS RESTRICTED TO <= 35 WEEKS
###############################################################


#**************************************************************
# Helper: inverse Box-Cox transformation
#**************************************************************

inverse_boxcox <- function(y, lambda) {
  
  if (abs(lambda) < 1e-6) {
    
    exp(y)
    
  } else {
    
    (lambda * y + 1)^(1 / lambda)
    
  }
}


#**************************************************************
# Helper: empirical centiles using 0.5-week bins
#**************************************************************

empirical_halfweek <- function(data, outcome) {
  
  data %>%
    filter(
      !is.na(ga_wks),
      !is.na(.data[[outcome]])
    ) %>%
    
    mutate(
      ga_bin = round(ga_wks * 2) / 2
    ) %>%
    
    group_by(ga_bin) %>%
    
    summarise(
      empirical_5 = quantile(
        .data[[outcome]],
        0.05,
        na.rm = TRUE
      ),
      
      empirical_50 = quantile(
        .data[[outcome]],
        0.50,
        na.rm = TRUE
      ),
      
      empirical_95 = quantile(
        .data[[outcome]],
        0.95,
        na.rm = TRUE
      ),
      
      n = n(),
      .groups = "drop"
    ) %>%
    
    rename(
      ga_wks = ga_bin
    )
}


#**************************************************************
# Helper: smooth empirical centiles
#**************************************************************

smooth_empirical <- function(emp_data, span = 0.5) {
  
  loess_5 <- loess(
    empirical_5 ~ ga_wks,
    data = emp_data,
    weights = n,
    span = span
  )
  
  loess_50 <- loess(
    empirical_50 ~ ga_wks,
    data = emp_data,
    weights = n,
    span = span
  )
  
  loess_95 <- loess(
    empirical_95 ~ ga_wks,
    data = emp_data,
    weights = n,
    span = span
  )
  
  
  emp_data %>%
    mutate(
      smooth_5 = predict(
        loess_5,
        newdata = data.frame(
          ga_wks = ga_wks
        )
      ),
      
      smooth_50 = predict(
        loess_50,
        newdata = data.frame(
          ga_wks = ga_wks
        )
      ),
      
      smooth_95 = predict(
        loess_95,
        newdata = data.frame(
          ga_wks = ga_wks
        )
      )
    )
}


#**************************************************************
# Restrict all analysis datasets to <= 35 weeks
#**************************************************************

df_fpr_agp_35 <- df_fpr_agp %>%
  filter(ga_wks <= 35)

df_fpr_nlr_35 <- df_fpr_nlr %>%
  filter(ga_wks <= 35)

df_fpr_lmr_35 <- df_fpr_lmr %>%
  filter(ga_wks <= 35)

df_fpr_plr_35 <- df_fpr_plr %>%
  filter(ga_wks <= 35)

df_fpr_sii_35 <- df_fpr_sii %>%
  filter(ga_wks <= 35)


# Prediction grid
ga_values <- seq(
  6,
  35,
  by = 0.1
)



###############################################################
# AGP BOX-COX
###############################################################


#**************************************************************
# Estimate Box-Cox lambda
#**************************************************************

bc_agp <- MASS::boxcox(
  lm(
    agp ~ ga_wks,
    data = df_fpr_agp_35
  ),
  plotit = FALSE
)

lambda_agp <- bc_agp$x[
  which.max(bc_agp$y)
]

lambda_agp


#**************************************************************
# Transform AGP
#**************************************************************

if (abs(lambda_agp) < 1e-6) {
  
  df_fpr_agp_35 <- df_fpr_agp_35 %>%
    mutate(
      agp_bc = log(agp)
    )
  
} else {
  
  df_fpr_agp_35 <- df_fpr_agp_35 %>%
    mutate(
      agp_bc = (
        agp^lambda_agp - 1
      ) / lambda_agp
    )
}


#**************************************************************
# Fit AGP Box-Cox GAMLSS
#**************************************************************

agp_boxcox <- gamlss(
  agp_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_agp_35,
  control = con
)


#**************************************************************
# AGP fitted centiles
#**************************************************************

centile_agp_bc <- centiles.pred(
  agp_boxcox,
  data = df_fpr_agp_35,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()


names(centile_agp_bc) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)


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


#**************************************************************
# AGP empirical centiles: 0.5-week bins
#**************************************************************

emp_agp <- empirical_halfweek(
  df_fpr_agp_35,
  "agp"
)

emp_agp <- smooth_empirical(
  emp_agp,
  span = 0.5
)



###############################################################
# NLR BOX-COX
###############################################################


#**************************************************************
# Estimate Box-Cox lambda
#**************************************************************

bc_nlr <- MASS::boxcox(
  lm(
    nlr_analysis ~ ga_wks,
    data = df_fpr_nlr_35
  ),
  plotit = FALSE
)

lambda_nlr <- bc_nlr$x[
  which.max(bc_nlr$y)
]

lambda_nlr


#**************************************************************
# Transform NLR
#**************************************************************

if (abs(lambda_nlr) < 1e-6) {
  
  df_fpr_nlr_35 <- df_fpr_nlr_35 %>%
    mutate(
      nlr_bc = log(nlr_analysis)
    )
  
} else {
  
  df_fpr_nlr_35 <- df_fpr_nlr_35 %>%
    mutate(
      nlr_bc = (
        nlr_analysis^lambda_nlr - 1
      ) / lambda_nlr
    )
}


#**************************************************************
# Fit NLR Box-Cox GAMLSS
#**************************************************************

nlr_boxcox <- gamlss(
  nlr_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_nlr_35,
  control = con
)


#**************************************************************
# NLR fitted centiles
#**************************************************************

centile_nlr_bc <- centiles.pred(
  nlr_boxcox,
  data = df_fpr_nlr_35,
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
    centile_5 = inverse_boxcox(
      centile_5,
      lambda_nlr
    ),
    
    centile_50 = inverse_boxcox(
      centile_50,
      lambda_nlr
    ),
    
    centile_95 = inverse_boxcox(
      centile_95,
      lambda_nlr
    )
  )


#**************************************************************
# NLR empirical centiles: 0.5-week bins
#**************************************************************

emp_nlr <- empirical_halfweek(
  df_fpr_nlr_35,
  "nlr_analysis"
)

emp_nlr <- smooth_empirical(
  emp_nlr,
  span = 0.5
)



###############################################################
# LMR BOX-COX
###############################################################


#**************************************************************
# Estimate Box-Cox lambda
#**************************************************************

bc_lmr <- MASS::boxcox(
  lm(
    lmr_analysis ~ ga_wks,
    data = df_fpr_lmr_35
  ),
  plotit = FALSE
)

lambda_lmr <- bc_lmr$x[
  which.max(bc_lmr$y)
]

lambda_lmr


#**************************************************************
# Transform LMR
#**************************************************************

if (abs(lambda_lmr) < 1e-6) {
  
  df_fpr_lmr_35 <- df_fpr_lmr_35 %>%
    mutate(
      lmr_bc = log(lmr_analysis)
    )
  
} else {
  
  df_fpr_lmr_35 <- df_fpr_lmr_35 %>%
    mutate(
      lmr_bc = (
        lmr_analysis^lambda_lmr - 1
      ) / lambda_lmr
    )
}


#**************************************************************
# Fit LMR Box-Cox GAMLSS
#**************************************************************

lmr_boxcox <- gamlss(
  lmr_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_lmr_35,
  control = con
)


#**************************************************************
# LMR fitted centiles
#**************************************************************

centile_lmr_bc <- centiles.pred(
  lmr_boxcox,
  data = df_fpr_lmr_35,
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
    centile_5 = inverse_boxcox(
      centile_5,
      lambda_lmr
    ),
    
    centile_50 = inverse_boxcox(
      centile_50,
      lambda_lmr
    ),
    
    centile_95 = inverse_boxcox(
      centile_95,
      lambda_lmr
    )
  )


#**************************************************************
# LMR empirical centiles: 0.5-week bins
#**************************************************************

emp_lmr <- empirical_halfweek(
  df_fpr_lmr_35,
  "lmr_analysis"
)

emp_lmr <- smooth_empirical(
  emp_lmr,
  span = 0.5
)



###############################################################
# PLR BOX-COX
###############################################################


#**************************************************************
# Estimate Box-Cox lambda
#**************************************************************

bc_plr <- MASS::boxcox(
  lm(
    plr_analysis ~ ga_wks,
    data = df_fpr_plr_35
  ),
  plotit = FALSE
)

lambda_plr <- bc_plr$x[
  which.max(bc_plr$y)
]

lambda_plr


#**************************************************************
# Transform PLR
#**************************************************************

if (abs(lambda_plr) < 1e-6) {
  
  df_fpr_plr_35 <- df_fpr_plr_35 %>%
    mutate(
      plr_bc = log(plr_analysis)
    )
  
} else {
  
  df_fpr_plr_35 <- df_fpr_plr_35 %>%
    mutate(
      plr_bc = (
        plr_analysis^lambda_plr - 1
      ) / lambda_plr
    )
}


#**************************************************************
# Fit PLR Box-Cox GAMLSS
#**************************************************************

plr_boxcox <- gamlss(
  plr_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_plr_35,
  control = con
)


#**************************************************************
# PLR fitted centiles
#**************************************************************

centile_plr_bc <- centiles.pred(
  plr_boxcox,
  data = df_fpr_plr_35,
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
    centile_5 = inverse_boxcox(
      centile_5,
      lambda_plr
    ),
    
    centile_50 = inverse_boxcox(
      centile_50,
      lambda_plr
    ),
    
    centile_95 = inverse_boxcox(
      centile_95,
      lambda_plr
    )
  )


#**************************************************************
# PLR empirical centiles: 0.5-week bins
#**************************************************************

emp_plr <- empirical_halfweek(
  df_fpr_plr_35,
  "plr_analysis"
)

emp_plr <- smooth_empirical(
  emp_plr,
  span = 0.5
)



###############################################################
# SII BOX-COX
###############################################################


#**************************************************************
# Estimate Box-Cox lambda
#**************************************************************

bc_sii <- MASS::boxcox(
  lm(
    sii_analysis ~ ga_wks,
    data = df_fpr_sii_35
  ),
  plotit = FALSE
)

lambda_sii <- bc_sii$x[
  which.max(bc_sii$y)
]

lambda_sii


#**************************************************************
# Transform SII
#**************************************************************

if (abs(lambda_sii) < 1e-6) {
  
  df_fpr_sii_35 <- df_fpr_sii_35 %>%
    mutate(
      sii_bc = log(sii_analysis)
    )
  
} else {
  
  df_fpr_sii_35 <- df_fpr_sii_35 %>%
    mutate(
      sii_bc = (
        sii_analysis^lambda_sii - 1
      ) / lambda_sii
    )
}


#**************************************************************
# Fit SII Box-Cox GAMLSS
#**************************************************************

sii_boxcox <- gamlss(
  sii_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = df_fpr_sii_35,
  control = con
)


#**************************************************************
# SII fitted centiles
#**************************************************************

centile_sii_bc <- centiles.pred(
  sii_boxcox,
  data = df_fpr_sii_35,
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
    centile_5 = inverse_boxcox(
      centile_5,
      lambda_sii
    ),
    
    centile_50 = inverse_boxcox(
      centile_50,
      lambda_sii
    ),
    
    centile_95 = inverse_boxcox(
      centile_95,
      lambda_sii
    )
  )


#**************************************************************
# SII empirical centiles: 0.5-week bins
#**************************************************************

emp_sii <- empirical_halfweek(
  df_fpr_sii_35,
  "sii_analysis"
)

emp_sii <- smooth_empirical(
  emp_sii,
  span = 0.5
)



###############################################################
# BOX-COX LAMBDA SUMMARY
###############################################################

tibble(
  Marker = c(
    "AGP",
    "NLR",
    "LMR",
    "PLR",
    "SII"
  ),
  
  Lambda = c(
    lambda_agp,
    lambda_nlr,
    lambda_lmr,
    lambda_plr,
    lambda_sii
  )
)



###############################################################
# PLOTS
###############################################################

#Ref value fx
#**************************************************************
# LITERATURE-BASED REFERENCE THRESHOLDS
#**************************************************************

#**************************************************************
# LITERATURE-BASED REFERENCE THRESHOLDS
# Gestational weeks 5–35 only
#**************************************************************

add_reference_lines <- function(primary,
                                secondary = NULL) {
  
  # Primary threshold
  segments(
    x0 = 5,
    y0 = primary,
    x1 = 35,
    y1 = primary,
    col = "firebrick",
    lty = 2,
    lwd = 2
  )
  
  # Secondary threshold
  if (!is.null(secondary)) {
    segments(
      x0 = 5,
      y0 = secondary,
      x1 = 35,
      y1 = secondary,
      col = "purple4",
      lty = 3,
      lwd = 2
    )
  }
}

#**************************************************************
# AGP PLOT
#**************************************************************
# Allow extra space beneath the plot
par(mar = c(9, 5, 4, 2))

plot(
  df_fpr_agp_35$ga_wks,
  df_fpr_agp_35$agp,
  xlab = "Gestational Age (weeks)",
  ylab = "AGP",
  main = "AGP Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  xlim = c(6, 35),
  ylim = c(0, 2),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_5,
  col = "red",
  lwd = 1.5
)

lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_50,
  col = "green",
  lwd = 1.5
)

lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_95,
  col = "blue",
  lwd = 1.5
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

legend(
  "topright",
  legend = c(
    "5th",
    "50th",
    "95th"
  ),
  col = c(
    "red",
    "green",
    "blue"
  ),
  lty = 1,
  lwd = 1.5,
  cex = 0.7
)


# Reference threshold
segments(
  x0 = 6, y0 = 1.0,
  x1 = 35, y1 = 1.0,
  col = "firebrick",
  lty = 2,
  lwd = 2
)



#**************************************************************
# NLR PLOT
#**************************************************************

plot(
  df_fpr_nlr_35$ga_wks,
  df_fpr_nlr_35$nlr_analysis,
  xlab = "Gestational Age (weeks)",
  ylab = "NLR",
  main = "NLR Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  xlim = c(6, 35),
  ylim = c(0, 10),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(
  centile_nlr_bc$ga_wks,
  centile_nlr_bc$centile_5,
  col = "red",
  lwd = 1.5
)

lines(
  centile_nlr_bc$ga_wks,
  centile_nlr_bc$centile_50,
  col = "green",
  lwd = 1.5
)

lines(
  centile_nlr_bc$ga_wks,
  centile_nlr_bc$centile_95,
  col = "blue",
  lwd = 1.5
)

lines(
  emp_nlr$ga_wks,
  emp_nlr$smooth_5,
  col = "red",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_nlr$ga_wks,
  emp_nlr$smooth_50,
  col = "green",
  lty = 2,
  lwd = 1.5
)

lines(
  emp_nlr$ga_wks,
  emp_nlr$smooth_95,
  col = "blue",
  lty = 2,
  lwd = 1.5
)

legend(
  "topright",
  legend = c(
    "5th",
    "50th",
    "95th"
  ),
  col = c(
    "red",
    "green",
    "blue"
  ),
  lty = 1,
  lwd = 1.5,
  cex = 0.7
)



#**************************************************************
# LMR PLOT
#**************************************************************

plot(
  df_fpr_lmr_35$ga_wks,
  df_fpr_lmr_35$lmr_analysis,
  xlab = "Gestational Age (weeks)",
  ylab = "LMR",
  main = "LMR Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  xlim = c(6, 35),
  ylim = c(0, 10),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(
  centile_lmr_bc$ga_wks,
  centile_lmr_bc$centile_5,
  col = "red",
  lwd = 1.5
)

lines(
  centile_lmr_bc$ga_wks,
  centile_lmr_bc$centile_50,
  col = "green",
  lwd = 1.5
)

lines(
  centile_lmr_bc$ga_wks,
  centile_lmr_bc$centile_95,
  col = "blue",
  lwd = 1.5
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

legend(
  "topright",
  legend = c(
    "5th",
    "50th",
    "95th"
  ),
  col = c(
    "red",
    "green",
    "blue"
  ),
  lty = 1,
  lwd = 1.5,
  cex = 0.7
)



#**************************************************************
# PLR PLOT
#**************************************************************

plot(
  df_fpr_plr_35$ga_wks,
  df_fpr_plr_35$plr_analysis,
  xlab = "Gestational Age (weeks)",
  ylab = "PLR",
  main = "PLR Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  xlim = c(6, 35),
  ylim = c(0, 400),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(
  centile_plr_bc$ga_wks,
  centile_plr_bc$centile_5,
  col = "red",
  lwd = 1.5
)

lines(
  centile_plr_bc$ga_wks,
  centile_plr_bc$centile_50,
  col = "green",
  lwd = 1.5
)

lines(
  centile_plr_bc$ga_wks,
  centile_plr_bc$centile_95,
  col = "blue",
  lwd = 1.5
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

legend(
  "topright",
  legend = c(
    "5th",
    "50th",
    "95th"
  ),
  col = c(
    "red",
    "green",
    "blue"
  ),
  lty = 1,
  lwd = 1.5,
  cex = 0.7
)



#**************************************************************
# SII PLOT
#**************************************************************

plot(
  df_fpr_sii_35$ga_wks,
  df_fpr_sii_35$sii_analysis,
  xlab = "Gestational Age (weeks)",
  ylab = "SII",
  main = "SII Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
  xlim = c(6, 35),
  ylim = c(0, 2500),
  pch = 15,
  cex = 0.2,
  col = "grey86"
)

lines(
  centile_sii_bc$ga_wks,
  centile_sii_bc$centile_5,
  col = "red",
  lwd = 1.5
)

lines(
  centile_sii_bc$ga_wks,
  centile_sii_bc$centile_50,
  col = "green",
  lwd = 1.5
)

lines(
  centile_sii_bc$ga_wks,
  centile_sii_bc$centile_95,
  col = "blue",
  lwd = 1.5
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

legend(
  "topright",
  legend = c(
    "5th",
    "50th",
    "95th"
  ),
  col = c(
    "red",
    "green",
    "blue"
  ),
  lty = 1,
  lwd = 1.5,
  cex = 0.7
)

overall_centiles <- tibble(
  Marker = c("AGP", "NLR", "LMR", "PLR", "SII"),
  
  P5 = c(
    quantile(df_fpr_agp_35$agp, 0.05, na.rm = TRUE),
    quantile(df_fpr_nlr_35$nlr_analysis, 0.05, na.rm = TRUE),
    quantile(df_fpr_lmr_35$lmr_analysis, 0.05, na.rm = TRUE),
    quantile(df_fpr_plr_35$plr_analysis, 0.05, na.rm = TRUE),
    quantile(df_fpr_sii_35$sii_analysis, 0.05, na.rm = TRUE)
  ),
  
  P50 = c(
    quantile(df_fpr_agp_35$agp, 0.50, na.rm = TRUE),
    quantile(df_fpr_nlr_35$nlr_analysis, 0.50, na.rm = TRUE),
    quantile(df_fpr_lmr_35$lmr_analysis, 0.50, na.rm = TRUE),
    quantile(df_fpr_plr_35$plr_analysis, 0.50, na.rm = TRUE),
    quantile(df_fpr_sii_35$sii_analysis, 0.50, na.rm = TRUE)
  ),
  
  P95 = c(
    quantile(df_fpr_agp_35$agp, 0.95, na.rm = TRUE),
    quantile(df_fpr_nlr_35$nlr_analysis, 0.95, na.rm = TRUE),
    quantile(df_fpr_lmr_35$lmr_analysis, 0.95, na.rm = TRUE),
    quantile(df_fpr_plr_35$plr_analysis, 0.95, na.rm = TRUE),
    quantile(df_fpr_sii_35$sii_analysis, 0.95, na.rm = TRUE)
  )
)

overall_centiles
library(gt)

overall_centiles %>%
  gt() %>%
  cols_label(
    Marker = "Marker",
    N = "N",
    P5 = "5th percentile",
    P50 = "50th percentile",
    P95 = "95th percentile"
  ) %>%
  fmt_number(
    columns = c(P5, P50, P95),
    decimals = 2
  ) %>%
  tab_header(
    title = "Overall Centiles by Marker",
    subtitle = "Gestational age ≤35 weeks"
  )

#**************************************************************
# OVERALL FPR-DERIVED REFERENCE VALUES
# Normal Approximation - All Biomarkers
#**************************************************************

library(dplyr)
library(tibble)

overall_fpr <- function(model, lambda, marker) {
  
  # Estimated mean and SD from Box-Cox FPR
  mu <- fitted(model, what = "mu")
  sigma <- fitted(model, what = "sigma")
  
  # Overall distribution across gestational ages
  overall_cdf <- function(x) {
    mean(pnorm(x, mean = mu, sd = sigma))
  }
  
  # Calculate overall percentiles on Box-Cox scale
  probs <- c(0.05, 0.50, 0.95)
  
  percentiles <- sapply(probs, function(p) {
    uniroot(
      function(x) overall_cdf(x) - p,
      interval = c(
        min(mu - 10 * sigma),
        max(mu + 10 * sigma)
      )
    )$root
  })
  
  # Back-transform to original biomarker scale
  tibble(
    Marker = marker,
    P5 = inverse_boxcox(percentiles[1], lambda),
    P50 = inverse_boxcox(percentiles[2], lambda),
    P95 = inverse_boxcox(percentiles[3], lambda)
  )
}

overall_centiles <- bind_rows(
  overall_fpr(agp_boxcox, lambda_agp, "AGP"),
  overall_fpr(nlr_boxcox, lambda_nlr, "NLR"),
  overall_fpr(lmr_boxcox, lambda_lmr, "LMR"),
  overall_fpr(plr_boxcox, lambda_plr, "PLR"),
  overall_fpr(sii_boxcox, lambda_sii, "SII")
) %>%
  mutate(across(c(P5, P50, P95), ~ round(.x, 2)))

overall_centiles

loss_summary <- tibble(
  Marker = c("AGP", "NLR", "LMR", "PLR", "SII"),
  
  Full_n = c(
    nrow(df_fpr_agp),
    nrow(df_fpr_nlr),
    nrow(df_fpr_lmr),
    nrow(df_fpr_plr),
    nrow(df_fpr_sii)
  ),
  
  N_le_35 = c(
    nrow(df_fpr_agp_35),
    nrow(df_fpr_nlr_35),
    nrow(df_fpr_lmr_35),
    nrow(df_fpr_plr_35),
    nrow(df_fpr_sii_35)
  )
) %>%
  mutate(
    N_lost = Full_n - N_le_35,
    Percent_lost = round(100 * N_lost / Full_n, 1)
  )

loss_summary

empirical_point_loss <- tibble(
  Marker = c("AGP", "NLR", "LMR", "PLR", "SII"),
  
  Full_empirical_points = c(
    n_distinct(df_fpr_agp$ga_wks),
    n_distinct(df_fpr_nlr$ga_wks),
    n_distinct(df_fpr_lmr$ga_wks),
    n_distinct(df_fpr_plr$ga_wks),
    n_distinct(df_fpr_sii$ga_wks)
  ),
  
  Empirical_points_le_35 = c(
    n_distinct(df_fpr_agp_35$ga_wks),
    n_distinct(df_fpr_nlr_35$ga_wks),
    n_distinct(df_fpr_lmr_35$ga_wks),
    n_distinct(df_fpr_plr_35$ga_wks),
    n_distinct(df_fpr_sii_35$ga_wks)
  )
) %>%
  mutate(
    Points_lost = Full_empirical_points - Empirical_points_le_35,
    Percent_lost = round(
      100 * Points_lost / Full_empirical_points,
      1
    )
  )

View(empirical_point_loss)


#################### PUB READY ############

#**************************************************************
# PUBLICATION-READY AGP REFERENCE CURVES
#**************************************************************

# Save high-resolution TIFF (600 DPI)
tiff(
  filename = "AGP_FPR_Reference_Curves.tiff",
  width = 7,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo"
)

# Plot settings
par(
  mar = c(5, 5, 2, 1),
  mgp = c(2.7, 0.8, 0),
  las = 1,
  cex.axis = 0.9,
  cex.lab = 1.1,
  bty = "l"
)

# Colors
col_5  <- "#D55E00"
col_50 <- "#009E73"
col_95 <- "#0072B2"

# Raw observations
plot(
  df_fpr_agp_35$ga_wks,
  df_fpr_agp_35$agp,
  
  xlim = c(5, 35),
  ylim = c(0, 2),
  
  xlab = "Gestational age (weeks)",
  ylab = "AGP (g/L)",
  
  main = "",
  
  pch = 16,
  cex = 0.25,
  col = adjustcolor("gray60", alpha.f = 0.15)
)

#**************************************************************
# FITTED CENTILES - SOLID
#**************************************************************

lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_5,
  col = col_5,
  lwd = 2.5
)

lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_50,
  col = col_50,
  lwd = 2.5
)

lines(
  centile_agp_bc$ga_wks,
  centile_agp_bc$centile_95,
  col = col_95,
  lwd = 2.5
)

#**************************************************************
# EMPIRICAL CENTILES - DASHED
#**************************************************************

lines(
  emp_agp$ga_wks,
  emp_agp$smooth_5,
  col = col_5,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_agp$ga_wks,
  emp_agp$smooth_50,
  col = col_50,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_agp$ga_wks,
  emp_agp$smooth_95,
  col = col_95,
  lty = 2,
  lwd = 1.5
)

#**************************************************************
# LITERATURE-BASED REFERENCE THRESHOLD
#**************************************************************

segments(
  x0 = 6,
  y0 = 1.0,
  x1 = 35,
  y1 = 1.0,
  col = "black",
  lty = 3,
  lwd = 2
)
# Label directly above reference line
text(
  x = 34.5,
  y = 1.04,
  labels = "AGP = 1.0 g/L",
  col = "black",
  cex = 0.8,
  adj = 1
)

# Compact legend within available plot space
legend(
  "topright",
  legend = c(
    "5th centile",
    "50th centile",
    "95th centile",
    "Fitted",
    "Empirical",
    "Reference threshold"
  ),
  col = c(
    col_5, col_50, col_95,
    "gray35", "gray35", "black"
  ),
  lty = c(1, 1, 1, 1, 2, 3),
  lwd = c(2, 2, 2, 2, 2, 2),
  cex = 0.7,
  bty = "n"
)

dev.off()

#**************************************************************
# PUBLICATION-READY NLR REFERENCE CURVES
# Box-Cox Fractional Polynomial Regression
#**************************************************************

# Save high-resolution TIFF
tiff(
  filename = "NLR_FPR_Reference_Curves.tiff",
  width = 7,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo"
)

# Plot settings
par(
  mar = c(5, 5, 2, 1),
  mgp = c(2.7, 0.8, 0),
  las = 1,
  cex.axis = 0.9,
  cex.lab = 1.1,
  bty = "l"
)

# Publication-friendly colors
col_5  <- "#D55E00"
col_50 <- "#009E73"
col_95 <- "#0072B2"

#**************************************************************
# OBSERVED NLR VALUES
#**************************************************************

plot(
  df_fpr_nlr_35$ga_wks,
  df_fpr_nlr_35$nlr_analysis,
  
  xlim = c(5, 35),
  ylim = c(0, 10),
  
  xlab = "Gestational age (weeks)",
  ylab = "NLR",
  
  main = "",
  
  pch = 16,
  cex = 0.25,
  col = adjustcolor("gray60", alpha.f = 0.15)
)

#**************************************************************
# FITTED CENTILES - SOLID LINES
#**************************************************************

lines(
  centile_nlr_bc$ga_wks,
  centile_nlr_bc$centile_5,
  col = col_5,
  lwd = 2.5
)

lines(
  centile_nlr_bc$ga_wks,
  centile_nlr_bc$centile_50,
  col = col_50,
  lwd = 2.5
)

lines(
  centile_nlr_bc$ga_wks,
  centile_nlr_bc$centile_95,
  col = col_95,
  lwd = 2.5
)

#**************************************************************
# EMPIRICAL CENTILES - DASHED LINES
#**************************************************************

lines(
  emp_nlr$ga_wks,
  emp_nlr$smooth_5,
  col = col_5,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_nlr$ga_wks,
  emp_nlr$smooth_50,
  col = col_50,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_nlr$ga_wks,
  emp_nlr$smooth_95,
  col = col_95,
  lty = 2,
  lwd = 1.5
)

#**************************************************************
# LITERATURE-BASED NLR REFERENCE THRESHOLDS
#**************************************************************

# Primary threshold: NLR >4.5
segments(
  x0 = 6,
  y0 = 4.5,
  x1 = 35,
  y1 = 4.5,
  col = "black",
  lty = 3,
  lwd = 2
)

# Secondary threshold: NLR >6
segments(
  x0 = 6,
  y0 = 6,
  x1 = 35,
  y1 = 6,
  col = "black",
  lty = 4,
  lwd = 2
)

# Reference threshold labels
text(
  x = 34.5,
  y = 4.7,
  labels = "NLR = 4.5",
  col = "black",
  cex = 0.75,
  adj = 1
)

text(
  x = 34.5,
  y = 6.2,
  labels = "NLR = 6.0",
  col = "black",
  cex = 0.75,
  adj = 1
)

#**************************************************************
# LEGEND
#**************************************************************

legend(
  "topright",
  legend = c(
    "5th centile",
    "50th centile",
    "95th centile",
    "Fitted",
    "Empirical",
    "Reference Thresholds"
  ),
  col = c(
    col_5,
    col_50,
    col_95,
    "gray35",
    "gray35",
    "black"
  ),
  lty = c(1, 1, 1, 1, 2),
  lwd = 2,
  cex = 0.7,
  bty = "n"
)

# Close and save
dev.off()

#**************************************************************
# PUBLICATION-READY FPR PLOTTING FUNCTION
#**************************************************************

plot_fpr_publication <- function(data,
                                 outcome,
                                 fitted,
                                 empirical,
                                 biomarker,
                                 y_max,
                                 primary,
                                 secondary = NULL,
                                 direction = ">") {

  # Save 600 DPI TIFF
  tiff(
    filename = paste0(
      biomarker,
      "_FPR_Reference_Curves.tiff"
    ),
    width = 7,
    height = 5.5,
    units = "in",
    res = 600,
    compression = "lzw",
    type = "cairo"
  )

  # Publication settings
  par(
    mar = c(8, 5, 2, 1),
    mgp = c(2.7, 0.8, 0),
    las = 1,
    cex.axis = 0.9,
    cex.lab = 1.1,
    bty = "l"
  )

  # Colors
  col_5  <- "#D55E00"
  col_50 <- "#009E73"
  col_95 <- "#0072B2"

  # Observed data
  plot(
    data$ga_wks,
    data[[outcome]],
    xlim = c(5, 35),
    ylim = c(0, y_max),
    xlab = "Gestational age (weeks)",
    ylab = biomarker,
    main = "",
    pch = 16,
    cex = 0.25,
    col = adjustcolor(
      "gray60",
      alpha.f = 0.15
    )
  )

  # Fitted centiles
  lines(
    fitted$ga_wks,
    fitted$centile_5,
    col = col_5,
    lwd = 2.5
  )

  lines(
    fitted$ga_wks,
    fitted$centile_50,
    col = col_50,
    lwd = 2.5
  )

  lines(
    fitted$ga_wks,
    fitted$centile_95,
    col = col_95,
    lwd = 2.5
  )

  # Empirical centiles
  lines(
    empirical$ga_wks,
    empirical$smooth_5,
    col = col_5,
    lty = 2,
    lwd = 1.5
  )

  lines(
    empirical$ga_wks,
    empirical$smooth_50,
    col = col_50,
    lty = 2,
    lwd = 1.5
  )

  lines(
    empirical$ga_wks,
    empirical$smooth_95,
    col = col_95,
    lty = 2,
    lwd = 1.5
  )

  # Primary reference threshold
  segments(
    x0 = 5,
    y0 = primary,
    x1 = 35,
    y1 = primary,
    col = "black",
    lty = 3,
    lwd = 2
  )

  text(
    x = 34.5,
    y = primary + y_max * 0.025,
    labels = paste0(
      biomarker, " ", direction, primary
    ),
    adj = 1,
    cex = 0.75,
    col = "black"
  )

  # Optional secondary reference threshold
  if (!is.null(secondary)) {

    segments(
      x0 = 5,
      y0 = secondary,
      x1 = 35,
      y1 = secondary,
      col = "gray35",
      lty = 4,
      lwd = 2
    )

    text(
      x = 34.5,
      y = secondary + y_max * 0.025,
      labels = paste0(
        biomarker, " ", direction, secondary
      ),
      adj = 1,
      cex = 0.75,
      col = "gray35"
    )
  }

  # Legend below the plot
  legend(
    "bottom",
    inset = c(0, -0.45),
    legend = c(
      "5th centile",
      "50th centile",
      "95th centile",
      "Fitted",
      "Empirical",
      "Reference thresholds"
    ),
    col = c(
      col_5,
      col_50,
      col_95,
      "gray35",
      "gray35",
      "black"
    ),
    lty = c(1, 1, 1, 1, 2, 3),
    lwd = 2,
    ncol = 3,
    cex = 0.7,
    bty = "n",
    xpd = TRUE
  )

  dev.off()
}

#**************************************************************
# PUBLICATION-READY LMR REFERENCE CURVES
# Box-Cox Fractional Polynomial Regression
#**************************************************************

tiff(
  filename = "LMR_FPR_Reference_Curves.tiff",
  width = 7,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo"
)

# Plot settings
par(
  mar = c(5, 5, 2, 1),
  mgp = c(2.7, 0.8, 0),
  las = 1,
  cex.axis = 0.9,
  cex.lab = 1.1,
  bty = "l"
)

# Colors
col_5  <- "#D55E00"
col_50 <- "#009E73"
col_95 <- "#0072B2"

#**************************************************************
# OBSERVED LMR VALUES
#**************************************************************

plot(
  df_fpr_lmr_35$ga_wks,
  df_fpr_lmr_35$lmr_analysis,
  xlim = c(5, 35),
  ylim = c(0, 10),
  xlab = "Gestational age (weeks)",
  ylab = "LMR",
  main = "",
  pch = 16,
  cex = 0.25,
  col = adjustcolor("gray60", alpha.f = 0.15)
)

#**************************************************************
# FITTED CENTILES - SOLID
#**************************************************************

lines(
  centile_lmr_bc$ga_wks,
  centile_lmr_bc$centile_5,
  col = col_5,
  lwd = 2.5
)

lines(
  centile_lmr_bc$ga_wks,
  centile_lmr_bc$centile_50,
  col = col_50,
  lwd = 2.5
)

lines(
  centile_lmr_bc$ga_wks,
  centile_lmr_bc$centile_95,
  col = col_95,
  lwd = 2.5
)

#**************************************************************
# EMPIRICAL CENTILES - DASHED
#**************************************************************

lines(
  emp_lmr$ga_wks,
  emp_lmr$smooth_5,
  col = col_5,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_lmr$ga_wks,
  emp_lmr$smooth_50,
  col = col_50,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_lmr$ga_wks,
  emp_lmr$smooth_95,
  col = col_95,
  lty = 2,
  lwd = 1.5
)

#**************************************************************
# LITERATURE-BASED LMR REFERENCE THRESHOLD
#**************************************************************

segments(
  x0 = 6,
  y0 = 2,
  x1 = 35,
  y1 = 2,
  col = "black",
  lty = 3,
  lwd = 2
)

# Label threshold
text(
  x = 34.5,
  y = 2.2,
  labels = "LMR <2",
  col = "black",
  cex = 0.8,
  adj = 1
)

#**************************************************************
# LEGEND
#**************************************************************

legend(
  "topright",
  legend = c(
    "5th centile",
    "50th centile",
    "95th centile",
    "Fitted",
    "Empirical",
    "Reference threshold"
  ),
  col = c(
    col_5,
    col_50,
    col_95,
    "gray35",
    "gray35",
    "black"
  ),
  lty = c(1, 1, 1, 1, 2, 3),
  lwd = 2,
  cex = 0.7,
  bty = "n"
)

# Save figure
dev.off()

#**************************************************************
# PUBLICATION-READY PLR REFERENCE CURVES
# Box-Cox Fractional Polynomial Regression
#**************************************************************

tiff(
  filename = "PLR_FPR_Reference_Curves.tiff",
  width = 7,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo"
)

# Plot settings
par(
  mar = c(9, 5, 2, 1),
  mgp = c(2.7, 0.8, 0),
  las = 1,
  cex.axis = 0.9,
  cex.lab = 1.1,
  bty = "l"
)

# Colors
col_5  <- "#D55E00"
col_50 <- "#009E73"
col_95 <- "#0072B2"

#**************************************************************
# OBSERVED PLR VALUES
#**************************************************************

plot(
  df_fpr_plr_35$ga_wks,
  df_fpr_plr_35$plr_analysis,
  xlim = c(5, 35),
  ylim = c(0, 400),
  xlab = "Gestational age (weeks)",
  ylab = "PLR",
  main = "",
  pch = 16,
  cex = 0.25,
  col = adjustcolor("gray60", alpha.f = 0.15)
)

#**************************************************************
# FITTED CENTILES - SOLID LINES
#**************************************************************

lines(
  centile_plr_bc$ga_wks,
  centile_plr_bc$centile_5,
  col = col_5,
  lwd = 2.5
)

lines(
  centile_plr_bc$ga_wks,
  centile_plr_bc$centile_50,
  col = col_50,
  lwd = 2.5
)

lines(
  centile_plr_bc$ga_wks,
  centile_plr_bc$centile_95,
  col = col_95,
  lwd = 2.5
)

#**************************************************************
# EMPIRICAL CENTILES - DASHED LINES
#**************************************************************

lines(
  emp_plr$ga_wks,
  emp_plr$smooth_5,
  col = col_5,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_plr$ga_wks,
  emp_plr$smooth_50,
  col = col_50,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_plr$ga_wks,
  emp_plr$smooth_95,
  col = col_95,
  lty = 2,
  lwd = 1.5
)

#**************************************************************
# LITERATURE-BASED PLR REFERENCE THRESHOLDS
#**************************************************************

# Primary threshold: PLR >250
segments(
  x0 = 6,
  y0 = 250,
  x1 = 35,
  y1 = 250,
  col = "black",
  lty = 3,
  lwd = 2
)

# Secondary threshold: PLR >350
segments(
  x0 = 6,
  y0 = 350,
  x1 = 35,
  y1 = 350,
  col = "gray35",
  lty = 4,
  lwd = 2
)

# Threshold labels
text(
  x = 34.5,
  y = 260,
  labels = "PLR >250",
  col = "black",
  cex = 0.8,
  adj = 1
)

text(
  x = 34.5,
  y = 360,
  labels = "PLR >350",
  col = "gray35",
  cex = 0.8,
  adj = 1
)

#**************************************************************
# LEGEND
#**************************************************************

legend(
  "topright",
  legend = c(
    "5th centile",
    "50th centile",
    "95th centile",
    "Fitted",
    "Empirical",
    "Reference thresholds"
  ),
  col = c(
    col_5,
    col_50,
    col_95,
    "gray35",
    "gray35",
    "black"
  ),
  lty = c(1, 1, 1, 1, 2, 3),
  lwd = 2,
  cex = 0.7,
  bty = "n"
)

# Save figure
dev.off()

#**************************************************************
# PUBLICATION-READY SII REFERENCE CURVES
# Box-Cox Fractional Polynomial Regression
#**************************************************************

tiff(
  filename = "SII_FPR_Reference_Curves.tiff",
  width = 9,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo"
)

par(
  mar = c(5, 5, 2, 10),  # Larger right margin
  mgp = c(2.7, 0.8, 0),
  las = 1,
  cex.axis = 0.9,
  cex.lab = 1.1,
  bty = "l",
  xpd = NA
)

# Colors
col_5  <- "#D55E00"
col_50 <- "#009E73"
col_95 <- "#0072B2"

#**************************************************************
# OBSERVED SII VALUES
#**************************************************************

plot(
  df_fpr_sii_35$ga_wks,
  df_fpr_sii_35$sii_analysis,
  xlim = c(5, 35),
  ylim = c(0, 2500),
  xlab = "Gestational age (weeks)",
  ylab = "SII",
  main = "",
  pch = 16,
  cex = 0.25,
  col = adjustcolor("gray60", alpha.f = 0.15)
)

#**************************************************************
# FITTED CENTILES - SOLID LINES
#**************************************************************

lines(
  centile_sii_bc$ga_wks,
  centile_sii_bc$centile_5,
  col = col_5,
  lwd = 2.5
)

lines(
  centile_sii_bc$ga_wks,
  centile_sii_bc$centile_50,
  col = col_50,
  lwd = 2.5
)

lines(
  centile_sii_bc$ga_wks,
  centile_sii_bc$centile_95,
  col = col_95,
  lwd = 2.5
)

#**************************************************************
# EMPIRICAL CENTILES - DASHED LINES
#**************************************************************

lines(
  emp_sii$ga_wks,
  emp_sii$smooth_5,
  col = col_5,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_sii$ga_wks,
  emp_sii$smooth_50,
  col = col_50,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_sii$ga_wks,
  emp_sii$smooth_95,
  col = col_95,
  lty = 2,
  lwd = 1.5
)

#**************************************************************
# LITERATURE-BASED SII REFERENCE THRESHOLDS
#**************************************************************

# Primary threshold: SII >1000
segments(
  x0 = 6,
  y0 = 1000,
  x1 = 35,
  y1 = 1000,
  col = "black",
  lty = 3,
  lwd = 2
)

# Secondary threshold: SII >2000
segments(
  x0 = 5,
  y0 = 2000,
  x1 = 35,
  y1 = 2000,
  col = "black",
  lty = 4,
  lwd = 2
)

#**************************************************************
# REFERENCE THRESHOLD LABELS
#**************************************************************

text(
  x = 5.5,
  y = 1070,
  labels = "SII >1,000",
  col = "black",
  cex = 0.8,
  adj = 0
)

text(
  x = 5.5,
  y = 2070,
  labels = "SII >2,000",
  col = "black",
  cex = 0.8,
  adj = 0
)

#**************************************************************
# LEGEND - TOP RIGHT
#**************************************************************

legend(
  "topright",
  inset = c(-0.28, 0),
  legend = c(
    "5th centile",
    "50th centile",
    "95th centile",
    "Fitted",
    "Empirical",
    "Reference thresholds"
  ),
  col = c(
    col_5,
    col_50,
    col_95,
    "gray35",
    "gray35",
    "black"
  ),
  lty = c(1, 1, 1, 1, 2, 3),
  lwd = 2,
  cex = 0.7,
  bty = "n",
  xpd = NA
)

# Save figure
dev.off()

#**************************************************************
# PUBLICATION-READY PLR REFERENCE CURVES
# Box-Cox Fractional Polynomial Regression
#**************************************************************

#**************************************************************
# SAVE HIGH-RESOLUTION TIFF
#**************************************************************

tiff(
  filename = "PLR_FPR_Reference_Curves.tiff",
  width = 9,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo"
)

#**************************************************************
# PLOT SETTINGS
#**************************************************************

par(
  mar = c(5, 5, 2, 10),
  mgp = c(2.7, 0.8, 0),
  las = 1,
  cex.axis = 0.9,
  cex.lab = 1.1,
  bty = "l"
)

# Publication-friendly colors
col_5  <- "#D55E00"
col_50 <- "#009E73"
col_95 <- "#0072B2"

#**************************************************************
# OBSERVED PLR VALUES
#**************************************************************

plot(
  df_fpr_plr_35$ga_wks,
  df_fpr_plr_35$plr_analysis,
  
  xlim = c(5, 35),
  ylim = c(0, 400),
  
  xlab = "Gestational age (weeks)",
  ylab = "PLR",
  
  main = "",
  
  pch = 16,
  cex = 0.25,
  col = adjustcolor("gray60", alpha.f = 0.15)
)

#**************************************************************
# FITTED CENTILES - SOLID LINES
#**************************************************************

lines(
  centile_plr_bc$ga_wks,
  centile_plr_bc$centile_5,
  col = col_5,
  lwd = 2.5
)

lines(
  centile_plr_bc$ga_wks,
  centile_plr_bc$centile_50,
  col = col_50,
  lwd = 2.5
)

lines(
  centile_plr_bc$ga_wks,
  centile_plr_bc$centile_95,
  col = col_95,
  lwd = 2.5
)

#**************************************************************
# EMPIRICAL CENTILES - DASHED LINES
#**************************************************************

lines(
  emp_plr$ga_wks,
  emp_plr$smooth_5,
  col = col_5,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_plr$ga_wks,
  emp_plr$smooth_50,
  col = col_50,
  lty = 2,
  lwd = 1.5
)

lines(
  emp_plr$ga_wks,
  emp_plr$smooth_95,
  col = col_95,
  lty = 2,
  lwd = 1.5
)

#**************************************************************
# LITERATURE-BASED PLR REFERENCE THRESHOLDS
#**************************************************************

# Primary threshold: PLR >250
segments(
  x0 = 5,
  y0 = 250,
  x1 = 35,
  y1 = 250,
  col = "black",
  lty = 3,
  lwd = 2
)

# Secondary threshold: PLR >350
segments(
  x0 = 5,
  y0 = 350,
  x1 = 35,
  y1 = 350,
  col = "black",
  lty = 4,
  lwd = 2
)

#**************************************************************
# REFERENCE THRESHOLD LABELS
#**************************************************************

# Primary threshold label
text(
  x = 34.5,
  y = 260,
  labels = "PLR >250",
  col = "black",
  cex = 0.8,
  adj = 1
)

# Secondary threshold label
text(
  x = 34.5,
  y = 360,
  labels = "PLR >350",
  col = "black",
  cex = 0.8,
  adj = 1
)

#**************************************************************
# LEGEND - JUST OUTSIDE TOP RIGHT
#**************************************************************

legend(
  "topright",
  inset = c(-0.28, 0),
  legend = c(
    "5th centile",
    "50th centile",
    "95th centile",
    "Fitted",
    "Empirical",
    "Reference thresholds"
  ),
  col = c(
    col_5,
    col_50,
    col_95,
    "gray35",
    "gray35",
    "black"
  ),
  lty = c(1, 1, 1, 1, 2, 3),
  lwd = 2,
  cex = 0.7,
  bty = "n",
  xpd = NA
)

#**************************************************************
# SAVE COMPLETED FIGURE
#**************************************************************

dev.off()


#**************************************************************
# EXPORT ALL REFERENCE VALUES AS CSV FILES
#**************************************************************

# 1. Gestational age-specific reference values
write.csv(
  reference_values,
  "FPR_Weekly_Reference_Values.csv",
  row.names = FALSE
)

# 2. Trimester-specific reference values (AUC)
write.csv(
  trimester_centiles,
  "FPR_Trimester_Reference_Values.csv",
  row.names = FALSE
)

# 3. Overall FPR-derived reference values
write.csv(
  overall_centiles,
  "FPR_Overall_Reference_Values.csv",
  row.names = FALSE
)

# Confirm where files were saved
getwd()
