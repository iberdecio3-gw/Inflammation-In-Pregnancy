#site specific fprs
#30 figures (6 sites, 5 biomarkers)

library(gamlss)

ga_values <- seq(6, 35, by = 0.1)

#### --- Kenya --- ####

#configure data
ke_df_fpr_agp_35 <- df_fpr_agp_35 %>%
  filter(site == "Kenya")
ke_df_fpr_nlr_35 <- df_fpr_nlr_35 %>%
  filter(site == "Kenya")
ke_df_fpr_plr_35 <- df_fpr_plr_35 %>%
  filter(site == "Kenya")
ke_df_fpr_lmr_35 <- df_fpr_lmr_35 %>%
  filter(site == "Kenya")
ke_df_fpr_sii_35 <- df_fpr_sii_35 %>%
  filter(site == "Kenya")

#AGP - Kenya 
{

#fpr on agp boxcox
ke_agp_boxcox <- gamlss(
  agp_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = ke_df_fpr_agp_35,
  control = con
)

# Predict centiles on Box-Cox transformed scale
ke_centile_agp_bc <- centiles.pred(
  ke_agp_boxcox,
  data = ke_df_fpr_agp_35,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()


names(ke_centile_agp_bc) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)

ke_centile_agp_bc <- ke_centile_agp_bc %>%
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
}

#NLR - Kenya 
{

#fpr on nlr boxcox
ke_nlr_boxcox <- gamlss(
  nlr_bc ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  family = NO(),
  data = ke_df_fpr_nlr_35,
  control = con
)

# Predict centiles on Box-Cox transformed scale
ke_centile_nlr_bc <- centiles.pred(
  ke_nlr_boxcox,
  data = ke_df_fpr_nlr_35,
  xname = "ga_wks",
  xvalues = ga_values,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame()


names(ke_centile_nlr_bc) <- c(
  "ga_wks",
  "centile_5",
  "centile_50",
  "centile_95"
)

ke_centile_nlr_bc <- ke_centile_nlr_bc %>%
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
}

#PLR - Kenya
{
  
  #fpr on plr boxcox
  ke_plr_boxcox <- gamlss(
    plr_bc ~ fp(ga_wks, 2),
    sigma.fo = ~ fp(ga_wks, 2),
    family = NO(),
    data = ke_df_fpr_plr_35,
    control = con
  )
  
  # Predict centiles on Box-Cox transformed scale
  ke_centile_plr_bc <- centiles.pred(
    ke_plr_boxcox,
    data = ke_df_fpr_plr_35,
    xname = "ga_wks",
    xvalues = ga_values,
    type = "centiles",
    cent = c(5, 50, 95)
  ) %>%
    as.data.frame()
  
  
  names(ke_centile_plr_bc) <- c(
    "ga_wks",
    "centile_5",
    "centile_50",
    "centile_95"
  )
  
  ke_centile_plr_bc <- ke_centile_plr_bc %>%
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
}

#LMR - Kenya
{
  
  #fpr on lmr boxcox
  ke_lmr_boxcox <- gamlss(
    lmr_bc ~ fp(ga_wks, 2),
    sigma.fo = ~ fp(ga_wks, 2),
    family = NO(),
    data = ke_df_fpr_lmr_35,
    control = con
  )
  
  # Predict centiles on Box-Cox transformed scale
  ke_centile_lmr_bc <- centiles.pred(
    ke_lmr_boxcox,
    data = ke_df_fpr_lmr_35,
    xname = "ga_wks",
    xvalues = ga_values,
    type = "centiles",
    cent = c(5, 50, 95)
  ) %>%
    as.data.frame()
  
  
  names(ke_centile_lmr_bc) <- c(
    "ga_wks",
    "centile_5",
    "centile_50",
    "centile_95"
  )
  
  ke_centile_lmr_bc <- ke_centile_lmr_bc %>%
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
}

#SII - Kenya
{
  
  #fpr on sii boxcox
  ke_sii_boxcox <- gamlss(
    sii_bc ~ fp(ga_wks, 2),
    sigma.fo = ~ fp(ga_wks, 2),
    family = NO(),
    data = ke_df_fpr_sii_35,
    control = con
  )
  
  # Predict centiles on Box-Cox transformed scale
  ke_centile_sii_bc <- centiles.pred(
    ke_sii_boxcox,
    data = ke_df_fpr_sii_35,
    xname = "ga_wks",
    xvalues = ga_values,
    type = "centiles",
    cent = c(5, 50, 95)
  ) %>%
    as.data.frame()
  
  
  names(ke_centile_sii_bc) <- c(
    "ga_wks",
    "centile_5",
    "centile_50",
    "centile_95"
  )
  
  ke_centile_sii_bc <- ke_centile_sii_bc %>%
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
}

#-------PLOTS (addition of empirical centiles)-------

#AGP

#Eempirical centiles
{
  
ke_emp_agp <- ke_df_fpr_agp_35 %>%
    filter(!is.na(ga_wks), !is.na(agp)) %>%
    mutate(ga_bin = round(ga_wks * 2) / 2) %>%
    group_by(ga_bin) %>%
    summarise(
      ga_wks = mean(ga_wks),
      empirical_5 = quantile(agp, 0.05),
      empirical_50 = quantile(agp, 0.50),
      empirical_95 = quantile(agp, 0.95),
      n = n(),
      .groups = "drop"
    )
  
  # Smooth empirical centiles
  ke_loess_5 <- loess(
    empirical_5 ~ ga_wks,
    data = ke_emp_agp,
    weights = n,
    span = 0.5
  )
  
  ke_loess_50 <- loess(
    empirical_50 ~ ga_wks,
    data = ke_emp_agp,
    weights = n,
    span = 0.5
  )
  
  ke_loess_95 <- loess(
    empirical_95 ~ ga_wks,
    data = ke_emp_agp,
    weights = n,
    span = 0.5
  )
  
  ke_emp_agp <- ke_emp_agp %>%
    mutate(
      smooth_5 = predict(ke_loess_5, newdata = data.frame(ga_wks = ga_wks)),
      smooth_50 = predict(ke_loess_50, newdata = data.frame(ga_wks = ga_wks)),
      smooth_95 = predict(ke_loess_95, newdata = data.frame(ga_wks = ga_wks))
    )

}
  
#Plot
ke_agp_plot <- function() {
plot(
  ke_df_fpr_agp_35$ga_wks,
  ke_df_fpr_agp_35$agp,
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
# 5th centile
lines(
  ke_centile_agp_bc$ga_wks,
  ke_centile_agp_bc$centile_5,
  col = "red",
  lwd = 1.5
)
# 50th centile
lines(
  ke_centile_agp_bc$ga_wks,
  ke_centile_agp_bc$centile_50,
  col = "green",
  lwd = 1.5
)
# 95th centile
lines(
  ke_centile_agp_bc$ga_wks,
  ke_centile_agp_bc$centile_95,
  col = "blue",
  lwd = 1.5
)

# Empirical centiles - dashed
lines(
  ke_emp_agp$ga_wks,
  ke_emp_agp$smooth_5,
  col = col_5,
  lty = 2,
  lwd = 2.5
)

lines(
  ke_emp_agp$ga_wks,
  ke_emp_agp$smooth_50,
  col = col_50,
  lty = 2,
  lwd = 2.5
)

lines(
  ke_emp_agp$ga_wks,
  ke_emp_agp$smooth_95,
  col = col_95,
  lty = 2,
  lwd = 2.5
)
}



#NLR
#Empirical centiles
{
  ke_emp_nlr <- ke_df_fpr_nlr_35 %>%
    filter(!is.na(ga_wks), !is.na(nlr_analysis)) %>%
    mutate(ga_bin = round(ga_wks * 2) / 2) %>%
    group_by(ga_bin) %>%
    summarise(
      ga_wks = mean(ga_wks),
      empirical_5 = quantile(nlr_analysis, 0.05),
      empirical_50 = quantile(nlr_analysis, 0.50),
      empirical_95 = quantile(nlr_analysis, 0.95),
      n = n(),
      .groups = "drop"
    )
  
  ke_loess_5 <- loess(
    empirical_5 ~ ga_wks,
    data = ke_emp_nlr,
    weights = n,
    span = 0.5
  )
  
  ke_loess_50 <- loess(
    empirical_50 ~ ga_wks,
    data = ke_emp_nlr,
    weights = n,
    span = 0.5
  )
  
  ke_loess_95 <- loess(
    empirical_95 ~ ga_wks,
    data = ke_emp_nlr,
    weights = n,
    span = 0.5
  )
  
  ke_emp_nlr <- ke_emp_nlr %>%
    mutate(
      smooth_5 = predict(ke_loess_5,
                         newdata = data.frame(ga_wks = ga_wks)),
      smooth_50 = predict(ke_loess_50,
                          newdata = data.frame(ga_wks = ga_wks)),
      smooth_95 = predict(ke_loess_95,
                          newdata = data.frame(ga_wks = ga_wks))
    )
}
#Plot
ke_nlr_plot <- function() {
  plot(
    ke_df_fpr_nlr_35$ga_wks,
    ke_df_fpr_nlr_35$nlr_analysis,
    xlab = "Gestational Age (weeks)",
    ylab = "NLR",
    main = "NLR Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
    xlim = c(6, 35),
    ylim = c(0, 10),
    pch = 15,
    cex = 0.2,
    col = "grey86"
  )
  
  # Fitted centiles
  lines(ke_centile_nlr_bc$ga_wks,
        ke_centile_nlr_bc$centile_5,
        col = col_5, lwd = 2.5)
  
  lines(ke_centile_nlr_bc$ga_wks,
        ke_centile_nlr_bc$centile_50,
        col = col_50, lwd = 2.5)
  
  lines(ke_centile_nlr_bc$ga_wks,
        ke_centile_nlr_bc$centile_95,
        col = col_95, lwd = 2.5)
  
  # Empirical centiles - dashed
  lines(ke_emp_nlr$ga_wks,
        ke_emp_nlr$smooth_5,
        col = col_5, lty = 2, lwd = 2.5)
  
  lines(ke_emp_nlr$ga_wks,
        ke_emp_nlr$smooth_50,
        col = col_50, lty = 2, lwd = 2.5)
  
  lines(ke_emp_nlr$ga_wks,
        ke_emp_nlr$smooth_95,
        col = col_95, lty = 2, lwd = 2.5)
}



#PLR
#Empirical centiles
{
  ke_emp_plr <- ke_df_fpr_plr_35 %>%
    filter(!is.na(ga_wks), !is.na(plr_analysis)) %>%
    mutate(ga_bin = round(ga_wks * 2) / 2) %>%
    group_by(ga_bin) %>%
    summarise(
      ga_wks = mean(ga_wks),
      empirical_5 = quantile(plr_analysis, 0.05),
      empirical_50 = quantile(plr_analysis, 0.50),
      empirical_95 = quantile(plr_analysis, 0.95),
      n = n(),
      .groups = "drop"
    )
  
  ke_loess_5 <- loess(
    empirical_5 ~ ga_wks,
    data = ke_emp_plr,
    weights = n,
    span = 0.5
  )
  
  ke_loess_50 <- loess(
    empirical_50 ~ ga_wks,
    data = ke_emp_plr,
    weights = n,
    span = 0.5
  )
  
  ke_loess_95 <- loess(
    empirical_95 ~ ga_wks,
    data = ke_emp_plr,
    weights = n,
    span = 0.5
  )
  
  ke_emp_plr <- ke_emp_plr %>%
    mutate(
      smooth_5 = predict(ke_loess_5,
                         newdata = data.frame(ga_wks = ga_wks)),
      smooth_50 = predict(ke_loess_50,
                          newdata = data.frame(ga_wks = ga_wks)),
      smooth_95 = predict(ke_loess_95,
                          newdata = data.frame(ga_wks = ga_wks))
    )
}
#Plot
ke_plr_plot <- function() {
  plot(
    ke_df_fpr_plr_35$ga_wks,
    ke_df_fpr_plr_35$plr_analysis,
    xlab = "Gestational Age (weeks)",
    ylab = "PLR",
    main = "PLR Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
    xlim = c(6, 35),
    ylim = c(0, 400),
    pch = 15,
    cex = 0.2,
    col = "grey86"
  )
  
  # Fitted centiles
  lines(ke_centile_plr_bc$ga_wks,
        ke_centile_plr_bc$centile_5,
        col = col_5, lwd = 2.5)
  
  lines(ke_centile_plr_bc$ga_wks,
        ke_centile_plr_bc$centile_50,
        col = col_50, lwd = 2.5)
  
  lines(ke_centile_plr_bc$ga_wks,
        ke_centile_plr_bc$centile_95,
        col = col_95, lwd = 2.5)
  
  # Empirical centiles - dashed
  lines(ke_emp_plr$ga_wks,
        ke_emp_plr$smooth_5,
        col = col_5, lty = 2, lwd = 2.5)
  
  lines(ke_emp_plr$ga_wks,
        ke_emp_plr$smooth_50,
        col = col_50, lty = 2, lwd = 2.5)
  
  lines(ke_emp_plr$ga_wks,
        ke_emp_plr$smooth_95,
        col = col_95, lty = 2, lwd = 2.5)
}


#LMR
#Empirical centiles
{
  ke_emp_lmr <- ke_df_fpr_lmr_35 %>%
    filter(!is.na(ga_wks), !is.na(lmr_analysis)) %>%
    mutate(ga_bin = round(ga_wks * 2) / 2) %>%
    group_by(ga_bin) %>%
    summarise(
      ga_wks = mean(ga_wks),
      empirical_5 = quantile(lmr_analysis, 0.05),
      empirical_50 = quantile(lmr_analysis, 0.50),
      empirical_95 = quantile(lmr_analysis, 0.95),
      n = n(),
      .groups = "drop"
    )
  
  ke_loess_5 <- loess(
    empirical_5 ~ ga_wks,
    data = ke_emp_lmr,
    weights = n,
    span = 0.5
  )
  
  ke_loess_50 <- loess(
    empirical_50 ~ ga_wks,
    data = ke_emp_lmr,
    weights = n,
    span = 0.5
  )
  
  ke_loess_95 <- loess(
    empirical_95 ~ ga_wks,
    data = ke_emp_lmr,
    weights = n,
    span = 0.5
  )
  
  ke_emp_lmr <- ke_emp_lmr %>%
    mutate(
      smooth_5 = predict(ke_loess_5,
                         newdata = data.frame(ga_wks = ga_wks)),
      smooth_50 = predict(ke_loess_50,
                          newdata = data.frame(ga_wks = ga_wks)),
      smooth_95 = predict(ke_loess_95,
                          newdata = data.frame(ga_wks = ga_wks))
    )
}
#Plot
ke_lmr_plot <- function() {
  plot(
    ke_df_fpr_lmr_35$ga_wks,
    ke_df_fpr_lmr_35$lmr_analysis,
    xlab = "Gestational Age (weeks)",
    ylab = "LMR",
    main = "LMR Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
    xlim = c(6, 35),
    ylim = c(0, 10),
    pch = 15,
    cex = 0.2,
    col = "grey86"
  )
  
  # Fitted centiles
  lines(ke_centile_lmr_bc$ga_wks,
        ke_centile_lmr_bc$centile_5,
        col = col_5, lwd = 2.5)
  
  lines(ke_centile_lmr_bc$ga_wks,
        ke_centile_lmr_bc$centile_50,
        col = col_50, lwd = 2.5)
  
  lines(ke_centile_lmr_bc$ga_wks,
        ke_centile_lmr_bc$centile_95,
        col = col_95, lwd = 2.5)
  
  # Empirical centiles - dashed
  lines(ke_emp_lmr$ga_wks,
        ke_emp_lmr$smooth_5,
        col = col_5, lty = 2, lwd = 2.5)
  
  lines(ke_emp_lmr$ga_wks,
        ke_emp_lmr$smooth_50,
        col = col_50, lty = 2, lwd = 2.5)
  
  lines(ke_emp_lmr$ga_wks,
        ke_emp_lmr$smooth_95,
        col = col_95, lty = 2, lwd = 2.5)
}



#SII
#Empirical centiles
{
  ke_emp_sii <- ke_df_fpr_sii_35 %>%
    filter(!is.na(ga_wks), !is.na(sii_analysis)) %>%
    mutate(ga_bin = round(ga_wks * 2) / 2) %>%
    group_by(ga_bin) %>%
    summarise(
      ga_wks = mean(ga_wks),
      empirical_5 = quantile(sii_analysis, 0.05),
      empirical_50 = quantile(sii_analysis, 0.50),
      empirical_95 = quantile(sii_analysis, 0.95),
      n = n(),
      .groups = "drop"
    )
  
  ke_loess_5 <- loess(
    empirical_5 ~ ga_wks,
    data = ke_emp_sii,
    weights = n,
    span = 0.5
  )
  
  ke_loess_50 <- loess(
    empirical_50 ~ ga_wks,
    data = ke_emp_sii,
    weights = n,
    span = 0.5
  )
  
  ke_loess_95 <- loess(
    empirical_95 ~ ga_wks,
    data = ke_emp_sii,
    weights = n,
    span = 0.5
  )
  
  ke_emp_sii <- ke_emp_sii %>%
    mutate(
      smooth_5 = predict(ke_loess_5,
                         newdata = data.frame(ga_wks = ga_wks)),
      smooth_50 = predict(ke_loess_50,
                          newdata = data.frame(ga_wks = ga_wks)),
      smooth_95 = predict(ke_loess_95,
                          newdata = data.frame(ga_wks = ga_wks))
    )
}

#Plot
ke_sii_plot <- function() {
  plot(
    ke_df_fpr_sii_35$ga_wks,
    ke_df_fpr_sii_35$sii_analysis,
    xlab = "Gestational Age (weeks)",
    ylab = "SII",
    main = "SII Box-Cox Fractional Polynomial Regression\n(5th, 50th, 95th Centiles)",
    xlim = c(6, 35),
    ylim = c(0, 2500),
    pch = 15,
    cex = 0.2,
    col = "grey86"
  )
  
  # Fitted centiles
  lines(ke_centile_sii_bc$ga_wks,
        ke_centile_sii_bc$centile_5,
        col = col_5, lwd = 2.5)
  
  lines(ke_centile_sii_bc$ga_wks,
        ke_centile_sii_bc$centile_50,
        col = col_50, lwd = 2.5)
  
  lines(ke_centile_sii_bc$ga_wks,
        ke_centile_sii_bc$centile_95,
        col = col_95, lwd = 2.5)
  
  # Empirical centiles - dashed
  lines(ke_emp_sii$ga_wks,
        ke_emp_sii$smooth_5,
        col = col_5, lty = 2, lwd = 2.5)
  
  lines(ke_emp_sii$ga_wks,
        ke_emp_sii$smooth_50,
        col = col_50, lty = 2, lwd = 2.5)
  
  lines(ke_emp_sii$ga_wks,
        ke_emp_sii$smooth_95,
        col = col_95, lty = 2, lwd = 2.5)
}

#*******************************
#* Publication Ready Grid AGP
#*******************************


ke_agp_plot <- function() {
  
  par(
    mar = c(3.5, 4, 1, 1),
    mgp = c(2.2, 0.6, 0),
    las = 1,
    cex.axis = 0.85,
    cex.lab = 0.95,
    bty = "l"
  )
  
  plot(
    ke_df_fpr_agp_35$ga_wks,
    ke_df_fpr_agp_35$agp,
    xlim = c(6, 35),
    ylim = c(0, 2),
    xlab = "Gestational age (weeks)",
    ylab = "AGP (g/L)",
    main = "",
    pch = 16,
    cex = 0.25,
    col = adjustcolor("gray60", alpha.f = 0.30)
  )
  
  # Fitted FPR centiles
  lines(ke_centile_agp_bc$ga_wks, ke_centile_agp_bc$centile_5,
        col = col_5, lwd = 3.5)
  lines(ke_centile_agp_bc$ga_wks, ke_centile_agp_bc$centile_50,
        col = col_50, lwd = 3.5)
  lines(ke_centile_agp_bc$ga_wks, ke_centile_agp_bc$centile_95,
        col = col_95, lwd = 3.5)
  
  # Empirical LOESS centiles
  lines(ke_emp_agp$ga_wks, ke_emp_agp$smooth_5,
        col = col_5, lty = 2, lwd = 2.5)
  lines(ke_emp_agp$ga_wks, ke_emp_agp$smooth_50,
        col = col_50, lty = 2, lwd = 2.5)
  lines(ke_emp_agp$ga_wks, ke_emp_agp$smooth_95,
        col = col_95, lty = 2, lwd = 2.5)
  
  # Literature reference threshold
  segments(6, 1.0, 35, 1.0,
           col = "black", lty = 3, lwd = 2)
  
  text(34.5, 1.04, "AGP = 1.0 g/L",
       adj = 1, cex = 0.75)
}

ke_nlr_plot <- function() {
  
  par(
    mar = c(3.5, 4, 1, 1),
    mgp = c(2.2, 0.6, 0),
    las = 1,
    cex.axis = 0.85,
    cex.lab = 0.95,
    bty = "l"
  )
  
  plot(
    ke_df_fpr_nlr_35$ga_wks,
    ke_df_fpr_nlr_35$nlr_analysis,
    xlim = c(6, 35),
    ylim = c(0, 10),
    xlab = "Gestational age (weeks)",
    ylab = "NLR",
    main = "",
    pch = 16,
    cex = 0.25,
    col = adjustcolor("gray60", alpha.f = 0.30)
  )
  
  # Fitted FPR centiles
  lines(ke_centile_nlr_bc$ga_wks, ke_centile_nlr_bc$centile_5,
        col = col_5, lwd = 3.5)
  lines(ke_centile_nlr_bc$ga_wks, ke_centile_nlr_bc$centile_50,
        col = col_50, lwd = 3.5)
  lines(ke_centile_nlr_bc$ga_wks, ke_centile_nlr_bc$centile_95,
        col = col_95, lwd = 3.5)
  
  # Empirical LOESS centiles
  lines(ke_emp_nlr$ga_wks, ke_emp_nlr$smooth_5,
        col = col_5, lty = 2, lwd = 2.5)
  lines(ke_emp_nlr$ga_wks, ke_emp_nlr$smooth_50,
        col = col_50, lty = 2, lwd = 2.5)
  lines(ke_emp_nlr$ga_wks, ke_emp_nlr$smooth_95,
        col = col_95, lty = 2, lwd = 2.5)
  
  # Literature reference thresholds
  segments(6, 4.5, 35, 4.5,
           col = "black", lty = 3, lwd = 2)
  segments(6, 6.0, 35, 6.0,
           col = "black", lty = 4, lwd = 2)
  
  text(34.5, 4.7, "NLR = 4.5",
       adj = 1, cex = 0.75)
  text(34.5, 6.2, "NLR = 6.0",
       adj = 1, cex = 0.75)
}

ke_lmr_plot <- function() {
  
  par(
    mar = c(3.5, 4, 1, 1),
    mgp = c(2.2, 0.6, 0),
    las = 1,
    cex.axis = 0.85,
    cex.lab = 0.95,
    bty = "l"
  )
  
  plot(
    ke_df_fpr_lmr_35$ga_wks,
    ke_df_fpr_lmr_35$lmr_analysis,
    xlim = c(6, 35),
    ylim = c(0, 20),
    xlab = "Gestational age (weeks)",
    ylab = "LMR",
    main = "",
    pch = 16,
    cex = 0.25,
    col = adjustcolor("gray60", alpha.f = 0.30)
  )
  
  # Fitted FPR centiles
  lines(ke_centile_lmr_bc$ga_wks, ke_centile_lmr_bc$centile_5,
        col = col_5, lwd = 3.5)
  lines(ke_centile_lmr_bc$ga_wks, ke_centile_lmr_bc$centile_50,
        col = col_50, lwd = 3.5)
  lines(ke_centile_lmr_bc$ga_wks, ke_centile_lmr_bc$centile_95,
        col = col_95, lwd = 3.5)
  
  # Empirical LOESS centiles
  lines(ke_emp_lmr$ga_wks, ke_emp_lmr$smooth_5,
        col = col_5, lty = 2, lwd = 2.5)
  lines(ke_emp_lmr$ga_wks, ke_emp_lmr$smooth_50,
        col = col_50, lty = 2, lwd = 2.5)
  lines(ke_emp_lmr$ga_wks, ke_emp_lmr$smooth_95,
        col = col_95, lty = 2, lwd = 2.5)
  
  # Literature reference threshold
  segments(6, 2, 35, 2,
           col = "black", lty = 3, lwd = 2)
  
  text(34.5, 2.2, "LMR <2",
       adj = 1, cex = 0.75)
}

ke_plr_plot <- function() {
  
  par(
    mar = c(3.5, 4, 1, 1),
    mgp = c(2.2, 0.6, 0),
    las = 1,
    cex.axis = 0.85,
    cex.lab = 0.95,
    bty = "l"
  )
  
  plot(
    ke_df_fpr_plr_35$ga_wks,
    ke_df_fpr_plr_35$plr_analysis,
    xlim = c(6, 35),
    ylim = c(0, 400),
    xlab = "Gestational age (weeks)",
    ylab = "PLR",
    main = "",
    pch = 16,
    cex = 0.25,
    col = adjustcolor("gray60", alpha.f = 0.30)
  )
  
  # Fitted FPR centiles
  lines(ke_centile_plr_bc$ga_wks, ke_centile_plr_bc$centile_5,
        col = col_5, lwd = 3.5)
  lines(ke_centile_plr_bc$ga_wks, ke_centile_plr_bc$centile_50,
        col = col_50, lwd = 3.5)
  lines(ke_centile_plr_bc$ga_wks, ke_centile_plr_bc$centile_95,
        col = col_95, lwd = 3.5)
  
  # Empirical LOESS centiles
  lines(ke_emp_plr$ga_wks, ke_emp_plr$smooth_5,
        col = col_5, lty = 2, lwd = 2.5)
  lines(ke_emp_plr$ga_wks, ke_emp_plr$smooth_50,
        col = col_50, lty = 2, lwd = 2.5)
  lines(ke_emp_plr$ga_wks, ke_emp_plr$smooth_95,
        col = col_95, lty = 2, lwd = 2.5)
  
  # Literature reference thresholds
  segments(6, 250, 35, 250,
           col = "black", lty = 3, lwd = 2)
  segments(6, 350, 35, 350,
           col = "black", lty = 4, lwd = 2)
  
  text(34.5, 260, "PLR >250",
       adj = 1, cex = 0.75)
  text(34.5, 360, "PLR >350",
       adj = 1, cex = 0.75)
}

ke_sii_plot <- function() {
  
  par(
    mar = c(3.5, 4, 1, 1),
    mgp = c(2.2, 0.6, 0),
    las = 1,
    cex.axis = 0.85,
    cex.lab = 0.95,
    bty = "l"
  )
  
  plot(
    ke_df_fpr_sii_35$ga_wks,
    ke_df_fpr_sii_35$sii_analysis,
    xlim = c(6, 35),
    ylim = c(0, 2500),
    xlab = "Gestational age (weeks)",
    ylab = "SII",
    main = "",
    pch = 16,
    cex = 0.25,
    col = adjustcolor("gray60", alpha.f = 0.30)
  )
  
  # Fitted FPR centiles
  lines(ke_centile_sii_bc$ga_wks, ke_centile_sii_bc$centile_5,
        col = col_5, lwd = 3.5)
  lines(ke_centile_sii_bc$ga_wks, ke_centile_sii_bc$centile_50,
        col = col_50, lwd = 3.5)
  lines(ke_centile_sii_bc$ga_wks, ke_centile_sii_bc$centile_95,
        col = col_95, lwd = 3.5)
  
  # Empirical LOESS centiles
  lines(ke_emp_sii$ga_wks, ke_emp_sii$smooth_5,
        col = col_5, lty = 2, lwd = 2.5)
  lines(ke_emp_sii$ga_wks, ke_emp_sii$smooth_50,
        col = col_50, lty = 2, lwd = 2.5)
  lines(ke_emp_sii$ga_wks, ke_emp_sii$smooth_95,
        col = col_95, lty = 2, lwd = 2.5)
  
  # Literature reference thresholds
  segments(6, 1000, 35, 1000,
           col = "black", lty = 3, lwd = 2)
  segments(6, 2000, 35, 2000,
           col = "black", lty = 4, lwd = 2)
  
  text(7, 1070, "SII >1,000",
       adj = 0, cex = 0.75)
  text(7, 2070, "SII >2,000",
       adj = 0, cex = 0.75)
}

ke_agp_grob <- as_grob(ke_agp_plot)
ke_nlr_grob <- as_grob(ke_nlr_plot)
ke_lmr_grob <- as_grob(ke_lmr_plot)
ke_plr_grob <- as_grob(ke_plr_plot)
ke_sii_grob <- as_grob(ke_sii_plot)

#grid 
{ 
  ke_combined <- plot_grid(
  ke_agp_grob,
  ke_nlr_grob,
  ke_lmr_grob,
  ke_plr_grob,
  ke_sii_grob,
  NULL,
  labels = c("A", "B", "C", "D", "E", ""),
  ncol = 2,
  align = "hv"
)

kenya_legend <- as_grob(function() {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  
  legend(
    "center",
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
    lwd = c(3.5, 3.5, 3.5, 3.5, 2.5, 2),
    ncol = 3,
    cex = 0.8,
    bty = "n"
  )
})

ke_final_figure <- plot_grid(
  ke_combined,
  kenya_legend,
  ncol = 1,
  rel_heights = c(1, 0.10)
)
}

print(ke_final_figure)

ggsave(
  "Kenya_FPR_Combined_Figure.png",
  plot = ke_final_figure,
  width = 10,
  height = 16,
  units = "in",
  dpi = 600,
  bg = "white"
)
