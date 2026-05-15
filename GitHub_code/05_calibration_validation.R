# ============================================================
# Internal validation and calibration of the extended Cox model
# Paranhos et al. - EOAD Prognostication
# ============================================================

library(survival)
library(rms)
library(dplyr)
library(ggplot2)

# -----------------------------
# 1) Read data
# -----------------------------
# Replace with the path to your dataset
data <- read.csv("path/to/your/data.csv", check.names = TRUE)

# Keep model variables and remove missing values
data2 <- data %>%
  select(Time_to_Event, Progressed_to_Mild,
         age_at_scan, gender, cdrsum, EOADsig_wscore) %>%
  filter(complete.cases(.))

data2$gender <- as.factor(data2$gender)
data2$Progressed_to_Mild <- as.integer(data2$Progressed_to_Mild)

# -----------------------------
# 2) Prediction horizon
# -----------------------------
u <- 24  # months

# -----------------------------
# 3) Fit model
# -----------------------------
# survival package version (used for Brier score and calibration slope)
cox_fit <- coxph(
  Surv(Time_to_Event, Progressed_to_Mild) ~ age_at_scan + gender + cdrsum + EOADsig_wscore,
  data = data2, x = TRUE, y = TRUE
)

# rms package version (used for bootstrap validation)
dd <- datadist(data2)
options(datadist = "dd")

cph_fit <- cph(
  Surv(Time_to_Event, Progressed_to_Mild) ~ age_at_scan + gender + cdrsum + EOADsig_wscore,
  data = data2, x = TRUE, y = TRUE, surv = TRUE, time.inc = u, units = "Month"
)

# -----------------------------
# 4) Apparent C-index
# -----------------------------
conc   <- summary(cox_fit)$concordance
c_index <- unname(conc[1])
c_se    <- unname(conc[2])

cat("\n--- Apparent C-index ---\n")
cat("C-index:", round(c_index, 3),
    " (95% CI:", round(c_index - 1.96 * c_se, 3), "to",
    round(c_index + 1.96 * c_se, 3), ")\n")

# -----------------------------
# 5) Apparent Brier and scaled Brier
# -----------------------------
brier_obj <- brier(cox_fit, times = u)

cat("\n--- Apparent Brier score ---\n")
cat("Brier score at", u, "months:", round(brier_obj$brier[1], 3), "\n")
cat("Scaled Brier score at", u, "months:", round(100 * brier_obj$rsquared[1], 1), "%\n")

# -----------------------------
# 6) Apparent calibration slope
# -----------------------------
lp <- predict(cox_fit, type = "lp")

cal_slope_model <- coxph(Surv(Time_to_Event, Progressed_to_Mild) ~ lp, data = data2)
cal_slope_sum   <- summary(cal_slope_model)

slope    <- cal_slope_sum$coef[1, "coef"]
slope_se <- cal_slope_sum$coef[1, "se(coef)"]

cat("\n--- Apparent calibration slope ---\n")
cat("Calibration slope:", round(slope, 3),
    " (95% CI:", round(slope - 1.96 * slope_se, 3), "to",
    round(slope + 1.96 * slope_se, 3), ")\n")

# -----------------------------
# 7) Observed-to-expected (O/E) ratio at u months
# -----------------------------
surv_func <- Survival(cph_fit)
pred_risk  <- 1 - surv_func(u, lp = lp)

km_fit   <- survfit(Surv(Time_to_Event, Progressed_to_Mild) ~ 1, data = data2)
obs_risk <- 1 - summary(km_fit, times = u, extend = TRUE)$surv
exp_risk <- mean(pred_risk)

cat("\n--- O/E ratio at", u, "months ---\n")
cat("Observed risk:", round(obs_risk, 3), "\n")
cat("Mean predicted risk:", round(exp_risk, 3), "\n")
cat("O/E ratio:", round(obs_risk / exp_risk, 3), "\n")

# -----------------------------
# 8) Bootstrap internal validation (1,000 iterations)
#    Harrell FE. Regression Modeling Strategies. 2nd ed. Springer; 2015.
# -----------------------------
set.seed(123)
val_boot <- validate(cph_fit, method = "boot", B = 1000, dxy = TRUE, u = u)

print(val_boot)

# Optimism-corrected C-index
dxy_corr <- val_boot["Dxy", "index.corrected"]
c_corr   <- dxy_corr / 2 + 0.5

# Optimism-corrected calibration slope
slope_corr <- val_boot["Slope", "index.corrected"]

cat("\n--- Bootstrap-corrected estimates ---\n")
cat("Optimism-corrected C-index          :", round(c_corr,     3), "\n")
cat("Optimism-corrected calibration slope:", round(slope_corr, 3), "\n")
cat("Shrinkage factor                    :", round(slope_corr, 3), "\n")

# -----------------------------
# 9) Bootstrap-corrected calibration plot
# -----------------------------
cal_matrix <- val_boot  # calibrate() output if used separately

# Convert from survival to event probability scale
cal_risk <- as.data.frame(cal_matrix) %>%
  mutate(
    pred_risk  = 1 - pred,
    obs_risk   = 1 - calibrated.corrected,
    lower_risk = pmax(1 - (calibrated.corrected + Upper), 0),
    upper_risk = pmin(1 - (calibrated.corrected + Lower), 1)
  )

ggplot(cal_risk, aes(x = pred_risk, y = obs_risk)) +
  geom_ribbon(aes(ymin = lower_risk, ymax = upper_risk),
              fill = "gray70", alpha = 0.4) +
  geom_line(color = "red", linewidth = 1.4) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    x     = "Predicted 24-month progression probability",
    y     = "Observed 24-month progression probability",
    title = "Bootstrap-corrected calibration (t = 2 years)"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_classic(base_size = 16)
