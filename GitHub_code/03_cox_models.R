# ============================================================
# Cox proportional hazards model - base and extended models
# Paranhos et al. - EOAD Prognostication
# ============================================================

library(survival)
library(broom)
library(dplyr)

# Read data
# Replace with the path to your dataset
data <- read.csv("path/to/your/data.csv", check.names = TRUE)

# Remove missing values for model variables
data <- data %>%
  select(Time_to_Event, Progressed_to_Mild,
         age_at_scan, gender, cdrsum, EOADsig_wscore) %>%
  filter(complete.cases(.))

data$gender <- as.factor(data$gender)
data$Progressed_to_Mild <- as.integer(data$Progressed_to_Mild)

# -----------------------------
# Base (clinical) model
# -----------------------------
cox_base <- coxph(
  Surv(Time_to_Event, Progressed_to_Mild) ~ age_at_scan + gender + cdrsum,
  data = data
)

# -----------------------------
# Extended model (base + EOAD-signature atrophy)
# -----------------------------
cox_ext <- coxph(
  Surv(Time_to_Event, Progressed_to_Mild) ~ age_at_scan + gender + cdrsum + EOADsig_wscore,
  data = data
)

# -----------------------------
# Hazard ratio table (extended model)
# -----------------------------
results_ext <- tidy(cox_ext, exponentiate = TRUE, conf.int = TRUE) %>%
  transmute(
    Variable  = term,
    HR        = round(estimate, 2),
    `95% CI`  = paste0("(", round(conf.low, 2), " - ", round(conf.high, 2), ")"),
    `p-value` = ifelse(p.value < 0.001, "< 0.001", format(round(p.value, 3), nsmall = 3))
  )

print(results_ext)

# -----------------------------
# AIC for each model
# -----------------------------
aic_base <- AIC(cox_base)
aic_ext  <- AIC(cox_ext)

cat("\nAIC (base model)    :", round(aic_base, 2), "\n")
cat("AIC (extended model):", round(aic_ext, 2), "\n")
