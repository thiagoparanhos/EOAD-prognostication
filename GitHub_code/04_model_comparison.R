# ============================================================
# Model comparison: Harrell's C, Likelihood Ratio Test, ΔAIC
# Paranhos et al. - EOAD Prognostication
# ============================================================

library(survival)
library(dplyr)

# Read data
# Replace with the path to your dataset
data <- read.csv("/Users/tp015/Partners HealthCare Dropbox/FTD Unit/Dickerson lab/Papers_chapters/Papers_submitted/Paranhos_EOAD_Prognostication/Analysis/R_Models/MCI_to_Mild/wscores_inverted_final_with_averages.csv", check.names = TRUE)

# Remove missing values for model variables
data <- data %>%
  select(Time_to_Event, Progressed_to_Mild,
         age_at_scan, gender, cdrsum, EOADsig_wscore) %>%
  filter(complete.cases(.))

data$gender <- as.factor(data$gender)
data$Progressed_to_Mild <- as.integer(data$Progressed_to_Mild)

# -----------------------------
# Fit nested models
# -----------------------------
cox_base <- coxph(
  Surv(Time_to_Event, Progressed_to_Mild) ~ age_at_scan + gender + cdrsum,
  data = data
)

cox_ext <- coxph(
  Surv(Time_to_Event, Progressed_to_Mild) ~ age_at_scan + gender + cdrsum + EOADsig_wscore,
  data = data
)

# -----------------------------
# Harrell's C-index
# -----------------------------
get_c <- function(fit) {
  s <- summary(fit)$concordance
  list(C = unname(s[1]), SE = unname(s[2]))
}

c_base <- get_c(cox_base)
c_ext  <- get_c(cox_ext)

cat("Harrell's C (base model)    :", round(c_base$C, 3), " (SE", round(c_base$SE, 3), ")\n")
cat("Harrell's C (extended model):", round(c_ext$C,  3), " (SE", round(c_ext$SE,  3), ")\n")
cat("Delta C (extended - base)   :", round(c_ext$C - c_base$C, 3), "\n\n")

# -----------------------------
# Likelihood Ratio Test
# -----------------------------
lrt <- anova(cox_base, cox_ext, test = "LRT")

lrt_chisq <- lrt$Chisq[2]
lrt_df    <- lrt$Df[2]
lrt_p     <- lrt$`Pr(>|Chi|)`[2]

cat("LRT (base vs extended): Chi-square =", round(lrt_chisq, 3),
    " df =", lrt_df,
    " p =", ifelse(lrt_p < 0.001, "< 0.001", format(round(lrt_p, 3), nsmall = 3)), "\n\n")

# -----------------------------
# AIC and ΔAIC
# -----------------------------
aic_base  <- AIC(cox_base)
aic_ext   <- AIC(cox_ext)
delta_aic <- aic_ext - aic_base

cat("AIC (base model)    :", round(aic_base, 2), "\n")
cat("AIC (extended model):", round(aic_ext,  2), "\n")
cat("Delta AIC (extended - base):", round(delta_aic, 2),
    ifelse(delta_aic < 0, " (favors extended model)", " (favors base model)"), "\n\n")
