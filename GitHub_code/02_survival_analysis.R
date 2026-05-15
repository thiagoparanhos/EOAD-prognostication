# ============================================================
# Kaplan-Meier survival curves by EOAD atrophy group
# Paranhos et al. - EOAD Prognostication
# ============================================================

library(survival)
library(survminer)

# Read data
# Replace with the path to your dataset
df <- read.csv("path/to/your/data.csv")

# Fit survival model stratified by EOAD atrophy group (1SD classification)
fit <- survfit(Surv(Time_to_Event, Progressed_to_Mild) ~ EOAD_atrophy_1SD, data = df)

# Plot Kaplan-Meier curves (displayed as cumulative event probability)
ggsurv <- ggsurvplot(
  fit,
  linetype      = c("solid", "solid", "dashed"),
  conf.int      = FALSE,
  palette       = c("grey", "black", "black"),
  fun           = "event",
  ylim          = c(0, 1),
  xlab          = "Time (months)",
  ylab          = "Dementia (proportion)",
  font.x        = 22,
  font.y        = 22,
  font.tickslab = 18,
  censor.size   = 0,
  size          = 1.0
)

ggsurv$plot <- ggsurv$plot +
  theme(legend.position = "none")

print(ggsurv)
