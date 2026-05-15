# ============================================================
# Sample size estimation for survival-based prognostic model
# Paranhos et al. - EOAD Prognostication
# ============================================================

library(pmsampsize)

# Minimum sample size for a survival model with:
#   - 4 candidate predictors
#   - anticipated R² (Cox-Snell) = 0.20
#   - shrinkage target = 0.90
#   - event rate = 3% per month
#   - prediction horizon = 24 months
#   - mean follow-up = 21.73 months

pmsampsize(
  type        = "s",
  csrsquared  = 0.25,
  shrinkage   = 0.9,
  parameters  = 4,
  rate        = 0.03,
  timepoint   = 24,
  meanfup     = 21.73
)
