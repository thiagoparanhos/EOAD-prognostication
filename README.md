# EOAD Prognostication — Statistical Analysis Code

This repository contains the R code used for the statistical analyses reported in:

> **EOAD-Signature Atrophy Predicts Dementia In Early-Onset MCI Due To Alzheimer's: An MRI-based prognostic biomarker**  
> Paranhos T, Katsumi Y, Brickhouse M, et al. *Neurology* (in press)

---

## Scripts

| File | Description |
|------|-------------|
| `01_sample_size_estimation.R` | Minimum sample size estimation for a survival-based prognostic model using the `pmsampsize` package |
| `02_survival_analysis.R` | Kaplan-Meier survival curves stratified by EOAD-signature atrophy group |
| `03_cox_models.R` | Cox proportional hazards models (base and extended) with hazard ratio tables and AIC |
| `04_model_comparison.R` | Model comparison using Harrell's C-index, Likelihood Ratio Test, and ΔAIC |
| `05_calibration_validation.R` | Internal validation via bootstrap resampling (1,000 iterations) with optimism-corrected discrimination, calibration slope, O/E ratio, scaled Brier score, and calibration plot |

---

## Requirements

Install the following R packages before running the scripts:

```r
install.packages(c("survival", "rms", "survminer", "broom", "dplyr", "ggplot2", "pmsampsize"))
```

---

## Data

Data used in this study are part of the LEADS (Longitudinal Early-Onset Alzheimer's Disease Study) cohort. Requests for data access should be directed to https://ncrad.org/access-samples/available-samples/leads.

The dataset should include the following variables:

| Variable | Description |
|----------|-------------|
| `Time_to_Event` | Time to progression or censoring (months) |
| `Progressed_to_Mild` | Event indicator (1 = progressed to dementia, 0 = censored) |
| `age_at_scan` | Age at baseline MRI (years) |
| `gender` | Sex (factor) |
| `cdrsum` | CDR Sum of Boxes at baseline |
| `EOADsig_wscore` | EOAD-signature cortical atrophy w-score |
| `EOAD_atrophy_1SD` | Atrophy group classification (used in survival curves) |


---

## Contact

Thiago Paranhos — thiago.s.paranhos@gmail.com

