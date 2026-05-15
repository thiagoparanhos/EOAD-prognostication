# EOAD Prognostic Model — Statistical Analysis Code

This repository contains the R code used for the statistical analyses reported in:
"EOAD-Signature Atrophy Predicts Dementia In Early-Onset MCI Due To Alzheimer's: An MRI-based prognostic biomarker"
Paranhos T, Katsumi Y, Brickhouse M, et al. (Neurology, under review)

What this code does:
- Fits Cox proportional hazards models estimating the association between baseline EOAD-signature cortical atrophy and time to progression from MCI to dementia
- Evaluates model discrimination (Harrell's C-index), calibration (calibration slope, O/E ratio), and overall performance (scaled Brier score)
- Performs internal validation using bootstrap resampling (1,000 iterations) with optimism correction
- Generates bootstrap-corrected calibration plots

Requirements:
R packages: survival, rms, dplyr, ggplot2

Data:
Data are not included in this repository. To replicate the analysis with your own dataset, replace the file path in the read.csv() call with the path to your data file. The dataset should include: time-to-event, event indicator, age, sex, CDR sum of boxes, and the EOAD-signature w-score.

Contact
Thiago Paranhos — thiago.s.paranhos@gmail.com
