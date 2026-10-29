# What Explains Systolic Blood Pressure in U.S. Adults?

Final project for STAT 5543 Applied Regression Analysis, Oklahoma State University (Spring 2025), redesigned in October 2026 around the cross-section workflow in Wooldridge's *Introductory Econometrics*.

The study regresses log systolic blood pressure on age, gender, race, smoking, diabetes and blood lipids for 2,126 fasting adults in NHANES 2011 to 2012. It tests two hypotheses:

- H1: blood pressure rises with age, and the age profile differs for men and women.
- H2: blood lipids (HDL, LDL, triglycerides) are jointly related to blood pressure.

The script follows the coding style and steps used in the course. Parts 1 to 6 are the course project; Part 7 is the extension with econometric checks.

1. Reading the data (Friedewald LDL check; pregnant women and triglycerides of 400 or more excluded)
2. Descriptive stats
3. Initial fit and transformation (log of blood pressure); blood pressure medication left out as a bad control
4. The model: Chow tests by gender to choose the interactions
5. Diagnostics: residual plots, VIF, studentized residuals, Cook's distance
6. Testing and prediction: extra sum of squares F-tests, percent effects, predicted blood pressure by age
7. Extension: RESET, sample comparison, omitted variable bias from obesity, Breusch-Pagan and White tests, robust standard errors, FGLS, median (LAD) regression, and a comparison table across models

The report is `report/STAT-5543-Project-Findings.pdf` (LaTeX source in the same folder).

## Files

| File | Contents |
|---|---|
| `stat5543-final-project-nowshed.R` | All analysis. Writes `output/results.txt` and `output/figure1.png` |
| `data/nhanes.csv` | NHANES extract from cardioStatsUSA (not in the repository; see Data) |
| `output/` | Printed results and figure |
| `report/` | Findings report: PDF, LaTeX source, `elsarticle.cls`, `elsarticle-harv.bst` and the figure |

## Report format

The findings report is written in LaTeX using Elsevier's elsarticle template (elsarticle.cls), the template for the Journal of Financial Economics and other Elsevier journals. The class file is kept in `report/` so the source compiles as is. `elsarticle.cls` and `elsarticle-harv.bst` are copyright Elsevier Ltd and are included unmodified under the LaTeX Project Public License; the MIT license of this repository does not cover them.

## Data

The data is the `nhanes_data` set from the cardioStatsUSA R package (Jaeger, Chen and Muntner, https://github.com/jhs-hwg/cardioStatsUSA, GPL-3.0), which is built from the public NHANES files of the CDC National Center for Health Statistics. The file is about 42 MB, so it is not kept in the repository. To create it, run once in R:

```r
install.packages("remotes")
remotes::install_github("jhs-hwg/cardioStatsUSA")
dir.create("data")
write.csv(cardioStatsUSA::nhanes_data, "data/nhanes.csv", row.names = FALSE)
```

The script keeps the 2011 to 2012 cycle and uses the columns `svy_year`, `bp_sys_mean`, `demo_gender`, `demo_pregnant`, `demo_race`, `demo_age_years`, `bp_med_use`, `cc_smoke`, `cc_diabetes`, `cc_bmi`, `chol_hdl`, `chol_ldl`, `chol_total` and `chol_trig`.

## How to run

Needs R with the packages `car` (3.1 or later, for `vif(type = "predictor")`), `MASS`, `lmtest`, `sandwich` and `quantreg`. Change the `setwd()` line at the top of the script to the folder where this project is kept, then run the script in RStudio or with `Rscript`. Printed results are also saved to `output/results.txt`.

The bootstrap for the median regression uses a fixed seed, so the results are reproducible.

## Code credit and AI assistance

Parts 1 to 6 (the main study) follow the code and conventions from the STAT 5543 class scripts by Dr. Pratyaydipta Rudra, Oklahoma State University. For Part 7 (the extension), I took help from Claude (Anthropic) to write the code. The study design, hypotheses, model choices and interpretation are my own. I have reviewed the code and would appreciate any suggestions for improvement.
