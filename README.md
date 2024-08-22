# R Scripts and Data for the manuscript 'XXXXX'

This repository contains R scripts and Data required to reproduce the analyses features in the following paper:

Rotunno,R., Crema, E.R. (20XX). XXXXX.

The repository is organised into the following six main directories:

- **data** ... Contains raw data and data preparation scripts.
- **runscripts** ... Contains R scripts for executing all analyses.
- **results_images** ... Contains R images files storing the output of all analyses.
- **post_process_scripts** ... Contains R scripts for processing R image files to generate figures and tables.
- **figures** ... Contains main and supplementary figures.
- **tables** ... Contains main and supplementary tables.

# Dataset

# Bayesian Analyses

Two sets of Bayesian analyses were carried out using presence data (quantile regression) and presence/absence data (binomial regression).Quantile regression analyses used the distance from a single putative origin as predictor of the observed radiocarbon dates, whilst the binomial model used the distance from one or more putative origin and the date associated with the sample as a predictor of the presence/absence of ceramic in each binned context. A total of three quantile regression models and seven binomial regression models were fitted. Origin point(s) used in the binomial models are shown in the table below:

| Model | Origin Sites                                |
|-------|---------------------------------------------|
| m1    | Bir Kiseiba                                 |
| m2    | Adrar Bous 10                               |
| m3    | Ounjougou Ravin de la Mouch                 |
| m4    | Bir Kiseiba + Adrar Bous 10                 |
| m5    | Bir Kiseiba + Ounjougou Ravin de la Mouch   |
| m6    | Adrar Bous 10 + Ounjougou Ravin de la Mouch |


## R Script Pipeline

The repository contains all R scripts as well as R image file containing R objects storing posterior samples for key parameters. The scripts should be executed in order (see number suffix) and each require different intermediate R image files or raw data. The table below lists all R scripts in order of execution, along with their requirement, output, and a brief description of what they do.


Model | R Script | Requires | Generates | Description | Runtime |
| --- | --- | --- | --- | --- | --- |
| NA | `runscripts/map.R` | `data/c14data.csv`;`data/shp/Research_area.*` | `figures/main/map.pdf` | Generates a distribution map of the sample sites | < 1 min |
| Bayesian Quantile Regression | `data/clean_and_prepare_data_quantreg_01.R` | `data/c14data.csv` | `data/input_quantreg.Rdata` | Preprocesses the 14C data (e.g. binning, definition of temporal constraints, etc.) for the quantile regression analysis.   | < 1 min |
| Bayesian Quantile Regression | `/runscripts/quantreg/prepare_quantreg_02b.R` | `data/input_quantreg.Rdata` | `runscripts/quantreg/readyrun_quantreg.Rdata` | Defines the nimble model along with data, constraints, initialisation lists. | < 1 min |
| Bayesian Quantile Regression | `/runscripts/quantreg/run_quantreg_a_03b.R` | `runscripts/quantreg/readyrun_quantreg.Rdata`; `data/input_quantreg.Rdata` | `results_images/quantreg/results_quantreg_a.Rdata` | Fits the quantile regression model using Adrar Bous 10 (ABUS10) as putative origin | ca. 24 hours |
| Bayesian Quantile Regression | `/runscripts/quantreg/run_quantreg_e_03b.R` | `runscripts/quantreg/readyrun_quantreg.Rdata`; `data/input_quantreg.Rdata` | `results_images/quantreg/results_quantreg_e.Rdata` | Fits the quantile regression model using Bir Kiseiba (E798) as putative origin | ca. 24 hours |
| Bayesian Quantile Regression | `/runscripts/quantreg/run_quantreg_o_03b.R` | `runscripts/quantreg/readyrun_quantreg.Rdata`; `data/input_quantreg.Rdata` | `results_images/quantreg/results_quantreg_o.Rdata` | Fits the quantile regression model using Ounjougo Ravin de la Mouche (ONJ_RDM) as putative origin | ca. 24 hours |
| Bayesian Quantile Regression | `post_process_scripts/quantile_reg_posterior_04b.R` | `results_images/quantreg/results_quantreg_a.Rdata`;`results_images/quantreg/results_quantreg_e.Rdata`;`results_images/quantreg/results_quantreg_o.Rdata`;`data/input_binom.Rdata` | `figures/si/quantile_regression.pdf`; `tables/si/quantile_regression_posterior.csv` | Generates a summary table with posteriors statistics and model diagnostic and figure with fitted model. | < 1 min |
| Bayesian Binomial Regression | `data/clean_and_prepare_data_binom_01.R` | `data/c14data.csv` | `data/input_binom.Rdata` | Preprocesses the 14C data (e.g. binning, definition of temporal constraints, etc.) for the binomial regression analysis.. | < 1 min |
| Bayesian Binomial Regression | `/runscripts/binom/prepare_binom_02b.R` | `data/input_binom.Rdata` | `runscripts/binom/readyrun_binom.Rdata` | Defines the nimble model along with data, constraints, initialisation lists. | < 1 min |
| Bayesian Binomial Regression | `/runscripts/binom/run_binom_m1_03a.R` | `runscripts/binom/readyrun_binom.Rdata`; `data/input_binom.Rdata` | `results_images/binom/results_binom_m1.Rdata` | Fits the binomial model m1 | ca. 24 hours |
| Bayesian Binomial Regression | `/runscripts/binom/run_binom_m2_03a.R` | `runscripts/binom/readyrun_binom.Rdata`; `data/input_binom.Rdata` | `results_images/binom/results_binom_m2.Rdata` | Fits the binomial model m2 | ca. 24 hours |
| Bayesian Binomial Regression | `/runscripts/binom/run_binom_m3_03a.R` | `runscripts/binom/readyrun_binom.Rdata`; `data/input_binom.Rdata` | `results_images/binom/results_binom_m3.Rdata` | Fits the binomial model m3 | ca. 24 hours |
| Bayesian Binomial Regression | `/runscripts/binom/run_binom_m4_03a.R` | `runscripts/binom/readyrun_binom.Rdata`; `data/input_binom.Rdata` | `results_images/binom/results_binom_m4.Rdata` | Fits the binomial model m4 | ca. 24 hours |
| Bayesian Binomial Regression | `/runscripts/binom/run_binom_m5_03a.R` | `runscripts/binom/readyrun_binom.Rdata`; `data/input_binom.Rdata` | `results_images/binom/results_binom_m5.Rdata` | Fits the binomial model m5 | ca. 24 hours |
| Bayesian Binomial Regression | `/runscripts/binom/run_binom_m6_03a.R` | `runscripts/binom/readyrun_binom.Rdata`; `data/input_binom.Rdata` | `results_images/binom/results_binom_m6.Rdata` | Fits the binomial model m6 | ca. 24 hours |
| Bayesian Binomial Regression | `/runscripts/binom/run_binom_m7_03a.R` | `runscripts/binom/readyrun_binom.Rdata`; `data/input_binom.Rdata` | `results_images/binom/results_binom_m7.Rdata` | Fits the binomial model m7 | ca. 24 hours |
| Bayesian Binomial Regression | `post_process_scripts/binomial_prediction_04a_1.R` | `/runscripts/binom/run_binom_m1_03a.R` ~ `/runscripts/binom/run_binom_m7_03a.R`; `data/input_binom.Rdata`; `data/shp/Research_area.*` | `figures/main/pred_m5.pdf`; `figures/main/pred_m7.pdf`; `figures/si/pred_m1.pdf`; figures/si/pred_m2.pdf`; figures/si/pred_m3.pdf`; figures/si/pred_m4.pdf`; figures/si/pred_m6.pdf`;  | Displays the posterior mean probability of pottery presence within the window of analyses for different time-slices | < 1 min |
| Bayesian Binomial Regression | `post_process_scripts/binomial_residuals_04a_2.R` | `/runscripts/binom/run_binom_m1_03a.R` ~ `/runscripts/binom/run_binom_m7_03a.R`; `data/input_binom.Rdata`; `data/shp/Research_area.*` | `figures/si/residuals_lisa_m5_m7_9k_7k.pdf` | Runs a local Getis-Ord G* statistics on the residuals for models m5 and m7 and generates a figure for key time-slices | < 1 min |
| Bayesian Binomial Regression | `post_process_scripts/binomial_posterior_04a_3.R` | `/runscripts/binom/run_binom_m1_03a.R` ~ `/runscripts/binom/run_binom_m7_03a.R`; `data/input_binom.Rdata` | `tables/si/binomial_posterior.csv`; `figures/si/marginal_predictions.pdf` | Generates a summary table with posteriors statistics and model diagnostic and figure with marginal predictions | < 1 min |
| Bayesian Binomial Regression | `post_process_scripts/binomial_model_comparison_04a_4.R` | `/runscripts/binom/run_binom_m1_03a.R` ~ `/runscripts/binom/run_binom_m7_03a.R`; `data/input_binom.Rdata` | `tables/main/table_model_comparison.csv` | Generates summary tables with the WAIC based model comparison | < 1 min |


# Figures and Tables


| Filename | Figure/Table |
| --- | --- |
| `figures/main/maps.pdf` | Fig.1 |
| `figures/main/pred_m5b.pdf` | Fig.2 |
| `figures/main/pred_m7b.pdf` | Fig.3 |
| `figures/si/quantile_regression.pdf | Fig.S1 |
| `figures/si/pred_m1b.pdf | Fig.S2 |
| `figures/si/pred_m2b.pdf | Fig.S3 |
| `figures/si/pred_m3b.pdf | Fig.S4 |
| `figures/si/pred_m4b.pdf | Fig.S5 |
| `figures/si/pred_m6b.pdf | Fig.S6 |
| `figures/si/marginal_predictions.pdf` | Fig.S7 |
| `figures/si/residuals_lisa_m5_m7_9k_7k.pdf` | Fig.S8 |
| `tables/main/table_model_comparison.csv` | Table 1 |
| `tables/main/quantile_regression_posterior.csv` | Table S1 |
| `tables/main/binomial_posterior.csv` | Table S2 |
# R Session Info

## Funding
* Philip Leverhulme Prize (#PLP-2019–304 Awarded to: E.Crema)
  




