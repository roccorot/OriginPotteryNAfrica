# R Scripts and Data for the manuscript 'Bayesian analyses of radiocarbon dates suggest multiple origins of ceramic technology in Early Holocene Africa.'

This repository contains R scripts and Data required to reproduce the analysis featured in the following manuscript:

Rotunno,R., Crema, E.R. (2025). Bayesian analyses of radiocarbon dates suggest multiple origins of ceramic technology in Early Holocene Africa.

The repository is organised into the following six main directories:

- **data** ... Contains raw data and data preparation scripts.
- **runscripts** ... Contains R scripts for executing all analyses.
- **results_images** ... Contains R image files storing the output of all analyses.
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
| Bayesian Quantile Regression | `/runscripts/quantreg/run_quantreg_a_03b.R` | `runscripts/quantreg/readyrun_quantreg.Rdata`; `data/input_quantreg.Rdata` | `results_images/quantreg/results_quantreg_a.Rdata` | Fits the quantile regression model using Adrar Bous 10 (ABUS10) as putative origin | ca. 5 hours |
| Bayesian Quantile Regression | `/runscripts/quantreg/run_quantreg_e_03b.R` | `runscripts/quantreg/readyrun_quantreg.Rdata`; `data/input_quantreg.Rdata` | `results_images/quantreg/results_quantreg_e.Rdata` | Fits the quantile regression model using Bir Kiseiba (E798) as putative origin | ca. 5 hours |
| Bayesian Quantile Regression | `/runscripts/quantreg/run_quantreg_o_03b.R` | `runscripts/quantreg/readyrun_quantreg.Rdata`; `data/input_quantreg.Rdata` | `results_images/quantreg/results_quantreg_o.Rdata` | Fits the quantile regression model using Ounjougo Ravin de la Mouche (ONJ_RDM) as putative origin | ca. 5 hours |
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
| `figures/si/quantile_regression.pdf` | Fig.S1 |
| `figures/si/pred_m1b.pdf` | Fig.S2 |
| `figures/si/pred_m2b.pdf` | Fig.S3 |
| `figures/si/pred_m3b.pdf` | Fig.S4 |
| `figures/si/pred_m4b.pdf` | Fig.S5 |
| `figures/si/pred_m6b.pdf` | Fig.S6 |
| `figures/si/marginal_predictions.pdf` | Fig.S7 |
| `figures/si/residuals_lisa_m5_m7_9k_7k.pdf` | Fig.S8 |
| `tables/main/table_model_comparison.csv` | Table 1 |
| `tables/main/quantile_regression_posterior.csv` | Table S1 |
| `tables/main/binomial_posterior.csv` | Table S2 |

# R Session Info
```
ttached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] viridis_0.6.5       viridisLite_0.4.2   tidyr_1.3.0        
 [4] terra_1.7-71        sfdep_0.2.4         spdep_1.3-5        
 [7] spData_2.3.1        sf_1.0-15           rnaturalearth_1.0.1
[10] raster_3.6-26       sp_2.1-3            purrr_1.0.2        
[13] patchwork_1.2.0     nimbleCarbon_0.2.5  nimble_1.2.0       
[16] ggplot2_3.5.1       coda_0.19-4.1       rcarbon_1.5.1      
[19] here_1.0.1         

loaded via a namespace (and not attached):
 [1] httr_1.4.7             jsonlite_1.8.8         splines_4.2.2         
 [4] foreach_1.5.2          spatstat.geom_3.3-2    numDeriv_2016.8-1.1   
 [7] pillar_1.9.0           lattice_0.20-45        glue_1.7.0            
[10] doSNOW_1.0.20          polyclip_1.10-6        colorspace_2.1-0      
[13] Matrix_1.6-1.1         spatstat.sparse_3.1-0  pkgconfig_2.0.3       
[16] s2_1.1.6               scales_1.3.0           spatstat.explore_3.3-2
[19] snow_0.4-4             tensor_1.5             spatstat.utils_3.1-0  
[22] pracma_2.4.4           tibble_3.2.1           proxy_0.4-27          
[25] spatstat.model_3.3-1   mgcv_1.8-41            generics_0.1.3        
[28] spatstat.random_3.3-1  withr_3.0.0            cli_3.6.3             
[31] magrittr_2.0.3         deldir_1.0-9           fansi_1.0.6           
[34] nlme_3.1-162           class_7.3-21           tools_4.2.2           
[37] lifecycle_1.0.4        munsell_0.5.0          spatstat.univar_3.0-0 
[40] compiler_4.2.2         e1071_1.7-14           rlang_1.1.4           
[43] classInt_0.4-10        units_0.8-5            grid_4.2.2            
[46] iterators_1.0.14       goftest_1.2-3          igraph_2.0.3          
[49] boot_1.3-28.1          spatstat.linnet_3.2-1  wk_0.9.1              
[52] gtable_0.3.4           codetools_0.2-19       abind_1.4-5           
[55] DBI_1.2.1              R6_2.5.1               gridExtra_2.3         
[58] knitr_1.45             dplyr_1.1.3            utf8_1.2.4            
[61] rprojroot_2.0.3        KernSmooth_2.23-20     spatstat.data_3.1-2   
[64] spatstat_3.1-1         Rcpp_1.0.13            vctrs_0.6.5           
[67] rpart_4.1.23           tidyselect_1.2.0       xfun_0.41  
```
## Funding
* Philip Leverhulme Prize (#PLP-2019–304 Awarded to: E.Crema)
* EHSCAN-Exploring Early Holocene Saharan Cultural Adaptation and social Networks through socio-ecological inferential modelling. (Engineering and Physical Sciences Research Council (EPSRC) #EP/Y028430/1; Awarded to: E.Crema & R.Rotunno)
  




