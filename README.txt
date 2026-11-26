Rodrigues, A., Marques, T., Mammola, S. (2026). No cool refuge below: European subterranean spiders will suffer large-scale impacts from climate change. Submitted.

This repository is archived in Zenodo (https://doi.org/10.5281/zenodo.24561290, DOI: 10.5281/zenodo.24561290) and in GitHub (https://github.com/AfonsoBMR/cave_spiders_SDM). Only the Zenodo version contains the necessary data to fully reproduce the analysis.


# Data sources:

- Predictors: See preparing_predictors.html for a detailed list of predictor variables and sources.

- Spiders_db_GEB_2025-UPDATE.xlsx: Occurrence dataset. For metadata and further information, see: Nicolosi, G. et al. An Open Occurrence Dataset for European Subterranean Spiders. Global Ecology and Biogeography (2026).

- Database_Mammola_et_al_2022_Figshare_V3.xlsx: Functional traits dataset. For metadata and further information, see: Mammola, S. et al. A trait database and updated checklist for European subterranean spiders. Sci. Data 9, (2022).


# Code: 

1. preparing_predictors.Rmd: This script documents the workflow used to prepare the predictor grids for modelling from the source rasters. Since its a time-consuming script of heavy spatial operations, we provide the already assembled predictor stacks resulting from this step in the 2.preparing_all_data/predictor_stacks folder (only in the Zenodo repository). As such, execution of preparing_predictors.Rmd and preparing_all_data.Rmd can be skipped, starting the pipeline from step 3.

2. preparing_all_data.Rmd: This script documents assembly of final predictor stacks for modelling, cleaning of the occurrence dataset, and data exploration.

3. range_bagging_v4.R: This script takes species distribution data + covariate data to model species niches and predict their distributions over multiple conditions, using a rarefaction range bagging method. It then computes various range change metrics for each species.

4. traits.R: This script models the relationship between functional traits and climate change vulnerability.

5. results_and_figures.R: This script prints the range bagging results and creates all figures.



# R session parameters used to run the code:

R version 4.4.3 (2025-02-28 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default


locale:
[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8    LC_MONETARY=English_United Kingdom.utf8
[4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.utf8    

time zone: Europe/Lisbon
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] plyr_1.8.9          gghalves_0.1.4      ggbeeswarm_0.7.3    ggpubr_1.0.0        investr_1.4.2       ggimage_0.3.5      
 [7] ggspatial_1.1.10    Matrix_1.7-6        phyloregion_1.0.9   BAT_2.11.1          cowplot_1.2.0       performance_0.17.1 
[13] DHARMa_0.5.0        ggeffects_2.3.2     glmmTMB_1.1.14      missForest_1.6.1    forcats_1.0.1       tibble_3.3.1       
[19] modEvA_3.45         future.apply_1.20.2 future_1.75.0       geometry_0.5.2      nls2_0.3-4          proto_1.0.0        
[25] stringr_1.6.0       ade4_1.7-24         usdm_2.1-7          GGally_2.4.0        tidyr_1.3.2         dplyr_1.2.1        
[31] ggplot2_4.0.3       readxl_1.5.0        dismo_1.3-16        raster_3.6-32       sp_2.2-3            here_1.0.2
[37] terra_1.9-34        knitr_1.52         

loaded via a namespace (and not attached):
  [1] maptpx_1.9-7            fs_2.1.0                sf_1.1-2                lubridate_1.9.5         RColorBrewer_1.1-3     
  [6] insight_1.5.2           doParallel_1.0.17       numDeriv_2016.8-1.1     backports_1.5.1         tools_4.4.3            
 [11] doRNG_1.8.6.3           R6_2.6.1                vegan_2.7-5             clustMixType_0.5-2      mgcv_1.9-4             
 [16] permute_0.9-10          withr_3.0.3             prettyunits_1.2.0       cli_3.6.6               sandwich_3.1-3         
 [21] slam_0.1-56             mvtnorm_1.4-2           S7_0.2.2                randomForest_4.7-1.2    proxy_0.4-29           
 [26] systemfonts_1.3.2       yulab.utils_0.2.4       DEoptim_2.2-8           parallelly_1.48.0       itertools_0.1-3        
 [31] maps_3.4.3              pdist_1.2.1             rstudioapi_0.19.0       optimParallel_1.0-3     generics_0.1.4         
 [36] gridGraphics_0.5-1      combinat_0.0-8          car_3.1-5               abind_1.4-8             lifecycle_1.0.5        
 [41] scatterplot3d_0.3-45    yaml_2.3.12             carData_3.0-6           clusterGeneration_1.3.8 recipes_1.3.3          
 [46] grid_4.4.3              PlotTools_0.4.0         crayon_1.5.3            predicts_0.2-2          lattice_0.23-1         
 [51] magick_2.9.1            pillar_1.11.1           boot_1.3-32             codetools_0.2-20        fastmatch_1.1-8        
 [56] glue_1.8.1              ggiraph_0.9.6           ggfun_0.2.1             fontLiberation_0.1.0    data.table_1.18.4      
 [61] vctrs_0.7.3             Rdpack_2.6.6            cellranger_1.1.0        gtable_0.3.6            datawizard_1.3.1       
 [66] ks_1.15.3               gower_1.0.2             xfun_0.61               rbibutils_2.4.1         prodlim_2026.03.11     
 [71] pracma_2.4.6            coda_0.19-4.1           reformulas_0.4.4        survival_3.8-9          timeDate_4052.112      
 [76] iterators_1.0.14        hardhat_1.4.3           units_1.0-1             lava_1.9.2              ipred_0.9-15           
 [81] nlme_3.1-170            smoothr_1.3.0           bit64_4.8.2             fontquiver_0.2.1        progress_1.2.3         
 [86] rprojroot_2.1.1         TMB_1.9.23              TreeTools_2.4.0         vipor_0.4.7             KernSmooth_2.23-27     
 [91] otel_0.2.0              rpart_4.1.27            colorspace_2.1-3        DBI_1.3.0               nnet_7.3-21            
 [96] phangorn_2.12.1         mnormt_2.1.2            tidyselect_1.2.1        bit_4.6.0               compiler_4.4.3         
[101] expm_1.0-0              fontBitstreamVera_0.1.1 bayestestR_0.18.1       scales_1.4.0            classInt_0.4-11        
[106] quadprog_1.5-8          rappdirs_0.3.4          palmerpenguins_0.1.1    digest_0.6.39           minqa_1.2.8            
[111] rmarkdown_2.32          htmltools_0.5.9         pkgconfig_2.0.3         lme4_2.0-6              fastmap_1.2.0          
[116] rlang_1.3.0             htmlwidgets_1.6.4       farver_2.1.2            zoo_1.9-0               mclust_6.1.3           
[121] ModelMetrics_1.2.2.2    magrittr_2.0.5          Formula_1.2-6           ggplotify_0.1.3         parameters_0.29.2      
[126] Rcpp_1.1.2              ape_5.8-1               gdtools_0.5.1           stringi_1.8.9           pROC_1.19.0.1          
[131] MASS_7.3-66             hypervolume_3.1.6       ggstats_0.13.0          parallel_4.4.3          listenv_1.0.0          
[136] splines_4.4.3           hms_1.1.4               igraph_2.3.3            ranger_0.18.0           fastcluster_1.3.0      
[141] ggsignif_0.6.4          rngtools_1.5.2          effectsize_1.0.3        reshape2_1.4.5          stats4_4.4.3           
[146] magic_1.6-1             evaluate_1.0.5          nloptr_2.2.1            foreach_1.5.2           purrr_1.2.2            
[151] broom_1.0.13            phytools_2.5-2          e1071_1.7-17            rstatix_1.1.0           class_7.3-24           
[156] beeswarm_0.4.0          cluster_2.1.8.3         timechange_0.4.0        globals_0.19.1          caret_7.0-1 