# CircRNA_pan-cancer
R version 3.6.0 (2019-04-26)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.936  LC_CTYPE=Chinese (Simplified)_China.936   
[3] LC_MONETARY=Chinese (Simplified)_China.936 LC_NUMERIC=C                              
[5] LC_TIME=Chinese (Simplified)_China.936    

attached base packages:
 [1] stats4    parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Rmisc_1.5                   plyr_1.8.6                  psych_2.0.9                 ggExtra_0.9                
 [5] devtools_2.3.2              usethis_1.6.3               ggfortify_0.4.11            DESeq2_1.24.0              
 [9] SummarizedExperiment_1.14.1 DelayedArray_0.10.0         BiocParallel_1.18.1         matrixStats_0.57.0         
[13] Biobase_2.44.0              GenomicRanges_1.36.1        GenomeInfoDb_1.20.0         IRanges_2.18.3             
[17] S4Vectors_0.22.1            BiocGenerics_0.30.0         tidyr_1.1.2                 dplyr_1.0.2                
[21] ggpubr_0.4.0                RColorBrewer_1.1-2          stringr_1.4.0               reshape2_1.4.4             
[25] WGCNA_1.69                  fastcluster_1.1.25          dynamicTreeCut_1.63-1       pheatmap_1.0.12            
[29] e1071_1.7-3                 caret_6.0-86                lattice_0.20-41             varSelRF_0.7-8             
[33] randomForest_4.6-14         pROC_1.16.2                 glmnet_4.0-2                Matrix_1.2-17              
[37] scatterplot3d_0.3-41        Rtsne_0.15                  VennDiagram_1.6.20          futile.logger_1.4.3        
[41] ggplot2_3.3.2              

loaded via a namespace (and not attached):
  [1] tidyselect_1.1.0       RSQLite_2.2.1          AnnotationDbi_1.46.1   htmlwidgets_1.5.2      munsell_0.5.0         
  [6] codetools_0.2-16       preprocessCore_1.46.0  miniUI_0.1.1.1         withr_2.3.0            colorspace_1.4-1      
 [11] knitr_1.30             rstudioapi_0.11        ggsignif_0.6.0         GenomeInfoDbData_1.2.1 mnormt_2.0.2          
 [16] bit64_4.0.5            rprojroot_1.3-2        vctrs_0.3.4            generics_0.1.0         lambda.r_1.2.4        
 [21] ipred_0.9-9            xfun_0.18              R6_2.4.1               doParallel_1.0.15      locfit_1.5-9.4        
 [26] bitops_1.0-6           assertthat_0.2.1       promises_1.1.1         scales_1.1.1           nnet_7.3-12           
 [31] gtable_0.3.0           processx_3.4.4         timeDate_3043.102      rlang_0.4.10           genefilter_1.66.0     
 [36] splines_3.6.0          rstatix_0.6.0          ModelMetrics_1.2.2.2   impute_1.58.0          broom_0.7.1           
 [41] checkmate_2.0.0        yaml_2.2.1             abind_1.4-5            backports_1.1.10       httpuv_1.5.4          
 [46] Hmisc_4.4-1            tools_3.6.0            lava_1.6.8.1           ellipsis_0.3.1         sessioninfo_1.1.1     
 [51] Rcpp_1.0.5             base64enc_0.1-3        zlibbioc_1.30.0        purrr_0.3.4            RCurl_1.98-1.2        
 [56] ps_1.3.4               prettyunits_1.1.1      rpart_4.1-15           haven_2.3.1            cluster_2.0.8         
 [61] fs_1.5.0               magrittr_1.5           data.table_1.13.0      futile.options_1.0.1   openxlsx_4.2.2        
 [66] tmvnsim_1.0-2          pkgload_1.1.0          hms_0.5.3              mime_0.9               xtable_1.8-4          
 [71] XML_3.99-0.3           rio_0.5.16             jpeg_0.1-8.1           readxl_1.3.1           gridExtra_2.3         
 [76] shape_1.4.5            testthat_3.0.2         compiler_3.6.0         tibble_3.0.3           crayon_1.3.4          
 [81] htmltools_0.5.0        later_1.1.0.1          Formula_1.2-3          geneplotter_1.62.0     lubridate_1.7.9       
 [86] DBI_1.1.0              formatR_1.7            MASS_7.3-51.4          car_3.0-10             cli_2.3.0             
 [91] gower_0.2.2            forcats_0.5.0          pkgconfig_2.0.3        foreign_0.8-71         recipes_0.1.15        
 [96] foreach_1.5.0          annotate_1.62.0        XVector_0.24.0         prodlim_2019.11.13     callr_3.5.1           
[101] digest_0.6.25          cellranger_1.1.0       htmlTable_2.1.0        curl_4.3               shiny_1.5.0           
[106] lifecycle_0.2.0        nlme_3.1-139           carData_3.0-4          desc_1.2.0             pillar_1.4.6          
[111] fastmap_1.0.1          pkgbuild_1.1.0         survival_3.2-7         GO.db_3.8.2            glue_1.4.2            
[116] remotes_2.2.0          zip_2.1.1              png_0.1-7              iterators_1.0.12       bit_4.0.4             
[121] class_7.3-15           stringi_1.4.6          blob_1.2.1             latticeExtra_0.6-29    memoise_1.1.0  
