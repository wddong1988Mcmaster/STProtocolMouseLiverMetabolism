### STProtocolMouseLiverMetabolism
An Integrated Pipeline for Cell-Type Annotation, Metabolic Pathway Profiling, and Spatial Communication Analysis in the Liver using Spatial Transcriptomics

##### Overview for downstream analysis #####
## This tutorial will cover the following tasks:

## 1. Quality Control (QC) and normalization
## 2. Dimension reduction
## 3. Clustering
## 4. Non-linear dimensional reduction (UMAP/t-SNE)
## 5. Identification of cluster biomarkers
## 6. Cell type annotation 
#      a. Automatic annotation using GPT; 
#      b. Annotation via deconvolution using a scRNA-seq reference 
#      c. Manual annotation
## 7. DEG analysis, pathway enrichment analysis and Spatially Variable Genes (SVGs) analysis
## 8. Integrative analysis across multiple samples or conditions
## 9. Pseudobulk analysis
## 10. Quantification of cell type composition 
## 11. Cellular communication
## 12. Metabolic activity analysis
#      a. Metabolic pathway activity
#      b. Metabolic interactions
#      c. Flux balance analysis (FBA)



The example files including FASTQ files, brightfield images, JSON files, HDF5 format files, a spatial folder (spatial/), a scRNA-seq reference for cell type deconvolution, RDS file for pseudobulk analysis and other RDS files for this protocol are available at:  https://www.dropbox.com/scl/fo/occj7x04tbkjrnqfrko1b/AOV196ag11OrAIUtWagY6g8?rlkey=qyzhi6f20x517znqovu4n8q22&st=j4mhkfpm&dl=0.




## > sessionInfo()
R version 4.4.3 (2025-02-28 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/Toronto
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rsvd_1.0.5                  scMetabolism_0.2.1          openai_0.4.1               
 [4] GPTCelltype_1.0.1           celldex_1.16.0              scRNAseq_2.20.0            
 [7] SingleCellExperiment_1.28.1 SingleR_2.8.0               SummarizedExperiment_1.36.0
[10] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3         MatrixGenerics_1.18.1      
[13] matrixStats_1.5.0           lubridate_1.9.4             forcats_1.0.0              
[16] stringr_1.5.1               purrr_1.0.4                 tidyr_1.3.1                
[19] tibble_3.2.1                ggplot2_3.5.2               tidyverse_2.0.0            
[22] PCAtest_0.0.2               export_0.3.0                tximportData_1.34.0        
[25] readr_2.1.5                 tximport_1.34.0             pathview_1.46.0            
[28] dplyr_1.1.4                 GEOquery_2.74.0             SeuratData_0.2.2.9002      
[31] patchwork_1.3.0             loupeR_1.1.4                hdf5r_1.3.12               
[34] Seurat_5.3.0                spatstat.utils_3.1-4        SeuratObject_5.1.0         
[37] sp_2.2-0                    VISION_2.1.0                AUCell_1.28.0              
[40] GSEABase_1.68.0             graph_1.84.1                annotate_1.84.0            
[43] XML_3.99-0.18               AnnotationDbi_1.68.0        IRanges_2.40.1             
[46] S4Vectors_0.44.0            Biobase_2.66.0              BiocGenerics_0.52.0        
[49] GSVA_2.0.7                 

loaded via a namespace (and not attached):
  [1] SpatialExperiment_1.16.0  R.methodsS3_1.8.2         dichromat_2.0-0.1        
  [4] goftest_1.2-3             Biostrings_2.74.1         HDF5Array_1.34.0         
  [7] vctrs_0.6.5               spatstat.random_3.4-1     plumber_1.3.0            
 [10] digest_0.6.37             png_0.1-8                 gypsum_1.2.0             
 [13] ggrepel_0.9.6             deldir_2.0-4              parallelly_1.45.0        
 [16] permute_0.9-7             alabaster.sce_1.6.0       magick_2.8.7             
 [19] MASS_7.3-64               fontLiberation_0.1.0      reshape2_1.4.4           
 [22] httpuv_1.6.16             withr_3.0.2               xfun_0.52                
 [25] survival_3.8-3            memoise_2.0.1             systemfonts_1.2.3        
 [28] ragg_1.4.0                zoo_1.8-14                KEGGgraph_1.66.0         
 [31] pbapply_1.7-2             R.oo_1.27.1               logging_0.10-108         
 [34] KEGGREST_1.46.0           promises_1.3.3            httr_1.4.7               
 [37] restfulr_0.0.15           globals_0.18.0            fitdistrplus_1.2-2       
 [40] rhdf5filters_1.18.1       rhdf5_2.50.2              rstudioapi_0.17.1        
 [43] UCSC.utils_1.2.0          miniUI_0.1.2              generics_0.1.4           
 [46] base64enc_0.1-3           sparsesvd_0.2-2           curl_6.3.0               
 [49] zlibbioc_1.52.0           ScaledMatrix_1.14.0       polyclip_1.10-7          
 [52] ExperimentHub_2.14.0      GenomeInfoDbData_1.2.13   SparseArray_1.6.2        
 [55] xtable_1.8-4              evaluate_1.0.3            S4Arrays_1.6.0           
 [58] BiocFileCache_2.14.0      hms_1.1.3                 webutils_1.2.2           
 [61] irlba_2.3.5.1             filelock_1.0.3            ROCR_1.0-11              
 [64] reticulate_1.42.0         spatstat.data_3.1-6       magrittr_2.0.3           
 [67] lmtest_0.9-40             Rgraphviz_2.50.0          later_1.4.2              
 [70] lattice_0.22-6            spatstat.geom_3.4-1       future.apply_1.20.0      
 [73] iotools_0.3-5             scattermore_1.2           cowplot_1.1.3            
 [76] RcppAnnoy_0.0.22          pillar_1.10.2             nlme_3.1-167             
 [79] compiler_4.4.3            beachmat_2.22.0           RSpectra_0.16-2          
 [82] stringi_1.8.7             tensor_1.5                GenomicAlignments_1.42.0 
 [85] plyr_1.8.9                BiocIO_1.16.0             crayon_1.5.3             
 [88] abind_1.4-8               org.Hs.eg.db_3.20.0       bit_4.6.0                
 [91] codetools_0.2-20          textshaping_1.0.1         BiocSingular_1.22.0      
 [94] openssl_2.3.3             flextable_0.9.9           alabaster.ranges_1.6.0   
 [97] plotly_4.10.4             mime_0.13                 splines_4.4.3            
[100] Rcpp_1.0.14               fastDummies_1.7.5         dbplyr_2.5.0             
[103] sparseMatrixStats_1.18.0  knitr_1.50                blob_1.2.4               
[106] BiocVersion_3.20.0        AnnotationFilter_1.30.0   listenv_0.9.1            
[109] DelayedMatrixStats_1.28.1 openxlsx_4.2.8            Matrix_1.6-5             
[112] statmod_1.5.0             tzdb_0.5.0                pkgconfig_2.0.3          
[115] tools_4.4.3               cachem_1.1.0              stargazer_5.2.3          
[118] RSQLite_2.4.1             viridisLite_0.4.2         DBI_1.2.3                
[121] wordspace_0.2-8           fastmap_1.2.0             rmarkdown_2.29           
[124] scales_1.4.0              grid_4.4.3                pbmcapply_1.5.1          
[127] ica_1.0-3                 Rsamtools_2.22.0          broom_1.0.8              
[130] AnnotationHub_3.14.0      officer_0.6.10            BiocManager_1.30.26      
[133] dotCall64_1.2             alabaster.schemas_1.6.0   RANN_2.6.2               
[136] farver_2.1.2              mgcv_1.9-1                yaml_2.3.10              
[139] rtracklayer_1.66.0        cli_3.6.4                 lifecycle_1.0.4          
[142] rsconnect_1.4.1           askpass_1.2.1             uwot_0.2.3               
[145] backports_1.5.0           BiocParallel_1.40.2       timechange_0.3.0         
[148] gtable_0.3.6              rjson_0.2.23              ggridges_0.5.6           
[151] devEMF_4.5-1              progressr_0.15.1          parallel_4.4.3           
[154] limma_3.62.2              jsonlite_2.0.0            RcppHNSW_0.6.0           
[157] bitops_1.0-9              bit64_4.6.0-1             loe_1.1                  
[160] Rtsne_0.17                vegan_2.7-1               alabaster.matrix_1.6.1   
[163] BiocNeighbors_2.0.1       zip_2.3.3                 alabaster.se_1.6.0       
[166] spatstat.univar_3.1-3     R.utils_2.13.0            alabaster.base_1.6.1     
[169] lazyeval_0.2.2            shiny_1.10.0              rgl_1.3.18               
[172] htmltools_0.5.8.1         rvg_0.3.5                 sctransform_0.4.2        
[175] rappdirs_0.3.3            ensembldb_2.30.0          glue_1.8.0               
[178] spam_2.11-1               httr2_1.1.2               XVector_0.46.0           
[181] gdtools_0.4.2             RCurl_1.98-1.17           mclust_6.1.1             
[184] gridExtra_2.3             igraph_2.1.4              R6_2.6.1                 
[187] GenomicFeatures_1.58.0    cluster_2.1.8             Rhdf5lib_1.28.0          
[190] swagger_5.17.14.1         ProtGenerics_1.38.0       DelayedArray_0.32.0      
[193] tidyselect_1.2.1          xml2_1.3.8                fontBitstreamVera_0.1.1  
[196] future_1.58.0             fastICA_1.2-7             KernSmooth_2.23-26       
[199] fontquiver_0.2.1          data.table_1.17.4         htmlwidgets_1.6.4        
[202] RColorBrewer_1.1-3        rlang_1.1.5               spatstat.sparse_3.1-0    
[205] spatstat.explore_3.4-3    uuid_1.2-1                rentrez_1.2.4
