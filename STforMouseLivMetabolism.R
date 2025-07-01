###### An Integrated Pipeline for Cell-Type Annotation, Metabolic Pathway Profiling, and Spatial Communication Analysis in the Liver using Spatial Transcriptomics.
###### Author: Dongdong Wang (McMaster University) 


##### Overview for downstream analysis #####
## This tutorial will cover the following tasks:

## 1. Quality Control (QC) and normalization
## 2. Dimension reduction
## 3. Clustering
## 4. Non-linear dimensional reduction (UMAP/t-SNE)
## 5. Identification of cluster biomarkers
## 6. Cell type annotation 
##     a. Automatic annotation using GPT; 
##     b. Annotation via deconvolution using a scRNA-seq reference 
##     c. Manual annotation
## 7. DEG analysis, pathway enrichment analysis and Spatially Variable Genes (SVGs) analysis
## 8. Integrative analysis across multiple samples or conditions
## 9. Pseudobulk analysis
## 10. Quantification of cell type composition 
## 11. Cellular communication
## 12. Metabolic activity analysis
##     a. Metabolic pathway activity
##     b. Metabolic interactions
##     c. Flux balance analysis (FBA)



## Environment initialization
rm(list = ls())
options(stringsAsFactor = F)

## Install packages and load necessary packages
install.packages("installr")
install.packages("remotes")
if (!require("hdf5r")) install.packages("hdf5r")
install.packages("Rtools")
install.packages('Seurat')
install.packages("SeuratData")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
install.packages("Signac")
BiocManager::install("SingleR")
BiocManager::install("scRNAseq")
BiocManager::install("celldex")
install.packages('SeuratObject')
install.packages("spatstat.utils")
devtools::install_github("thomasp85/patchwork")
devtools::install_version("spatstat", version = "3.0.5")
remotes::install_github('satijalab/azimuth', ref = 'master')
BiocManager::install("tximport")
install.packages("tidyverse")
BiocManager::install("tximportData")
install.packages("devtools")
devtools::install_github("arleyc/PCAtest")
remotes::install_github("10XGenomics/loupeR")
loupeR::setup()
install.packages("openai")
remotes::install_github("Winnie09/GPTCelltype")
install.packages("data.table")
install.packages("wesanderson")
install.packages("AUCell")
install.packages("GSEABase")
install.packages("GSVA")
install.packages("VISION")
remove.packages("promises")
install.packages("promises")
install.packages("devtools")
remove.packages("VISION")
remove.packages("AUCell")
remove.packages("GSVA")
devtools::install_github("YosefLab/VISION@v2.1.0") #Please note that the version would be v2.1.0
devtools::install_github("wu-yc/scMetabolism")
remotes::install_github("10XGenomics/loupeR")
loupeR::setup()
devtools::install_github('satijalab/seurat-data')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("pathview")
BiocManager::install("tximport")
BiocManager::install("tximportData")
devtools::install_github("arleyc/PCAtest")
BiocManager::install("scRNAseq")
BiocManager::install("celldex")
remotes::install_github("Winnie09/GPTCelltype")
BiocManager::install("GSVA")
BiocManager::install("GSEABase")
BiocManager::install("AUCell")
install.packages("promises")

library(GSVA)
library(GSEABase)
library(AUCell)
library(VISION)
library('SeuratObject')
library("spatstat.utils")
library(spatstat.utils)
library(Seurat)
library("hdf5r")
library(loupeR)
library(patchwork)
library(SeuratData)
library(GEOquery)
library(dplyr)
library(pathview)
library("tximport")
library("readr")
library("tximportData")
library("export")
library(PCAtest)
library(patchwork)
library(tidyverse)
library(SingleR)
library(scRNAseq)
library(celldex)
library(GPTCelltype)
library(openai)
library(scMetabolism)
library(ggplot2)
library(rsvd)

sessionInfo()




############# Load data sets from 10X Space Ranger
############# Working with multiple slices in Seurat
### read it in and perform the same initial normalization
# liver1
liver1=Load10X_Spatial(
  data_dir <- 'C:/YourFolder/Sam19',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "condition1a",
  filter.matrix = TRUE,
  to.upper = FALSE
)

# liver2
liver2=Load10X_Spatial(
  data_dir <- 'C:/YourFolder/Sam10',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "condition2a",
  filter.matrix = TRUE,
  to.upper = FALSE
)

# liver3
liver3=Load10X_Spatial(
  data_dir <- 'C:/YourFolder/Sam18',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "condition1b",
  filter.matrix = TRUE,
  to.upper = FALSE
)

# liver4
liver4=Load10X_Spatial(
  data_dir <- 'C:/YourFolder/Sam11',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "condition2b",
  filter.matrix = TRUE,
  to.upper = FALSE
)

# Change orig.ident name and 
head(liver1[[]])
head(liver2[[]])
head(liver3[[]])
head(liver4[[]])

liver1[[]]$orig.ident='condition1a'
liver2[[]]$orig.ident='condition2a'

liver3[[]]$orig.ident='condition1b'
liver4[[]]$orig.ident='condition2b'

#Rename cell identity classes
levels(liver1)
levels(liver2)
levels(liver3)
levels(liver4)
liver1 <- RenameIdents(liver1, 'SeuratProject' = 'condition1a')
liver2 <- RenameIdents(liver2, 'SeuratProject' = 'condition2a')

liver3 <- RenameIdents(liver3, 'SeuratProject' = 'condition1b')
liver4 <- RenameIdents(liver4, 'SeuratProject' = 'condition2b')

#Name project.name
liver1@project.name = "condition1a_liver"
liver2@project.name = "condition2a_liver"
liver3@project.name = "condition1b_liver"
liver4@project.name = "condition2b_liver"


## 1. Quality Control (QC) and normalization
# for liver1
plot1 <- VlnPlot(liver1, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver1, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv1_gene', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

plot1 <- VlnPlot(liver1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv1_counts', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

# for liver2
plot1 <- VlnPlot(liver2, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver2, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv2_gene', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

plot1 <- VlnPlot(liver2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver2, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv2_counts', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

# for liver3
plot1 <- VlnPlot(liver3, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver3, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv3_gene', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

plot1 <- VlnPlot(liver3, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver3, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv3_counts', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

# for liver4
plot1 <- VlnPlot(liver4, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver4, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv4_gene', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

plot1 <- VlnPlot(liver4, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver4, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv4_counts', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

# We filter cells that have unique feature counts over 7,500 or less than 2000
liver1 <- subset(liver1, subset = nFeature_Spatial > 1000 & nFeature_Spatial < 7500)
liver2 <- subset(liver2, subset = nFeature_Spatial > 1000 & nFeature_Spatial < 7500)
liver3 <- subset(liver3, subset = nFeature_Spatial > 1000 & nFeature_Spatial < 7500)
liver4 <- subset(liver4, subset = nFeature_Spatial > 1000 & nFeature_Spatial < 7500)


# Apply sctransform normalization
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform.
# Replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
options(future.globals.maxSize = 4 * 1024^3)
liver1SCT <- SCTransform(liver1, method = "glmGamPoi", assay = "Spatial", verbose = TRUE)
liver2SCT <- SCTransform(liver2, method = "glmGamPoi", assay = "Spatial", verbose = TRUE)
liver3SCT <- SCTransform(liver3, method = "glmGamPoi", assay = "Spatial", verbose = TRUE)
liver4SCT <- SCTransform(liver4, method = "glmGamPoi", assay = "Spatial", verbose = TRUE)

# for liver2: show counts and features after SCTransform
plot1 <- VlnPlot(liver2SCT, features = "nFeature_SCT", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver2SCT, features = "nFeature_SCT") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv2_gene_SCT', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

plot1 <- VlnPlot(liver2SCT, features = "nCount_SCT", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver2SCT, features = "nCount_SCT") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv2_counts_SCT', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

## save the object at this point
saveRDS(liver1SCT, file = "liver1SCT_Prtc.rds")
saveRDS(liver2SCT, file = "liver2SCT_Prtc.rds")
saveRDS(liver3SCT, file = "liver3SCT_Prtc.rds")
saveRDS(liver4SCT, file = "liver4SCT_Prtc.rds")

## Load the data
liver1SCT <- readRDS("liver1SCT_Prtc.rds")
liver2SCT <- readRDS("liver2SCT_Prtc.rds")
liver3SCT <- readRDS("liver3SCT_Prtc.rds")
liver4SCT <- readRDS("liver4SCT_Prtc.rds")


## 2. Dimension reduction
## 3. Clustering
## 4. Non-linear dimensional reduction (UMAP/t-SNE)
obj=liver1SCT
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 1) #resolution between 0.4-1.2 typically returns good results
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
# Look at cluster IDs of the first 5 cells
head(Idents(obj), 5)
liver1SCT1= obj

# visualize the results of the clustering
p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
p1 + p2
graph2svg(x = NULL, file='2-cluster_liv1', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")


# liver2
obj=liver2SCT
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
# PlotPCA results
p1 = DimPlot(obj, reduction = "pca") + NoLegend()
p1
graph2svg(x = NULL, file='2-1-PCADimecluster_liv2', font = "Arial", cairo = TRUE,   
          width = 4, height = 4, bg = "transparent")
p1 = VizDimLoadings(obj, dims = 1:2, nfeatures = 15, reduction = "pca")
p1
graph2svg(x = NULL, file='2-2-PCADimecluster_liv2', font = "Arial", cairo = TRUE,   
          width = 6, height = 5, bg = "transparent")

# Clustering and UMAP
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 1) #resolution between 0.4-1.2 typically returns good results
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
head(Idents(obj), 5)
liver2SCT1= obj

# visualize the results of the clustering
p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
p1 + p2
graph2svg(x = NULL, file='2-cluster_liv2', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")

# Clustering and t-SNE
obj=liver2SCT
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 1) #resolution between 0.4-1.2 typically returns good results
obj <- RunTSNE(obj, reduction = "pca", dims = 1:30)
head(Idents(obj), 5)
liver2SCT1_tsne= obj

# visualize the results of the clustering
p1 <- DimPlot(obj, reduction = "tsne", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
p1 + p2
graph2svg(x = NULL, file='2-cluster_liv2_tSNE', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")

# liver3
obj=liver3SCT
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 1) #resolution between 0.4-1.2 typically returns good results
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
head(Idents(obj), 5)
liver3SCT1= obj

# visualize the results of the clustering
p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
p1 + p2
graph2svg(x = NULL, file='2-cluster_liv3', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")

# liver4
obj=liver4SCT
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 1) #resolution between 0.4-1.2 typically returns good results
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
head(Idents(obj), 5)
liver4SCT1= obj

# visualize the results of the clustering
p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
p1 + p2
graph2svg(x = NULL, file='2-cluster_liv4', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")



## DimHeatmap(): easy exploration of the primary sources of heterogeneity in a dataset,
# useful when trying to decide which PCs to include for further downstream analyses
# liver1
DimHeatmap(liver1SCT1, dims = 1:6, cells = 500, balanced = TRUE)
graph2svg(x = NULL, file='3-DimHeatmap_liv1', font = "Arial", cairo = TRUE,   
          width = 8, height = 7, bg = "transparent")
# liver2
DimHeatmap(liver2SCT1, dims = 1:6, cells = 500, balanced = TRUE)
graph2svg(x = NULL, file='3-DimHeatmap_liv2', font = "Arial", cairo = TRUE,   
          width = 8, height = 7, bg = "transparent")

# liver3
DimHeatmap(liver3SCT1, dims = 1:6, cells = 500, balanced = TRUE)
graph2svg(x = NULL, file='3-DimHeatmap_liv3', font = "Arial", cairo = TRUE,   
          width = 8, height = 7, bg = "transparent")
# liver4
DimHeatmap(liver4SCT1, dims = 1:6, cells = 500, balanced = TRUE)
graph2svg(x = NULL, file='3-DimHeatmap_liv4', font = "Arial", cairo = TRUE,   
          width = 8, height = 7, bg = "transparent")


## Determine the ‘dimensionality’ of the dataset
# If an ‘elbow’ is around PC9-10, it suggests the majority of true signal is captured in the first 10 PCs.
# liver1
ElbowPlot(liver1SCT1)
graph2svg(x = NULL, file='3_1-ElbowPlot_liv1', font = "Arial", cairo = TRUE,   
          width = 4, height = 3.5, bg = "transparent")
# liver2
ElbowPlot(liver2SCT1)
graph2svg(x = NULL, file='3_1-ElbowPlot_liv2', font = "Arial", cairo = TRUE,   
          width = 4, height = 3.5, bg = "transparent")
# liver3
ElbowPlot(liver3SCT1)
graph2svg(x = NULL, file='3_1-ElbowPlot_liv3', font = "Arial", cairo = TRUE,   
          width = 4, height = 3.5, bg = "transparent")
# liver4
ElbowPlot(liver4SCT1)
graph2svg(x = NULL, file='3_1-ElbowPlot_liv4', font = "Arial", cairo = TRUE,   
          width = 4, height = 3.5, bg = "transparent")

## save the object at this point
saveRDS(liver1SCT1, file = "Liver1_AftCluster_Protc.rds")
saveRDS(liver2SCT1, file = "Liver2_AftCluster_Protc.rds")
saveRDS(liver3SCT1, file = "Liver3_AftCluster_Protc.rds")
saveRDS(liver4SCT1, file = "Liver4_AftCluster_Protc.rds")

## Load the data
liver1SCT1 <- readRDS("Liver1_AftCluster_Protc.rds")
liver2SCT1 <- readRDS("Liver2_AftCluster_Protc.rds")
liver3SCT1 <- readRDS("Liver3_AftCluster_Protc.rds")
liver4SCT1 <- readRDS("Liver4_AftCluster_Protc.rds")



## 5. Identification of cluster biomarkers
# Finding differentially expressed features (cluster biomarker) for liver1
# And generates an expression heatmap for given cells and features (10)
# liver1
obj=liver1SCT1
# find markers for every cluster compared to all remaining cells,
# report only the positive ones
obj.markers <- FindAllMarkers(obj, only.pos = TRUE)
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(obj, features = top10$gene) + NoLegend()
# save it in png file

#### save the file
liver1SCT1=obj
saveRDS(liver1SCT1, file = "liver1_final.rds")

## Load the data
liver1SCT1 <- readRDS("liver1_final.rds")



#### Finding differentially expressed features (cluster biomarkers) for liver2
# liver2
obj=liver2SCT1
# find markers for every cluster compared to all remaining cells,
# report only the positive ones
obj.markers <- FindAllMarkers(obj, only.pos = TRUE)
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(obj, features = top10$gene) + NoLegend()
# save it in png file

## Plot violin fig for biomarkers of clusters
p1=VlnPlot(obj, features = c("Oat", "Cyp2e1", "Slc1a2", "Cyp2c29"), ncol = 2)
p1
graph2svg(x = NULL, file='5_1-VlnClus4&5_liv2', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 5, bg = "transparent")

## Plot in images
p1=SpatialFeaturePlot(obj, features = c("Oat", "Cyp2e1", "Slc1a2", "Cyp2c29"))
p1
graph2svg(x = NULL, file='5_2-SpaClus4&5_liv2', font = "Arial", cairo = TRUE,   
          width = 5, height = 6, bg = "transparent")

## Plot in UMAP
FeaturePlot(obj, features = c("Oat", "Cyp2e1", "Slc1a2", "Cyp2c29"))
graph2svg(x = NULL, file='5_3-FeatinUMAPClus4&5_liv2', font = "Arial", cairo = TRUE,   
          width = 5, height = 6, bg = "transparent")

#### save the file
liver2SCT1=obj
saveRDS(liver2SCT1, file = "liver2_prtc_final.rds")

## Load the data
liver2SCT1 <- readRDS("iver2_prtc_final.rds")



#### Finding differentially expressed features (cluster biomarkers) for liver3
# liver3
obj=liver3SCT1
# find markers for every cluster compared to all remaining cells,
# report only the positive ones
obj.markers <- FindAllMarkers(obj, only.pos = TRUE)
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
liver3_marker=obj.markers
DoHeatmap(obj, features = top10$gene) + NoLegend()
# save it in png file

#### save the file
liver3SCT1=obj
saveRDS(liver3SCT1, file = "liver3_prtc_final.rds")

## Load the data
liver3SCT1 <- readRDS("liver3_prtc_final.rds")



#### Finding differentially expressed features (cluster biomarkers) for liver4
# liver4
obj=liver4SCT1
# find markers for every cluster compared to all remaining cells,
# report only the positive ones
obj.markers <- FindAllMarkers(obj, only.pos = TRUE)
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
liver4_marker=obj.markers
DoHeatmap(obj, features = top10$gene) + NoLegend()
# save it in png file

#### save the file
liver4SCT1=obj
saveRDS(liver4SCT1, file = "liver4_final.rds")

## Load the data
liver4SCT1 <- readRDS("liver4_final.rds")



## 6. Cell type annotation 
##    a. Automatic annotation using GPT

# IMPORTANT! Assign your OpenAI API key. See Vignette for details
Sys.setenv(OPENAI_API_KEY = 'xx')

# Assume you have already run the Seurat pipeline https://satijalab.org/seurat/
# "obj" is the Seurat object; "markers" is the output from FindAllMarkers(obj)
# Cell type annotation by GPT-4

# Find Markers
#liver.merge=PrepSCTFindMarkers(liver.merge, assay = "SCT", verbose = TRUE)
all_markers = FindAllMarkers(liver2SCT1, assay = "SCT")
markers = all_markers

# GPT-4 annotation
res <- gptcelltype(markers, tissuename = 'mouse liver', model = 'gpt-4')
res
write.csv(res, 'Annotation_liv2.csv')
celltype <- c("0: Hepatocytes1", "1: Hepatocytes2", "2: Hepatocytes3", "3: Kuffer cells", "4: Hepatocytes4","5: Hepatocytes5",
              "6: Hepatocytes6", "7: Hepatic Stellate Cells", "8: Hepatocytes7", "9: Macrophages", "10: Erythrocytes",
              "11: B cells", "12: Immune Cells", "13: Hepatocytes8")
liver2SCT2=liver2SCT1
names(celltype) <- levels(liver2SCT2)
liver2SCT2 <- RenameIdents(liver2SCT2, celltype)

# Visualize cell type annotation on UMAP
DimPlot(liver2SCT2, label = TRUE)
graph2svg(x = NULL, file='10-3-Liv2_DimPlot_labeling_scRNA', font = "Arial", cairo = TRUE,   
          width = 7, height = 5, bg = "transparent")
head(liver.merge1[[]])

### visualized in a SpatialDimPlot
SpatialDimPlot(liver2SCT2, label = TRUE, label.size = 4)
graph2svg(x = NULL, file='9-liv2_SpatialDimPlot', font = "Arial", cairo = TRUE,   
          width = 8, height = 6, bg = "transparent")

###### save file
saveRDS(liver2SCT2, file = "liv2_SCT2_Prtc.rds")

## Load the data
liver2SCT2 <- readRDS("liv2_SCT2_Prtc.rds")




##  b. Annotation via deconvolution using a scRNA-seq reference 
##  Methods: DWLS
#Install Giotto Suite
devtools::install_version("Matrix", version = "1.6-5")
install.packages("terra")
library("terra")
devtools::install_github("drieslab/Giotto@suite_dev")
devtools::install_github("drieslab/GiottoVisuals@dev")
devtools::install_github("drieslab/GiottoClass@dev")
devtools::install_github("drieslab/GiottoUtils@dev")
remotes::install_github("drieslab/GiottoData")
pak::pkg_install("drieslab/Giotto")

### Conversion of Seurat V5 to Giotto
# Load required Libraries
library(Giotto)
installGiottoEnvironment()
library(data.table)
library(GiottoData)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(SpatialExperiment)
library(SummarizedExperiment)
library(spacexr)

#### Load the scRNAseq Reference
combined <- readRDS("RefLivcombined.rds")
ref = combined

# Downsample 200 cells per identity group (cell type) (This step is to save time)
seurat_obj = ref
Idents(seurat_obj) <- "celltype"
seurat_downsampled <- subset(seurat_obj, downsample = 200)

# Convert Seurat object to Giotto
ref = seurat_downsampled
head(ref@assays)
ref1 <- seuratToGiottoV5(sobject = ref, spatial_assay = "SCT")
table(ref1@feat_metadata$cell$rna$celltype)
giotto_SC = ref1
giotto_SC <- normalizeGiotto(giotto_SC)

## Calculate the cell type markers
showGiottoCellMetadata(giotto_SC)
markers_scran <- findMarkers_one_vs_all(gobject = giotto_SC, 
                                        method = "scran",
                                        expression_values = "normalized",
                                        cluster_column = "celltype", 
                                        min_feats = 3)
top_markers <- markers_scran[, head(.SD, 10), by = "cluster"]$feats

## Create the signature matrix
sign_matrix <- makeSignMatrixDWLSfromMatrix(
  matrix = getExpression(giotto_SC,
                         values = "normalized",
                         output = "matrix"),
  cell_type = pDataDT(giotto_SC)$celltype,
  sign_gene = top_markers)
head(sign_matrix)

### Load Seurat object for annotation
liver2SCT2 <- readRDS("liv2_SCT2_prtc.rds")
Obj1 <- liver2SCT2

# Convert Seurat object to Giotto
Obj1$cell_type <- Obj1$seurat_clusters
liverG <- seuratToGiottoV5(sobject = Obj1, spatial_assay = "Spatial")

# Print Giotto object
liverG

# load the object
g <- liverG 
showGiottoCellMetadata(g)

## Run the DWLS Deconvolution
g <- runDWLSDeconv(g,
                   sign_matrix = sign_matrix,
                   cluster_column = "cell_type",
                   n_cell = 10, #n_cell = 10,
                   return_gobject = TRUE)

## Visualize
spatDeconvPlot(g, 
               show_image = FALSE,
               radius = 50
)

### Conversion of Giotto to Seurat V5
# Reverting GiottoObj to Seurat
#install.packages("geometry")
library("geometry")
G_labled <- giottoToSeuratV5(g)
G_labled



##  b. Annotation via deconvolution using a scRNA-seq reference 
##  Methods: Robust Cell Type Decomposition (RCTD)
if (!requireNamespace("spacexr", quietly = TRUE)) {
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
}
library(spacexr)

# load in the reference scRNA-seq dataset
combined <- readRDS("RefLivcombined.rds")
ref = combined
Idents(ref) <- "celltype"

counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$celltype)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)
head(reference@cell_types)

# get rid of the number of cells < 25, 
class(reference)
str(reference)
slotNames(reference)
# extract cell type
cell_types <- reference@cell_types

# count number of each cell type
celltype_counts <- table(cell_types)

# find out number > 25  cell type 
valid_celltypes <- names(celltype_counts[celltype_counts > 25])

# get cell ID
selected_cells <- names(cell_types)[cell_types %in% valid_celltypes]

# selected_cells is the cells we want
length(selected_cells)  # see how many
head(selected_cells)    # see context

# subset counts、nUMI、cell_types
counts <- reference@counts[, selected_cells]
nUMI <- reference@nUMI[selected_cells]
cell_types <- reference@cell_types[selected_cells]
table(cell_types)
cell_types <- droplevels(cell_types)

# get new Reference subject
reference_subset <- new("Reference",
                        counts = counts,
                        nUMI = nUMI,
                        cell_types = cell_types)

# Load Seurat object and generate query object
liver2SCT2 <- readRDS("liv2_SCT2_prtc.rds")
Obj1 <- liver2SCT2
liverRCTD = Obj1

counts_hd <- liverRCTD[["SCT"]]$counts
liverRCTD_cells_hd <- colnames(liverRCTD[["SCT"]])
coords <- GetTissueCoordinates(liverRCTD)[liverRCTD_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
RCTD <- create.RCTD(query, reference_subset, max_cores = 4)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
liverRCTD <- AddMetaData(liverRCTD, metadata = RCTD@results$results_df)

# generate figures
myRCTD = liverRCTD
results <- myRCTD@results

# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)
#> Warning in dir.create(resultsdir): 'RCTD_Plots' already exists
# make the plots 
# Plots the confident weights for each cell type as in full_mode (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots the weights for each cell type as in doublet_mode. (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 
# Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
# save file
graph2tif(x = NULL, file='RCTDcelltypes', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")

# generate plots
liverRCTD$celltype <- Idents(liverRCTD)
head(liverRCTD@meta.data)
Idents(liverRCTD) <- "first_type"
DimPlot(liverRCTD, reduction = "umap", label = F)

###### save file
saveRDS(liverRCTD, file = "liverRCTD_Prtc.rds")

## Load the data
liver2SCT2 <- readRDS("liv2_SCT2_Prtc.rds")



## 6. Cell type annotation 
##     c. Manual annotation
## loading packages
library(clusterProfiler)

## load seurat object
liver2SCT2 <- readRDS("liv2_SCT2.rds")
obj=liver2SCT2

## load the biomarkers of cell types
cellmarker = read.gmt("mouse_liver_celltype_markers_MASLD.gmt")

## Plot biomarkers in UMAP or spatial slide
SpatialFeaturePlot(obj, features = c("Cyp2f2", "Hal", "Gls2"))
graph2tif(x = NULL, file='5_x-HepCV_liv2', font = "Arial", cairo = TRUE,   
          width = 8, height = 5, bg = "transparent")

FeaturePlot(obj, features = c("Cyp2f2", "Hal", "Gls2"))
graph2svg(x = NULL, file='5_3-FeatinUMAPClus4&5_liv2', font = "Arial", cairo = TRUE,   
          width = 5, height = 6, bg = "transparent")

VlnPlot(obj, features = c("Cyp2f2", "Hal", "Gls2"))
graph2tif(x = NULL, file='5_x1-HepPV_liv2', font = "Arial", cairo = TRUE,   
                  width = 11, height = 5, bg = "transparent")



## 7. DEG analysis, pathway enrichment analysis and Spatially Variable Genes (SVGs) analysis
## Find difference between hepatocytes (0, 1) , here we compare cluster1 to cluster0
liver2SCT2 <- readRDS("liv2_SCT2_prtc.rds")
liver.mergeIntegr=liver2SCT2
de_markers <- FindMarkers(liver.mergeIntegr, ident.2 = "0: Hepatocytes1", ident.1 = "1: Hepatocytes2")
head(de_markers)
write.csv(de_markers, 'DEhep1vs0.csv')


##### Volcano Plot
### GDF15 vs Vehicle
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(de_markers, plot(avg_log2FC, -log10(p_val_adj), pch=20, main="Volcano plot for hepatocytes 1 vs 0 (blue)", xlim=c(-3,3)))
# my setting: padj<.05&abs(log2FC)>1
with(subset(de_markers, p_val_adj<.05 & abs(avg_log2FC)>0.5), points(avg_log2FC, -log10(p_val_adj), pch=20, col="blue"))
graph2svg(x = NULL, file='10-4-Volcano_Hep1vs0', font = "Arial", cairo = TRUE,   
          width = 3, height = 3.5, bg = "transparent")

# Pathway enrichment analysis: GO enrichment analysis for clusters 1 vs 0
BiocManager::install("GEOquery")
BiocManager::install("pathview")
BiocManager::install("KEGG.db")
BiocManager::install("enrichplot")
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("org.Mm.eg.db")

library(GEOquery)
library(dplyr)
library(pathview)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(GOstats)
library ("org.Mm.eg.db")
library(Ipaper)
library(KEGG.db)

columns(org.Mm.eg.db)
keytypes(org.Mm.eg.db)

###Prapation dataset for enrichment analysis:log2FC >0.8 &padj < 0.1
res_gvv= de_markers
head(res_gvv)

colnames(res_gvv)[2] <- "log2FoldChange"
colnames(res_gvv)[5] <- "padj"
head(res_gvv)
summary(res_gvv, res_gvv$padj < 0.05 & res_gvv$log2FoldChange>0)

### cnetplot for GO log2FC >1&padj < 0.01
###Prapation dataset for enrichment analysis:abs(log2FoldChange) >0 & padj < 0.05
summary(res_gvv, res_gvv$padj < 0.05 & res_gvv$log2FoldChange>0)
res_gvvOrdered1 <- res_gvv[order(res_gvv$padj),]
res_gvvSig1 <- subset(res_gvvOrdered1, padj < 0.05)
res_gvvSig1
res_gvvOrdered2 <- res_gvvSig1[order(res_gvvSig1$log2FoldChange),]
res_gvvSig2 <- subset(res_gvvOrdered2, abs(log2FoldChange) >0)
res_gvvSig2
write.csv(as.data.frame(res_gvvSig2), "DE_HSC_Padj05logFC1_Circle.csv")
res_gvvSig2$SYMBOL=rownames(res_gvvSig2)

### Get dataform
library(ggnewscale)
gene <- res_gvvSig2 [,6]
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
gene.df
ego_all <- enrichGO(gene = gene.df$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "all",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)


res_gvvSig2=as.data.frame(res_gvvSig2)
res_gvvSig2mg=merge(res_gvvSig2, gene.df, by="SYMBOL")
res_gvvSig2mg
geneList = res_gvvSig2mg[,3]
head(geneList)
names(geneList) = as.character(res_gvvSig2mg[,1])
geneList = sort(geneList, decreasing = TRUE)

oragnx <- setReadable(ego_all, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(oragnx, foldChange=geneList, 
         circular = TRUE, 
         colorEdge = TRUE,
         node_label='gene',
         showCategory = 8)

graph2svg(x = NULL, file='10-5-circle_Hepa1vs0', font = "Arial", cairo = TRUE,   
          width = 9, height = 5, bg = "transparent")




# spatially variable genes (SVGs) analysis--pre-annotated anatomical regions
# Spatial heterogeneity
## load packages
library(GSVA)
library(msigdbr)
library(escape)
library(GSEABase)
library(dittoSeq)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)
library(org.Dm.eg.db)
library(clusterProfiler)

## Load the data
liver2SCT2 <- readRDS("liv2_SCT2_prtc.rds")

## Manually select spots interactively
# Launch an interactive plot to select spots
seurat_obj <- liver2SCT2
seurat_obj$celltype = Idents(seurat_obj)

## use condition2a as example 
SpatialDimPlot(seurat_obj, group.by = "celltype", images = "condition2a", alpha = c(0.5, 0.5))
graph2tif(x = NULL, file='11_1-SpatialMetab', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent") 

# Prepare input data and Coordinates
data.input1 = Seurat::GetAssayData(seurat_obj, layer = "data", assay = "SCT") # normalized data matrix
spatial.locs1 = Seurat::GetTissueCoordinates(seurat_obj, scale = NULL, cols = c("imagerow", "imagecol")) 
coords=spatial.locs1

# select region 1_central vein
library(ggplot2)
p <- ggplot(coords, aes(x = x, y = y)) +
  geom_point() +
  #scale_y_reverse() +
  #scale_x_reverse() +
  theme_void()
selected <- CellSelector(p)
str(selected)
head(colnames(seurat_obj))

# Label selected spots
# give names for all spots
seurat_obj$region <- "other"
# give names for selected region
seurat_obj$region[selected] <- "Zone3CV"
table(seurat_obj$region)
SpatialDimPlot(seurat_obj, group.by = "region", images = "condition2a", alpha = c(0.5, 0.5))
graph2tif(x = NULL, file='11_2-SpatialMetab_CV', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent") 

# select region 2_portal vein
library(ggplot2)
p <- ggplot(coords, aes(x = x, y = y)) +
  geom_point() +
  #scale_y_reverse() +
  #scale_x_reverse() +
  theme_void()
selected <- CellSelector(p)
str(selected)
head(colnames(seurat_obj))

# Label selected spots
# give names for selected region
seurat_obj$region[selected] <- "Zone1PV"
table(seurat_obj$region)
SpatialDimPlot(seurat_obj, group.by = "region", images = "condition2a", alpha = c(0.5, 0.5))
graph2tif(x = NULL, file='5_x-SpatialMetab_CV_PV', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent") 


###### save file
saveRDS(seurat_obj, file = "liver2_CV-PV_prtc.rds")

## Load the data
liver2_CV-PV <- readRDS("liver2_CV-PV_prtc.rds")



####### DEG analysis between PV and CV 
##### Find the DE btw ctl and treatment in the same cell types
table(seurat_obj$region)
head(seurat_obj@meta.data)

### add region as Idents
Idents(seurat_obj) <- "region"
de_markers <- FindMarkers(seurat_obj, 
                              ident.1 = "Zone3CV", 
                              ident.2 = "Zone1PV",
                              test.use = "wilcox",   # or choose DESeq2, MAST
                              logfc.threshold = 0.25, 
                              min.pct = 0.1)

##### Volcano Plot
### Zone3CV vs Zone1PV
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(de_markers, plot(avg_log2FC, -log10(p_val_adj), pch=20, main="SVGs_Zone3CV vs Zone1PV (blue)", xlim=c(-3,3)))
# my setting: padj<.05&abs(log2FC)>1
with(subset(de_markers, p_val_adj<.05 & abs(avg_log2FC)>0.5), points(avg_log2FC, -log10(p_val_adj), pch=20, col="blue"))
graph2tif(x = NULL, file='5_x-Volcano_SVGs', font = "Arial", cairo = TRUE,   
          width = 3.8, height = 3.5, bg = "transparent")




# Spatially Variable Genes (SVGs) analysis--annotation-free (***takes 6 hours***)
###### whose are expressed only in specific spatial locations
liver2SCT2 <- readRDS("liv2_SCT2_prtc.rds")
liver2SCT3 <- FindSpatiallyVariableFeatures(liver2SCT2, assay = "SCT", selection.method = "moransi",
                                        features = rownames(liver2SCT2), r.metric = 5, slot = "data")
top.features <- head(rownames(liver2SCT3), 6)
SpatialFeaturePlot(liver2SCT3, features = top.features, ncol = 3, alpha = c(0.1, 1))
graph2tif(x = NULL, file='5_x-SVGs', font = "Arial", cairo = TRUE,   
          width = 8, height = 6, bg = "transparent") 

###### save file
saveRDS(liver2SCT3, file = "liver2SCT3SVGs_prtc.rds")

## Load the data
liver2SCT3 <- readRDS("liver2SCT3SVGs_prtc.rds")




## 8. Integrative analysis across multiple samples or conditions
## Load the data
liver1SCT <- readRDS("liver1SCT_Prtc.rds")
liver2SCT <- readRDS("liver2SCT_Prtc.rds")
liver3SCT <- readRDS("liver3SCT_Prtc.rds")
liver4SCT <- readRDS("liver4SCT_Prtc.rds")

### we provide the merge function.
liver.merge <- merge(liver1SCT, y=c(liver2SCT, liver3SCT, liver4SCT), add.cell.ids = c("condition1a", "condition2a", "condition1b", "condition2b"), project = "STprotocol")
levels(liver.merge)

# this function uses minimum of the median UMI (calculated using the raw UMI counts) of individual objects
# to reverse the individual SCT regression model using minimum of median UMI as the sequencing depth covariate. 
liver.merge=PrepSCTFindMarkers(liver.merge, assay = "SCT", verbose = TRUE)

# integrate data from the two conditions (control1 and control2)
# When aligning two genome sequences together,
# identification of shared/homologous regions can help to interpret differences
# between the sequences as well.
liver.mergeIntegr = liver.merge
liver.mergeIntegr
liver.mergeIntegr[["SCT"]]

# run standard analysis workflow
DefaultAssay(liver.mergeIntegr) <- "SCT"
VariableFeatures(liver.mergeIntegr) <- c(VariableFeatures(liver1SCT), VariableFeatures(liver2SCT),
                                         VariableFeatures(liver3SCT), VariableFeatures(liver4SCT))
liver.mergeIntegr <- RunPCA(liver.mergeIntegr, npcs = 30, verbose = FALSE)

# integration (HarmonyIntegration)
liver.mergeIntegr <- IntegrateLayers(object = liver.mergeIntegr, method = HarmonyIntegration, orig.reduction = "pca",
                                     normalization.method = "SCT", new.reduction = "harmony", verbose = T)
# we can now visualize and cluster the datasets.
liver.mergeIntegr <- FindNeighbors(liver.mergeIntegr, reduction = "harmony", dims = 1:30)
liver.mergeIntegr <- FindClusters(liver.mergeIntegr, verbose = FALSE, resolution = 0.8, cluster.name = "harmony_clusters") # 17 cluster
liver.mergeIntegr <- RunUMAP(liver.mergeIntegr, reduction = "harmony",
                             dims = 1:30, reduction.name = "umap.harmony")
# Visualization
DimPlot(liver.mergeIntegr, reduction = "umap.harmony", label = TRUE, group.by = c("ident", "orig.ident"))
graph2svg(x = NULL, file='17-clusterMergeIntegr_DimPlot', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")
#without labeling
DimPlot(liver.mergeIntegr, reduction = "umap.harmony", label = F, group.by = c("ident", "orig.ident"))
graph2svg(x = NULL, file='17_1-clusterMergeIntegr_DimPlotwoLabel_2samp', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")
graph2svg(x = NULL, file='9_1-clusterMergeIntegr_DimPlotwoLabel_2samp', font = "Arial", cairo = TRUE,   
          width = 7, height = 4, bg = "transparent")

SpatialDimPlot(liver.mergeIntegr, label = TRUE, label.size = 3)
graph2svg(x = NULL, file='17_2-clusterMergeIntegr_SpatialDimPlot', font = "Arial", cairo = TRUE,   
          width = 15, height = 6, bg = "transparent")


### Annotation by GPT-4: the name is same as the merge before integration
# IMPORTANT! Assign your OpenAI API key. See Vignette for details
Sys.setenv(OPENAI_API_KEY = 'xx')

# Load packages
library(GPTCelltype)
library(openai)

# Find Markers
#liver.merge=PrepSCTFindMarkers(liver.merge, assay = "SCT", verbose = TRUE)
all_markers = FindAllMarkers(liver.mergeIntegr, assay = "SCT")
markers = all_markers
# GPT-4 annotation
#res <- gptcelltype(markers, tissuename = 'mouse liver', model = 'gpt-4')
res <- gptcelltype(markers, tissuename = 'liver', model = 'gpt-4')
res

celltype <- c("0: Hepatocytes1", "1: Hepatocytes2", "2: Hepatocytes3", "3: Hepatocytes4", "4: Hepatocytes5","5: Hepatocytes, HSC",
              "6: HSC", "7: Hepatocytes6", "8: Erythrocytes", "9: Kupffer Cells, MΦ", "10: B cells",
              "11: Hepatocytes7", "12: Immune Cells (NK cells, IFN-stimu cells, DC)")

liver.mergeIntegr1=liver.mergeIntegr
names(celltype) <- levels(liver.mergeIntegr1)
liver.mergeIntegr1 <- RenameIdents(liver.mergeIntegr1, celltype)

# Visualize cell type annotation on UMAP
DimPlot(liver.mergeIntegr1)
graph2svg(x = NULL, file='18-Integrat_DimPlot_labeling', font = "Arial", cairo = TRUE,   
          width = 8, height = 5, bg = "transparent")

head(liver.merge1[[]])

### save file
saveRDS(liver.mergeIntegr1, file = "liver_merge_integration1_protc.rds")
## Load the data
liver.mergeIntegr1 <- readRDS("liver_merge_integration1_protc.rds")


## Compare HSC cell features between ctl vs treatment ###########
## Load the data
liver.mergeIntegr1 <- readRDS("liver_merge_integration1_protc.rds")
head(liver.mergeIntegr1[[]])

##### Find the DE btw condition1 and condition2 in HSC
immune.combined=liver.mergeIntegr1
head(immune.combined@meta.data)
tail(immune.combined@meta.data)
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$orig.ident, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
table(immune.combined$celltype.stim)
demarkerHSC <- FindMarkers(immune.combined, ident.1 = "6: HSC_GDF15_1", ident.2 = "6: HSC_PairFed1", verbose = FALSE)
head(demarkerHSC, n = 15)
write.csv(demarkerHSC, 'DE_HSC_ctlTrt.csv')
demarkerHSC <- read.csv("DE_HSC_ctlTrt.csv", header=TRUE)
head(demarkerHSC)
row.names(demarkerHSC)=demarkerHSC$X
demarkerHSC=demarkerHSC[,-1]
de_markers=demarkerHSC


## Volcano Plot with label gene names.
install.packages("ggpubr")
install.packages("ggthemes")
library(ggpubr)
library(ggthemes)

de_markers=as.data.frame(de_markers)
deg.data <-de_markers
head(deg.data)
deg.data$logQ <- -log10(deg.data$p_val_adj)
deg.data$log2FC=deg.data$avg_log2FC
deg.data$FDR=deg.data$p_val_adj
deg.data$ID=rownames(deg.data)

# plot basic volcano
ggscatter(deg.data, x = "log2FC", y = "logQ") + theme_base()
# set a subcolumn and set up and down regulated genes
deg.data$Group = "normal"
deg.data$Group[which( (deg.data$FDR < 0.01) & (deg.data$log2FC > 0.5) )] = "up"
deg.data$Group[which( (deg.data$FDR < 0.01) & (deg.data$log2FC < -0.5) )] = "down"

# see how many genes up and down
table(deg.data$Group)

# plot new volcano
ggscatter(deg.data, x = "log2FC", y = "logQ",
          color = "Group") + theme_base()


# add a new collumn "Label"
deg.data$Label = ""
# get gene logQ>3, |log2FC|>2
res_gvv=deg.data
head(res_gvv)
res_gvvOrdered1 <- res_gvv[order(res_gvv$logQ),]
res_gvvSig1 <- subset(res_gvvOrdered1, logQ > 3)
res_gvvSig1
res_gvvOrdered2 <- res_gvvSig1[order(res_gvvSig1$log2FC),]
res_gvvSig2 <- subset(res_gvvOrdered2, abs(log2FC) >2)
res_gvvSig2
dim(res_gvvSig2)
write.csv(as.data.frame(res_gvvSig2), "Clus4vs3logQ3logFC2_results_wt_Fibr.csv")
res_gvvSig3 <- res_gvvSig2[order(abs(res_gvvSig2$log2FC), decreasing = T),]
res_gvvSig3
res_gvvSig3 <- res_gvvSig3[1:50, ]
write.csv(as.data.frame(res_gvvSig3), "Top50Hcc4Vs3WT_Fibr.csv")

log2FC.genes <- head(res_gvvSig3$ID, 40)
deg.top20.genes <- log2FC.genes
deg.top20.genes <- deg.top20.genes[1:14]

#put them in 'Label'
deg.data$Label[match(deg.top20.genes, deg.data$ID)] <- deg.top20.genes

print (deg.data$Label)
table (deg.data$Label)

# change color and axis labelling to make figure looking better
ggscatter(deg.data, x = "log2FC", y = "logQ",
          color = "Group", 
          palette = c("#00BA38", "#BBBBBB", "#F8766D"),
          size = 2,
          label = deg.data$Label, 
          font.label = 8, 
          repel = T,
          xlab = "log2FC", 
          ylab = "-log10(FDR)") + 
  theme_base() + 
  geom_hline(yintercept = 3, linetype="dashed") +
  geom_vline(xintercept = c(-1,1), linetype="dashed")
library(export)
graph2svg(x = NULL, file='10-5_Volcano_HSC_label_Fibr', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 3.5, bg = "transparent")



# Enrichment analysis--GO enrichment analysis
library(GEOquery)
library(dplyr)
library(pathview)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(GOstats)
library ("org.Mm.eg.db")
library(Ipaper)
library(KEGG.db)

columns(org.Mm.eg.db)
keytypes(org.Mm.eg.db)

###Prapation dataset for enrichment analysis:log2FC >0.8 &padj < 0.1
res_gvv=deg.data
head(res_gvv)
res_gvvOrdered1 <- res_gvv[order(res_gvv$logQ),]
res_gvvSig1 <- subset(res_gvvOrdered1, logQ > 0.05)
res_gvvSig1
res_gvvOrdered2 <- res_gvvSig1[order(res_gvvSig1$log2FC),]
res_gvvSig2 <- subset(res_gvvOrdered2, abs(log2FC) >0.5)
res_gvvSig2
dim(res_gvvSig2)
write.csv(as.data.frame(res_gvvSig2), "Clus4vs3logQ3logFC2_results_wt_Fibr.csv")
res_gvvSig2$SYMBOL=rownames(res_gvvSig2)

### Get dataform
gene <- res_gvvSig2 [,12]
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
gene.df

### GO_all enrichment
ego_all <- enrichGO(gene = gene.df$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01)

head(ego_all)
dim(ego_all)
write.csv(ego_all,'GO_ALL_HSC.csv')
dev.off()
dotplot(ego_all, showCategory=15)+ggtitle("GO enrichment analysis")
dotplot1=dotplot(ego_all, showCategory=15)+ggtitle("GO enrichment analysis")
write_fig(dotplot1, "10-6_GO enrichment analysis.tif", show = FALSE, devices = "tif", width = 8,
          height = 5.5)
dev.off()


## 9. Pseudobulk analysis
# pseudobulking and Perform DE analysis after pseudobulking
### load data (we offer another dataset example which includes 3 replicates per condition)
liver.mergeIntegr1 <- readRDS("liver.mergeIntegr1Meta.rds")
LivPseu = liver.mergeIntegr1

### Rebuilt stim in metadata 
head(LivPseu@meta.data)
LivPseu$stim = LivPseu$group
table(LivPseu$stim)
LivPseu$stim[LivPseu$stim == "control1"] <- "control"
LivPseu$stim[LivPseu$stim == "control2"] <- "control"
LivPseu$stim[LivPseu$stim == "control3"] <- "control"
LivPseu$stim[LivPseu$stim == "treat1"] <- "treat"
LivPseu$stim[LivPseu$stim == "treat2"] <- "treat"
LivPseu$stim[LivPseu$stim == "treat3"] <- "treat"

# pseudobulk the counts based on stim-repli-celltype
pseudo_liv <- AggregateExpression(LivPseu, assays = "SCT", return.seurat = T, 
                                  group.by = c("stim","group", "celltype") #specify a vector, such as c('ident', 'replicate', 'celltype')
)

# each 'cell' is a stim-repli-celltype pseudobulk profile
tail(Cells(pseudo_liv))
pseudo_liv$celltype.stim <- paste(pseudo_liv$celltype, pseudo_liv$stim, sep = "_")
table(pseudo_liv$celltype.stim)

###### save file
saveRDS(pseudo_liv, file = "pseudo_liv.rds")
## Load the data
pseudo_liv <- readRDS("pseudo_liv.rds")



# Next, we perform DE testing on the pseudobulk level for Hepatocytes
head(pseudo_liv@meta.data)
head(pseudo_liv$celltype.stim)
table(pseudo_liv$celltype.stim)
Idents(pseudo_liv) <- "celltype.stim"
bulk.H.de <- FindMarkers(object = pseudo_liv, 
                         ident.1 = "HSCs_treat", 
                         ident.2 = "HSCs_control",
                         test.use = "DESeq2")
head(bulk.H.de, n = 15)
de_markers = bulk.H.de


##### Volcano Plot
### Hepatocytes1_treat vs Hepatocytes1_control
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(de_markers, plot(avg_log2FC, -log10(p_val_adj), pch=20, main="Volcano plot for hepatocyte1 treat vs Ctl (blue)", xlim=c(-3,3)))

# my setting: padj<.05&abs(log2FC)>1
with(subset(de_markers, p_val_adj<.05 & abs(avg_log2FC)>0.5), points(avg_log2FC, -log10(p_val_adj), pch=20, col="blue"))

graph2tif(x = NULL, file='5_x-Volcano_HSCPsudo', font = "Arial", cairo = TRUE,   
          width = 3, height = 3.5, bg = "transparent")



# GO enrichment analysis for tumor clusters 1 vs 0
library(GEOquery)
library(dplyr)
library(pathview)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(GOstats)
library ("org.Mm.eg.db")
library(Ipaper)
library(KEGG.db)

columns(org.Mm.eg.db)
keytypes(org.Mm.eg.db)

###Prapation dataset for enrichment analysis:log2FC >0.8 &padj < 0.1
res_gvv= de_markers
head(res_gvv)
colnames(res_gvv)[2] <- "log2FoldChange"
colnames(res_gvv)[5] <- "padj"
head(res_gvv)
summary(res_gvv, res_gvv$padj < 0.05 & res_gvv$log2FoldChange>0)


### cnetplot for GO log2FC >1&padj < 0.01
###Prapation dataset for enrichment analysis:log2FC >1&padj < 0.01
summary(res_gvv, res_gvv$padj < 0.1 & res_gvv$log2FoldChange>0)
res_gvvOrdered1 <- res_gvv[order(res_gvv$padj),]
res_gvvSig1 <- subset(res_gvvOrdered1, padj < 0.1)
res_gvvSig1
res_gvvOrdered2 <- res_gvvSig1[order(res_gvvSig1$log2FoldChange),]
res_gvvSig2 <- subset(res_gvvOrdered2, abs(log2FoldChange) >0)
res_gvvSig2
write.csv(as.data.frame(res_gvvSig2), "DE_HSC_Padj05logFC1_Circle.csv")
res_gvvSig2$SYMBOL=rownames(res_gvvSig2)

### Get dataform
BiocManager::install("ggnewscale")
library(ggnewscale)

gene <- res_gvvSig2 [,6]
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
gene.df
ego_all <- enrichGO(gene = gene.df$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "all",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

write.csv(ego_all,'GO_ALL_HSC_protc.csv')
dev.off()

dotplot(ego_all, showCategory=15)+ggtitle("GO enrichment analysis")
graph2tif(x = NULL, file='5_x_GO enrichment analysis_HepaPsudo.tif', font = "Arial", cairo = TRUE,   
          width = 8, height = 5.5, bg = "transparent")


## 10. Quantification of cell type composition 
# extract meta data
liver.mergeIntegr1 <- readRDS("liver_merge_integration1_protc.rds")
liver.merge1 = liver.mergeIntegr1
head(liver.merge1@meta.data)
library(data.table)
md = liver.merge1@meta.data %>% as.data.table()

## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md_count = md[, .N, by = c("orig.ident", "seurat_clusters")]
md_count = as.data.frame(md_count)
md_count = md_count[order(md_count$seurat_clusters),]
write.csv(md_count, '10_0_celltype_num.csv')

celltypeAnno = levels(liver.merge1)
write.csv(celltypeAnno, '10_1_celltype_Anno.csv')

grp=c("condition1a", "condition2a", "condition1b", "condition2b")
for (i in grp){
  CellTyp=md_count[md_count$orig.ident==i,]
  totl=sum(CellTyp$N)
  CellTyp$Percent=CellTyp$N/totl
  CellTyp
  nam=paste0('CellTypPerc_protc',i,'.csv')
  write.csv(CellTyp, nam)
  
}

##### Percentage accumulative figure
a=read.csv("CellTypPerc_protc4fig.csv", header=TRUE)
a
ggplot(a, aes(x=CellType, y=Percent, fill=Group)) +
  geom_bar(stat="identity", colour="black") +
  guides(fill=guide_legend(reverse=TRUE)) +
  scale_fill_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept= 50, linetype="dotted", color="blue")

graph2tif(x = NULL, file='5_x-CellTypePercWithNum', font = "Arial", cairo = TRUE,   
          width = 6, height = 5, bg = "transparent")





## 11. Cellular communication
# using ‘CellChat’
install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
install.packages("parallelly")
devtools::install_github("jinworks/CellChat")
ptm = Sys.time()
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)

###### Integration: In order to work with multiple slices (with integration)
### we provide the merge function.
liver.merge <- merge(liver1SCT1, y= liver3SCT1, add.cell.ids = c("condition1a", "condition1b"), project = "STProtocol")
levels(liver.merge)

# this function uses minimum of the median UMI (calculated using the raw UMI counts) of individual objects
# to reverse the individual SCT regression model using minimum of median UMI as the sequencing depth covariate. 
liver.merge=PrepSCTFindMarkers(liver.merge, assay = "SCT", verbose = TRUE)

########## integrate data from the two conditions (control and treatment)
# When aligning two genome sequences together,
# identification of shared/homologous regions can help to interpret differences
# between the sequences as well.
liver.mergeIntegr = liver.merge
liver.mergeIntegr
liver.mergeIntegr[["SCT"]]

# run standard analysis workflow
DefaultAssay(liver.mergeIntegr) <- "SCT"
VariableFeatures(liver.mergeIntegr) <- c(VariableFeatures(liver1SCT1), 
                                         VariableFeatures(liver3SCT1))
liver.mergeIntegr <- RunPCA(liver.mergeIntegr, npcs = 30, verbose = FALSE)

# integration (HarmonyIntegration)
liver.mergeIntegr <- IntegrateLayers(object = liver.mergeIntegr, method = HarmonyIntegration, orig.reduction = "pca",
                                     normalization.method = "SCT", new.reduction = "harmony", verbose = T)
# we can now visualize and cluster the datasets.
liver.mergeIntegr <- FindNeighbors(liver.mergeIntegr, reduction = "harmony", dims = 1:30)
liver.mergeIntegr <- FindClusters(liver.mergeIntegr, verbose = FALSE, resolution = 0.8, cluster.name = "harmony_clusters") # 17 cluster
liver.mergeIntegr <- RunUMAP(liver.mergeIntegr, reduction = "harmony",
                             dims = 1:30, reduction.name = "umap.harmony")

# Visualization
DimPlot(liver.mergeIntegr, reduction = "umap.harmony", label = TRUE, group.by = c("ident", "orig.ident"))
graph2svg(x = NULL, file='17-clusterMergeIntegr_DimPlot', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")
#without labeling
DimPlot(liver.mergeIntegr, reduction = "umap.harmony", label = F, group.by = c("ident", "orig.ident"))
graph2svg(x = NULL, file='17_1-clusterMergeIntegr_DimPlotwoLabel_2samp', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")
graph2svg(x = NULL, file='9_1-clusterMergeIntegr_DimPlotwoLabel_2samp', font = "Arial", cairo = TRUE,   
          width = 7, height = 4, bg = "transparent")

SpatialDimPlot(liver.mergeIntegr, label = TRUE, label.size = 3)
graph2svg(x = NULL, file='17_2-clusterMergeIntegr_SpatialDimPlot', font = "Arial", cairo = TRUE,   
          width = 15, height = 6, bg = "transparent")


### Annotation: the name is same as the merge before integration
# IMPORTANT! Assign your OpenAI API key. See Vignette for details
Sys.setenv(OPENAI_API_KEY = 'xx')
# Load packages
library(GPTCelltype)
library(openai)
# Find Markers
all_markers = FindAllMarkers(liver.mergeIntegr, assay = "SCT")
markers = all_markers

# GPT-4 annotation
#res <- gptcelltype(markers, tissuename = 'mouse liver', model = 'gpt-4')
res <- gptcelltype(markers, tissuename = 'liver', model = 'gpt-4')
res
write.csv(res, 'Celltype_Integ_v10.csv')
celltype <- c("0: Hepatocytes1", "1: Hepatocytes2", "2: Cholangiocytes", "3: Hepatocytes3", "4: Hepatic Stellate Cells","5: Hepatocytes4",
              "6: Hepatocytes5", "7: Erythroblasts", "8: B Cells", "9: Hepatocytes6", "10: Kupffer Cells")

liver.mergeIntegr1=liver.mergeIntegr
names(celltype) <- levels(liver.mergeIntegr1)
liver.mergeIntegr1 <- RenameIdents(liver.mergeIntegr1, celltype)

# Visualize cell type annotation on UMAP
DimPlot(liver.mergeIntegr1)
graph2svg(x = NULL, file='18-Integrat_DimPlot_labeling', font = "Arial", cairo = TRUE,   
          width = 8, height = 5, bg = "transparent")

head(liver.merge1[[]])

### save file
saveRDS(liver.mergeIntegr1, file = "CellCommuniLiv1a3.rds")

## Load the data
liver.mergeIntegr1 <- readRDS("CellCommuniLiv1a3.rds")
#liver.mergeIntegr1 = seurat_object

# Split two datasets from integrative file
object= liver.mergeIntegr1
table(object$orig.ident)
Liv1 <- subset(object, orig.ident == 'condition1a')
table(Liv1$orig.ident)
Liv2 <- subset(object, orig.ident == 'condition1b')
table(Liv1$orig.ident)

### add data
seu1 = Liv1
seu2 = Liv2

# show the image and annotated spots
color.use <- scPalette(nlevels(seu1)); names(color.use) <- levels(seu1)
p1 <- Seurat::SpatialDimPlot(seu1, label = F, label.size = 3, cols = color.use)
color.use <- scPalette(nlevels(seu2)); names(color.use) <- levels(seu2)
p2 <- Seurat::SpatialDimPlot(seu2, label = F, label.size = 3, cols = color.use) + NoLegend()
p1 + p2
graph2svg(x = NULL, file='10-9-ImageAnnoPlots', font = "Arial", cairo = TRUE,   
          width = 8, height = 5, bg = "transparent")

# Prepare input data for CelChat analysis
data.input1 = Seurat::GetAssayData(seu1, slot = "data", assay = "SCT") # normalized data matrix
data.input2 = Seurat::GetAssayData(seu2, slot = "data", assay = "SCT") 
genes.common <- intersect(rownames(data.input1), rownames(data.input2))
colnames(data.input1) <- paste0("A1_", colnames(data.input1))
colnames(data.input2) <- paste0("A2_", colnames(data.input2))
data.input <- cbind(data.input1[genes.common, ], data.input2[genes.common, ])

# define the meta data
# a column named `samples` should be provided for spatial transcriptomics analysis,
# which is useful for analyzing cell-cell communication by aggregating multiple samples/replicates.
# Of note, for comparison analysis across different conditions, 
# users still need to create a CellChat object seperately for each condition.
meta1 = data.frame(labels = Idents(seu1), samples = "A1") # manually create a dataframe consisting of the cell labels
meta2 = data.frame(labels = Idents(seu2), samples = "A2") 
meta <- rbind(meta1, meta2)
rownames(meta) <- colnames(data.input)

# a factor level should be defined for the `meta$labels` and `meta$samples`
meta$labels <- factor(meta$labels, levels = levels(Idents(seu1)))
meta$samples <- factor(meta$samples, levels = c("A1", "A2"))
unique(meta$labels) # check the cell labels
#meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
unique(meta$samples) # check the sample labels

# load spatial transcriptomics information
# Spatial locations of spots from full (NOT high/low) resolution images are required. 
# For 10X Visium, this information is in `tissue_positions.csv`. 
spatial.locs1 = Seurat::GetTissueCoordinates(seu1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs2 = Seurat::GetTissueCoordinates(seu2, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs <- rbind(spatial.locs1, spatial.locs2)
rownames(spatial.locs) <- colnames(data.input)

# Scale factors of spatial coordinates
# For 10X Visium, the conversion factor of converting spatial coordinates
# from Pixels to Micrometers can be computed as the ratio of the theoretical spot size (i.e., 65um)
# over the number of pixels that span the diameter of a theoretical spot size
# in the full-resolution image (i.e., 'spot_diameter_fullres' in pixels in the 'scalefactors_json.json' file).
scalefactors1 = jsonlite::fromJSON(txt = file.path("C:/Yourfolder/Sam19/spatial", 'scalefactors_json.json'))
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2)
scalefactors2 = jsonlite::fromJSON(txt = file.path("C:/Yourfolder/Sam18/spatial", 'scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2)
spatial.factors <- rbind(spatial.factors1, spatial.factors2)
rownames(spatial.factors) <- c("A1", "A2")

#### Create a CellChat object
head(spatial.locs)
spatial.locs=spatial.locs[,-3]
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
cellchat

#### Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat@DB <- CellChatDB.use


#### Preprocessing the expression data for cell-cell communication analysis
# To infer the cell state-specific communications, we identify over-expressed ligands
# or receptors in one cell group and then identify over-expressed
# ligand-receptor interactions if either ligand or receptor is over-expressed.
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#### Part II: Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
#For certain functions, each worker needs access to certain global variables.
#If these are larger than the default limit, you will see this error. 
#To get around this, you can set options(future.globals.maxSize = X),where X is the maximum allowed size in bytes.
#So to set it to 1GB, you would run options(future.globals.maxSize = 1000 * 1024^2). 
#Note that this will increase your RAM usage so set this number mindfully.
options(future.globals.maxSize = 1050 * 1024^2)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)

#Users can filter out the cell-cell communication if there are only few cells in certain cell groups. 
#By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Infer the cell-cell communication at a signaling pathway level
#CellChat computes the communication probability on signaling pathway level by summarizing the communication
#probabilities of all ligands-receptors interactions associated with each signaling pathway.

#NB: The inferred intercellular communication network of each ligand-receptor pair
#and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
cellchat <- computeCommunProbPathway(cellchat)

##### Calculate the aggregated cell-cell communication network
#We can calculate the aggregated cell-cell communication network by
#counting the number of links or summarizing the communication probability.
cellchat <- aggregateNet(cellchat)

#We can also visualize the aggregated cell-cell communication network. 
#For example, showing the number of interactions or the total interaction strength (weights)
#between any two cell groups using circle plot or heatmap plot.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count),
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
graph2svg(x = NULL, file='10-10-CCCNumInte', font = "Arial", cairo = TRUE,   
          width = 8, height = 8, bg = "transparent")


netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight),
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
graph2svg(x = NULL, file='10-11-CCCInteWeig', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")


netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
graph2svg(x = NULL, file='10-12-CCCHeatmap', font = "Arial", cairo = TRUE,   
          width = 4, height = 4, bg = "transparent")


#### Part III: Visualization of cell-cell communication network
#Upon infering the cell-cell communication network, CellChat provides various functionality
#for further data exploration, analysis, and visualization. 
#Here we only showcase the circle plot and the new spatial plot.
#All the signaling pathways showing significant communications can be accessed by following
cellchat@netP$pathways
pathways.show <- c("MIF") 

# Circle plot
par(mfrow=c(1,1), xpd=TRUE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
graph2svg(x = NULL, file='10-13-CCCMacrophagePathway', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")


# Spatial plot
#for Liv1
par(mfrow=c(1,1))
# Setting `vertex.label.cex = 0` to hide the labels on the spatial plot
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A1",
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 0)
graph2svg(x = NULL, file='10-14-CCCMacrophagePathwaySpat_liv1', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")

#for Liv3
par(mfrow=c(1,1))
# Setting `vertex.label.cex = 0` to hide the labels on the spatial plot
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A2",
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 0)
graph2svg(x = NULL, file='10-15-CCCMacrophagePathwaySpat_liv3', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")

#### Compute and visualize the network centrality scores:
# Compute the network centrality scores
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Visualize the computed centrality scores using heatmap, allowing ready
# identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

graph2svg(x = NULL, file='10-16-CCCcentrality', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4, bg = "transparent")

# USER can show this information on the spatial transcriptomics when
# visualizing a signaling network, e.g., bigger circle indicates larger incoming signaling
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A1", layout = "spatial",
                    edge.width.max = 2, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 6, vertex.label.cex = 0)
graph2svg(x = NULL, file='10-17-CCCDSignalSpa_liv1', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4, bg = "transparent")

#### Compute the contribution of each ligand-receptor pair to the overall signaling pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)
graph2svg(x = NULL, file='10-18-CCClig-rec pairContri', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4, bg = "transparent")


## When visualizing gene expression distribution on tissue using spatialFeaturePlot,
## users also need to provide the sample.use as an input.
# Take an input of a few genes
spatialFeaturePlot(cellchat, features = c("Cd74","Cd44"),
                   sample.use = "A1", point.size = 0.8, color.heatmap = "Reds", direction = 1)
graph2svg(x = NULL, file='10-19-CCClig-recGeneExpSpa', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")

#### Part V: Save the CellChat object
saveRDS(cellchat, file = "cellchat_mouse_Liv1a3_PF.rds")




## 12. Metabolic activity analysis
## This part includes following:
# a. Metabolic pathway activity
# b. Metabolic interactions
# c. Metabolites analysis_Flux balance analysis (FBA)



# a. Metabolic pathway activity -- Compare metabolic activities between condition1 vs condition2
# We test activitis using 85 metabolism pathways gene sets
# using irGSEA
library(GSVA)
library(msigdbr)
library(escape)
library(GSEABase)
library(dittoSeq)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)
library(org.Dm.eg.db)
library(clusterProfiler)
library(presto)
library(msigdbr)
library(doMC)
library(dplyr)
library(Seurat)
library(irGSEA)

### load data of Metabolic pathways
irGsealiv <- readRDS("irGsealiv_protc.rds")

# generate stim-celltype metadata
head(irGsealiv@meta.data)
table(irGsealiv$group)
irGsealiv$celltype.stim <- paste(irGsealiv$celltype, irGsealiv$stim, sep = "_")
table(irGsealiv$celltype.stim)

###### save file
saveRDS(irGsealiv, file = "irGsealiv_protc1.rds")

# Load genesets
MetaPW = read.gmt("MetabolicPathway.gmt")
MetaPW
str(MetaPW)
table(MetaPW$term)

# Convert your df (MetaPW) into a list format
MetaPW_list <- split(MetaPW$gene, MetaPW$term)
str(MetaPW_list)
head(MetaPW_list)

# Running irGSEA using own metabolic programes
irGsea_livfinal <- irGSEA.score(object = irGsealiv, assay = "SCT", 
                                slot = "data", custom = T, geneset = MetaPW_list, 
                                method = c("AUCell", "UCell", "singscore", 
                                           "ssgsea", "JASMINE"), #, "viper"
                                aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                                kcdf = 'Gaussian')

###### save file
saveRDS(irGsea_livfinal, file = "irGsealivFinal_protc.rds")

## Load data set
irGsea_livfinal <- readRDS("irGsealivFinal_protc.rds")

## compare metabolic activity in different cell types
sub_obj = irGsea_livfinal
result.dge <- irGSEA.integrate(object = sub_obj, 
                               group.by = "celltype",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell", "singscore", 
                                          "ssgsea", "JASMINE"))
## Plot figures
# heatmap plot
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge,
                                      method = "RRA",
                                      top = 50,
                                      show.geneset = NULL)
irGSEA.heatmap.plot
graph2tif(x = NULL, file='5_x_irGSEA_all_heatPlot_MetGenrGrp1r1', font = "Arial", cairo = TRUE,
          width = 10, height = 10, bg = "transparent")

## bubble.plot
irGSEA.bubble.plot <- irGSEA.bubble(object = result.dge, 
                                    method = "RRA", 
                                    top = 50)
irGSEA.bubble.plot
graph2tif(x = NULL, file='5_x-irGSEA_LivMeta_r1', font = "Arial", cairo = TRUE,   
          width = 8, height = 8, bg = "transparent")

## density heatmap
#Show the expression and distribution of “FAO” in Ucell among clusters.
densityheatmap <- irGSEA.densityheatmap(object = sub_obj,
                                        method = "UCell",
                                        group.by = "celltype",
                                        show.geneset = "Fatty acid oxidation")
densityheatmap
graph2tif(x = NULL, file='5_x-irGSEA_LivMetaDensity_r1', font = "Arial", cairo = TRUE,   
          width = 6, height = 5, bg = "transparent")

## half vlnplot
# Show the expression and distribution of “TCA-Cycle” in Ucell among clusters.
head (sub_obj$celltype)
halfvlnplot <- irGSEA.halfvlnplot(object = sub_obj,
                                  method = "AUCell", group.by = "celltype",
                                  show.geneset = "Citrate cycle (TCA cycle)")
halfvlnplot
graph2tif(x = NULL, file='5_x-irGSEA_LivMetahalfvln_r1', font = "Arial", cairo = TRUE,
          width = 6, height = 4, bg = "transparent")




####### Assess activity of predefined Metabolic pathways (GSEA + cnetplot)
# Load data set
irGsealiv <- readRDS("irGsealiv_protc1.rds")

### Find the DE btw conditions in the same cell types
library(DESeq2)
table(irGsealiv$celltype.stim)

head(irGsealiv@meta.data)
Idents(irGsealiv) <- "celltype.stim"
table (Idents(irGsealiv))
CellDE <- FindMarkers(object = irGsealiv, 
                      ident.1 = "Hepatocytes1_condition1", 
                      ident.2 = "Hepatocytes1_condition2")

deg = CellDE
genelist = deg$avg_log2FC
names(genelist)=rownames(deg)
genelist=sort(genelist, decreasing = T)
head(genelist)


# Load genesets
MetaPW = read.gmt("MetabolicPathway.gmt")
MetaPW
str(MetaPW)
table(MetaPW$term)

# Convert your df (MetaPW) into a list format
MetaPW_list <- split(MetaPW$gene, MetaPW$term)
str(MetaPW_list)
head(MetaPW_list)

## GSEA
egmt = GSEA(genelist, TERM2GENE=MetaPW,
            minGSSize = 1,
            pvalueCutoff = 1,
            verbose = FALSE)
head(egmt)
egmt@result
gsea_results_df = egmt@result
rownames(gsea_results_df)
write.csv(gsea_results_df, file="gsea_DetailKEGG_dfHep3r1.csv")

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(egmt, categorySize="pvalue", foldChange=genelist, showCategory = 10)
graph2tif(x = NULL, file='5_x_GSEA_MetaActivityCnetPlot_3r1', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 5, bg = "transparent")


## 12. Metabolic activity analysis
##     b. Metabolic interactions
# ssGSEA analysis and pathway interaction (correlation) analysis #######

# Load genesets
MetaPW = read.gmt("MetabolicPathway.gmt")
str(MetaPW)
table(MetaPW$term)

# Convert your df (MetaPW) into a list format
MetaPW_list <- split(MetaPW$gene, MetaPW$term)
str(MetaPW_list)
head(MetaPW_list)

## Load data set
irGsealiv <- readRDS("irGsealiv_protc1.rds")
seurat_obj=irGsealiv

## setting analysis
table(seurat_obj$group)
DefaultAssay(seurat_obj) <- "SCT"

# extract data frame
expr <- as.matrix(GetAssayData(seurat_obj, slot = "data", assay = "SCT"))
head(expr)

# define metabolic genesets
head(MetaPW)
head(MetaPW_list)
gene_sets <- MetaPW_list

# ssGSEA
library(GSVA)
library(limma)
scores <- gsva(ssgseaParam(expr, gene_sets), verbose=T)  # results are: gene_set × cell/spot
head(scores)

# reverse and add them to meta.data
scr_ssgsea <- t(scores)
seurat_object <- AddMetaData(seurat_obj, metadata = scr_ssgsea)
head(seurat_object@meta.data)
seurat_object
saveRDS(seurat_object, file = "seurat_object_ssGSEA_Protc.rds")

# Now calculate Spearman correlations (for spatial data) in enrichment scores between different pathways
library(corrplot)
cor_matrix <- cor(scr_ssgsea, method = "spearman")
write.csv(cor_matrix, 'CorrSsgsea_Protc.csv')

# Calculate P value
Pval <- cor.mtest(scr_ssgsea)
write.csv(Pval, 'CorrSsgsea_Pval_Protc.csv')

### plot pathway correlation in ssGSEA analysis
dittoScatterHex(seurat_object, x.var = "Glycolysis / Gluconeogenesis",
                y.var = "Citrate cycle (TCA cycle)", color.var = "celltype",
                colors = c(1:30), max.density = 1.5)
graph2tif(x = NULL, file='5_x_Cor_celltype_r1', font = "Arial", cairo = TRUE,   
          width = 7, height = 6, bg = "transparent")

### compare different stimulations in hepatocytes
library(corrplot)
table(seurat_object$celltype)
table(seurat_object$group)
table(seurat_object$stim)
Idents(seurat_object)
seurat_object_sub <- subset(seurat_object, subset = celltype %in% c("Hepatocytes1","Hepatocytes2",
                                                                    "Hepatocytes3", "Hepatocytes4", "Hepatocytes5",
                                                                    "Hepatocytes6", "Hepatocytes7","Hepatocytes8"))
# run correlation analysis
scr_ssgseaFrame1 <- FetchData(seurat_object_sub, vars = c("Fatty acid oxidation", "Citrate cycle (TCA cycle)", "stim"))
head(scr_ssgseaFrame1)
tail(scr_ssgseaFrame1)
scr_ssgseaFrame1$group <- sub("_.*", "", rownames(scr_ssgseaFrame1))
head(scr_ssgseaFrame1)
library(ggplot2)
ggplot(scr_ssgseaFrame1, aes(x = scr_ssgseaFrame1$`Fatty acid oxidation`, y = scr_ssgseaFrame1$`Citrate cycle (TCA cycle)`, color = stim)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = stim)) +
  labs(title = "Hepatocytes", x = "FAO", 
       y = "TCA cycle") +
  ggpubr::stat_cor(method = "spearman")+
  theme_minimal()

graph2tif(x = NULL, file='5_x-ssGSEA_LivMeta_conditions', font = "Arial", cairo = TRUE,   
          width = 5, height = 4, bg = "transparent")



# cell-cell communication using ‘CellChat’
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)

###### Integration: In order to work with multiple slices (with integration)
### we provide the merge function.
liver1SCT1 <- readRDS("liver1_final.rds")
liver2SCT1 <- readRDS("liver2_final.rds")
liver.merge <- merge(liver1SCT1, y= liver2SCT1, add.cell.ids = c("control1", "treat1"), project = "CCCmeta")
levels(liver.merge)

# this function uses minimum of the median UMI (calculated using the raw UMI counts) of individual objects
# to reverse the individual SCT regression model using minimum of median UMI as the sequencing depth covariate.
liver.merge=PrepSCTFindMarkers(liver.merge, assay = "SCT", verbose = TRUE)

########## integrate data from the two conditions (control and treatment)
# When aligning two genome sequences together,
# identification of shared/homologous regions can help to interpret differences
# between the sequences as well.
liver.mergeIntegr = liver.merge
liver.mergeIntegr
liver.mergeIntegr[["SCT"]]

# run standard analysis workflow
DefaultAssay(liver.mergeIntegr) <- "SCT"
VariableFeatures(liver.mergeIntegr) <- c(VariableFeatures(liver1SCT1),
                                         VariableFeatures(liver2SCT1))
liver.mergeIntegr <- RunPCA(liver.mergeIntegr, npcs = 30, verbose = FALSE)

# integration (HarmonyIntegration)
liver.mergeIntegr <- IntegrateLayers(object = liver.mergeIntegr, method = HarmonyIntegration, orig.reduction = "pca",
                                     normalization.method = "SCT", new.reduction = "harmony", verbose = T)

# we can now visualize and cluster the datasets.
liver.mergeIntegr <- FindNeighbors(liver.mergeIntegr, reduction = "harmony", dims = 1:30)
liver.mergeIntegr <- FindClusters(liver.mergeIntegr, verbose = FALSE, resolution = 0.8, cluster.name = "harmony_clusters") # 17 cluster
liver.mergeIntegr <- RunUMAP(liver.mergeIntegr, reduction = "harmony",
                             dims = 1:30, reduction.name = "umap.harmony")

# Visualization
DimPlot(liver.mergeIntegr, reduction = "umap.harmony", label = TRUE, group.by = c("ident", "orig.ident"))
graph2svg(x = NULL, file='17-clusterMergeIntegr_DimPlotr1', font = "Arial", cairo = TRUE,
          width = 15, height = 5, bg = "transparent")
#without labeling
DimPlot(liver.mergeIntegr, reduction = "umap.harmony", label = F, group.by = c("ident", "orig.ident"))
graph2svg(x = NULL, file='17_1-clusterMergeIntegr_DimPlotwoLabel_2sampr1', font = "Arial", cairo = TRUE,
          width = 15, height = 5, bg = "transparent")
graph2svg(x = NULL, file='9_1-clusterMergeIntegr_DimPlotwoLabel_2sampr1', font = "Arial", cairo = TRUE,
          width = 7, height = 4, bg = "transparent")

SpatialDimPlot(liver.mergeIntegr, label = TRUE, label.size = 3)
graph2svg(x = NULL, file='17_2-clusterMergeIntegr_SpatialDimPlotr1', font = "Arial", cairo = TRUE,
          width = 15, height = 6, bg = "transparent")


### Annotation: the name is same as the merge before integration
# IMPORTANT! Assign your OpenAI API key. See Vignette for details
Sys.setenv(OPENAI_API_KEY = 'xx')
# Load packages
library(GPTCelltype)
library(openai)

# Find Markers
#liver.merge=PrepSCTFindMarkers(liver.merge, assay = "SCT", verbose = TRUE)
all_markers = FindAllMarkers(liver.mergeIntegr, assay = "SCT")
markers = all_markers
# GPT-4 annotation
res <- gptcelltype(markers, tissuename = 'liver', model = 'gpt-4')
res
celltype <- c("Hepatocytes1", "Hepatocytes2", "Hepatocytes3", "Hepatocytes4", "Hepatocytes5","Hepatocytes6",
              "HSCs", "Macrophages", "Hepatocytes7", "Hepatocytes8", "B Cells", "Erythroid cells", "Immune Cells")

liver.mergeIntegr1=liver.mergeIntegr
names(celltype) <- levels(liver.mergeIntegr1)
liver.mergeIntegr1 <- RenameIdents(liver.mergeIntegr1, celltype)

# Assign cell type annotation back to Seurat object
# liver.merge1=liver.merge
# liver.merge1@meta.data$celltype <- as.factor(res[as.character(Idents(liver.merge1))])
# Visualize cell type annotation on UMAP
DimPlot(liver.mergeIntegr1)
graph2svg(x = NULL, file='18-Integrat_DimPlot_labeling', font = "Arial", cairo = TRUE,
          width = 8, height = 5, bg = "transparent")

head(liver.merge1[[]])

### save file
saveRDS(liver.mergeIntegr1, file = "CCCmeta.rds")
## Load the data
liver.mergeIntegr1 <- readRDS("CCCmeta.rds")

# Change name and Split two datasets from integrative file
head(liver.mergeIntegr1@meta.data)
Idents(liver.mergeIntegr1)
levels(liver.mergeIntegr1)
object= liver.mergeIntegr1
table(object$orig.ident)
Idents(object)
object$celltype <- Idents(object)
head(object@meta.data)
object$orig.ident[object$orig.ident == "PairFed1"] <- "control1"
object$orig.ident[object$orig.ident == "GDF15_1"] <- "treat1"
object$stim = object$orig.ident
object$stim[object$stim == "control1"] <- "control"
object$stim[object$stim == "treat1"] <- "treat"
head(object@meta.data)

# split samples
Liv1 <- subset(object, orig.ident == 'control1')
table(Liv1$orig.ident)

Liv2 <- subset(object, orig.ident == 'treat1')
table(Liv2$orig.ident)

### add data
seu1 = Liv1
seu2 = Liv2

# show the image and annotated spots
color.use <- scPalette(nlevels(seu1)); names(color.use) <- levels(seu1)
p1 <- Seurat::SpatialDimPlot(seu1, label = F, label.size = 3, cols = color.use)
color.use <- scPalette(nlevels(seu2)); names(color.use) <- levels(seu2)
p2 <- Seurat::SpatialDimPlot(seu2, label = F, label.size = 3, cols = color.use) + NoLegend()
p1 + p2
graph2svg(x = NULL, file='10-9-ImageAnnoPlotsr1', font = "Arial", cairo = TRUE,   
          width = 8, height = 5, bg = "transparent")


# Prepare input data for CellChat analysis
data.input1 = Seurat::GetAssayData(seu1, layer = "data", assay = "SCT") # normalized data matrix
data.input2 = Seurat::GetAssayData(seu2, layer = "data", assay = "SCT") 
genes.common <- intersect(rownames(data.input1), rownames(data.input2))
colnames(data.input1) <- paste0("A1_", colnames(data.input1))
colnames(data.input2) <- paste0("A2_", colnames(data.input2))
data.input <- cbind(data.input1[genes.common, ], data.input2[genes.common, ])


# define the meta data
# a column named `samples` should be provided for spatial transcriptomics analysis,
# which is useful for analyzing cell-cell communication by aggregating multiple samples/replicates.
# Of note, for comparison analysis across different conditions, 
# users still need to create a CellChat object seperately for each condition.
meta1 = data.frame(labels = Idents(seu1), samples = "A1") # manually create a dataframe consisting of the cell labels
meta2 = data.frame(labels = Idents(seu2), samples = "A2") 
meta <- rbind(meta1, meta2)
rownames(meta) <- colnames(data.input)

# a factor level should be defined for the `meta$labels` and `meta$samples`
meta$labels <- factor(meta$labels, levels = levels(Idents(seu1)))
meta$samples <- factor(meta$samples, levels = c("A1", "A2"))
unique(meta$labels) # check the cell labels
#meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
unique(meta$samples) # check the sample labels

# load spatial transcriptomics information
# Spatial locations of spots from full (NOT high/low) resolution images are required. 
# For 10X Visium, this information is in `tissue_positions.csv`. 
spatial.locs1 = Seurat::GetTissueCoordinates(seu1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs2 = Seurat::GetTissueCoordinates(seu2, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs <- rbind(spatial.locs1, spatial.locs2)
rownames(spatial.locs) <- colnames(data.input)

nrow(spatial.locs)
length(colnames(data.input))
dim(spatial.locs)
dim(data.input)

# Scale factors of spatial coordinates
# For 10X Visium, the conversion factor of converting spatial coordinates
# from Pixels to Micrometers can be computed as the ratio of the theoretical spot size (i.e., 65um)
# over the number of pixels that span the diameter of a theoretical spot size
# in the full-resolution image (i.e., 'spot_diameter_fullres' in pixels in the 'scalefactors_json.json' file).
scalefactors1 = jsonlite::fromJSON(txt = file.path("C:/YourFolder/Sam19/spatial", 'scalefactors_json.json'))
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2)
scalefactors2 = jsonlite::fromJSON(txt = file.path("C:/YourFolder/Sam10/spatial", 'scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2)
spatial.factors <- rbind(spatial.factors1, spatial.factors2)
rownames(spatial.factors) <- c("A1", "A2")

#### Create a CellChat object
head(spatial.locs)
spatial.locs=spatial.locs[,-3]
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
cellchat

#### Set the ligand-receptor interaction database
# Use the CellChatDB with metabolic signaling

# $$$$ only for metabolic signaling
# get pathways related glucose metabolism
CellChatDB <- CellChatDB.mouse
metabolic_pathways <- c("INSULIN", "IGF", "GIPR", "GCG", 
                        "LEP", "ADIPONECTIN", "RESISTIN", "ApoE")
CellChatDB.use <- subsetDB(CellChatDB, search = metabolic_pathways, key = "pathway_name")
cellchat@DB <- CellChatDB.use
unique(CellChatDB.use$interaction$annotation)
unique(CellChatDB.use$interaction$pathway_name)

###$$$ pathways for lipid metabolism
CellChatDB <- CellChatDB.mouse
lipid_pathways <- c("LEP", "ADIPONECTIN", "RESISTIN", "ANGPTL", "ApoE", "ApoA", "ApoB",
                    "LXA4", "27HC", "Cholesterol", "Calcitriol", "Desmosterol",
                    "DHEA", "DHT", "Estradiol", "Progesterone", "Testosterone")
CellChatDB.use <- subsetDB(CellChatDB, search = lipid_pathways, key = "pathway_name")
cellchat@DB <- CellChatDB.use
unique(CellChatDB.use$interaction$annotation)
unique(CellChatDB.use$interaction$pathway_name)


###$$$ pathways for amino acid metabolism
CellChatDB <- CellChatDB.mouse
aa_pathways <- c("GABA-A", "GABA-B", "Glutamate", "Glycine", "SerotoninDopamin", 
                 "Histamine", "IGFBP", "NMU", "NPY", "NTS", "VIP", 
                 "PACAP", "SOMATOSTATIN", "TAFA", "PROK")

CellChatDB.use <- subsetDB(CellChatDB, search = aa_pathways, key = "pathway_name")
cellchat@DB <- CellChatDB.use
unique(CellChatDB.use$interaction$annotation)
unique(CellChatDB.use$interaction$pathway_name)


#### Preprocessing the expression data for cell-cell communication analysis
# To infer the cell state-specific communications, we identify over-expressed ligands
# or receptors in one cell group and then identify over-expressed
# ligand-receptor interactions if either ligand or receptor is over-expressed.

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#### Part II: Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
#For certain functions, each worker needs access to certain global variables.
#If these are larger than the default limit, you will see this error. 
#To get around this, you can set options(future.globals.maxSize = X),where X is the maximum allowed size in bytes.
#So to set it to 1GB, you would run options(future.globals.maxSize = 1000 * 1024^2). 
#Note that this will increase your RAM usage so set this number mindfully.
options(future.globals.maxSize = 1050 * 1024^2)
cellchat@images$coordinates <- as.matrix(cellchat@images$coordinates)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)

#Users can filter out the cell-cell communication if there are only few cells in certain cell groups. 
#By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Infer the cell-cell communication at a signaling pathway level
#CellChat computes the communication probability on signaling pathway level by summarizing the communication
#probabilities of all ligands-receptors interactions associated with each signaling pathway.
#NB: The inferred intercellular communication network of each ligand-receptor pair
#and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
cellchat <- computeCommunProbPathway(cellchat)


##### Calculate the aggregated cell-cell communication network
#We can calculate the aggregated cell-cell communication network by
#counting the number of links or summarizing the communication probability.
cellchat <- aggregateNet(cellchat)

#We can also visualize the aggregated cell-cell communication network. 
#For example, showing the number of interactions or the total interaction strength (weights)
#between any two cell groups using circle plot or heatmap plot.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count),
                 weight.scale = T, label.edge= F, title.name = "Number of interactions (Glu Met)")
graph2tif(x = NULL, file='10-10-CCCNumInte_r1', font = "Arial", cairo = TRUE,   
          width = 8, height = 8, bg = "transparent")


par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight),
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength (Glu Met)")
graph2tif(x = NULL, file='10-11-CCCInteWeig_r1', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")


netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
graph2tif(x = NULL, file='10-12-CCCHeatmap_r1', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")


#### Part III: Visualization of cell-cell communication network
#Upon infering the cell-cell communication network, CellChat provides various functionality
#for further data exploration, analysis, and visualization. 
#Here we only showcase the circle plot and the new spatial plot.

#All the signaling pathways showing significant communications can be accessed by following
cellchat@netP$pathways

pathways.show <- c("IGF") 

# Circle plot
par(mfrow=c(1,1), xpd=TRUE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
graph2tif(x = NULL, file='10-13-CCCIGF_r1', font = "Arial", cairo = TRUE,   
          width = 4, height = 4, bg = "transparent")

# Spatial plot
#for Liv1
par(mfrow=c(1,1))
# Setting `vertex.label.cex = 0` to hide the labels on the spatial plot
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A1",
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 0)
graph2tif(x = NULL, file='10-14-CCCMacrophagePathwaySpat_liv1_r1', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")

#for Liv2
par(mfrow=c(1,1))
# Setting `vertex.label.cex = 0` to hide the labels on the spatial plot
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A2",
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 0)
graph2svg(x = NULL, file='10-15-CCCMacrophagePathwaySpat_liv3', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")

#### Compute and visualize the network centrality scores:
# Compute the network centrality scores
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Visualize the computed centrality scores using heatmap, allowing ready
# identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

graph2tif(x = NULL, file='10-16-CCCcentrality_r1', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4, bg = "transparent")


#### Compute the contribution of each ligand-receptor pair to the overall signaling pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)
graph2tif(x = NULL, file='10-18-CCClig-rec pairContri_r1', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 2.5, bg = "transparent")

## When visualizing gene expression distribution on tissue using spatialFeaturePlot,
## users also need to provide the sample.use as an input.
# Take an input of a few genes
spatialFeaturePlot(cellchat, features = c("Itgav","Itgb3"),
                   sample.use = "A1", point.size = 0.8, color.heatmap = "Reds", direction = 1)
graph2tif(x = NULL, file='10-19-CCClig-recGeneExpSpa_s1_r1', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")

spatialFeaturePlot(cellchat, features = c("Itgav","Itgb3"),
                   sample.use = "A2", point.size = 0.8, color.heatmap = "Reds", direction = 1)
graph2tif(x = NULL, file='10-19-CCClig-recGeneExpSpa_s2_r1', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")

#### Part V: Save the CellChat object
saveRDS(cellchat, file = "cellchat_mouse_Liv1a2_Met.rds")







# Comparism between control and treat using lipid metabolism
## Load the data
liver.mergeIntegr1 <- readRDS("CCCmeta.rds")

# Chnage name and Split two datasets from integrative file
object= liver.mergeIntegr1
table(object$orig.ident)
Idents(object)
object$celltype <- Idents(object)
head(object@meta.data)
object$orig.ident[object$orig.ident == "PairFed1"] <- "control1"
object$orig.ident[object$orig.ident == "GDF15_1"] <- "treat1"
object$stim = object$orig.ident
object$stim[object$stim == "control1"] <- "control"
object$stim[object$stim == "treat1"] <- "treat"

# split samples
Liv1 <- subset(object, orig.ident == 'control1')
table(Liv1$orig.ident)

Liv2 <- subset(object, orig.ident == 'treat1')
table(Liv2$orig.ident)


### add data
seu1 = Liv1
seu2 = Liv2

# Prepare input data for CellChat analysis
data.input1 = Seurat::GetAssayData(seu1, layer = "data", assay = "SCT") # normalized data matrix
data.input2 = Seurat::GetAssayData(seu2, layer = "data", assay = "SCT") 

# define the meta data
# a column named `samples` should be provided for spatial transcriptomics analysis,
# which is useful for analyzing cell-cell communication by aggregating multiple samples/replicates.
# Of note, for comparison analysis across different conditions, 
# users still need to create a CellChat object seperately for each condition.
meta1 = data.frame(labels = Idents(seu1), samples = "control") # manually create a dataframe consisting of the cell labels
meta2 = data.frame(labels = Idents(seu2), samples = "treat") 

# a factor level should be defined for the `meta$labels` and `meta$samples`
meta1$labels <- factor(meta1$labels, levels = levels(Idents(seu1)))
meta1$samples <- factor(meta1$samples, levels = "control")

meta2$labels <- factor(meta2$labels, levels = levels(Idents(seu2)))
meta2$samples <- factor(meta2$samples, levels = "treat")
unique(meta$labels) # check the cell labels
#meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))

unique(meta1$samples) # check the sample labels
unique(meta2$samples)

# load spatial transcriptomics information
# Spatial locations of spots from full (NOT high/low) resolution images are required. 
# For 10X Visium, this information is in `tissue_positions.csv`. 
spatial.locs1 = Seurat::GetTissueCoordinates(seu1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs2 = Seurat::GetTissueCoordinates(seu2, scale = NULL, cols = c("imagerow", "imagecol")) 

nrow(spatial.locs1)
length(colnames(data.input1))
dim(spatial.locs1)
dim(data.input1)

# Scale factors of spatial coordinates
# For 10X Visium, the conversion factor of converting spatial coordinates
# from Pixels to Micrometers can be computed as the ratio of the theoretical spot size (i.e., 65um)
# over the number of pixels that span the diameter of a theoretical spot size
# in the full-resolution image (i.e., 'spot_diameter_fullres' in pixels in the 'scalefactors_json.json' file).
scalefactors1 = jsonlite::fromJSON(txt = file.path("C:/YourFolder/Sam19/spatial", 'scalefactors_json.json'))
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2)

scalefactors2 = jsonlite::fromJSON(txt = file.path("C:/YourFolder/Sam10/spatial", 'scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2)


#### Create a CellChat object
head(spatial.locs1)
head(spatial.locs2)
spatial.locs1=spatial.locs1[,-3]
spatial.locs2=spatial.locs2[,-3]

cellchat1 <- createCellChat(object = data.input1, meta = meta1, group.by = "labels",
                            datatype = "spatial", coordinates = spatial.locs1, spatial.factors = spatial.factors1)
cellchat2 <- createCellChat(object = data.input2, meta = meta2, group.by = "labels",
                            datatype = "spatial", coordinates = spatial.locs2, spatial.factors = spatial.factors2)

cellchat1
cellchat2

#### Set the ligand-receptor interaction database
# Use the CellChatDB with metabolic signaling

###$$$ pathways for lipid metabolism
CellChatDB <- CellChatDB.mouse
lipid_pathways <- c("LEP", "ADIPONECTIN", "RESISTIN", "ANGPTL", "ApoE", "ApoA", "ApoB",
                    "LXA4", "27HC", "Cholesterol", "Calcitriol", "Desmosterol",
                    "DHEA", "DHT", "Estradiol", "Progesterone", "Testosterone")
CellChatDB.use <- subsetDB(CellChatDB, search = lipid_pathways, key = "pathway_name")
cellchat1@DB <- CellChatDB.use #important
cellchat2@DB <- CellChatDB.use #important
unique(CellChatDB.use$interaction$annotation)
unique(CellChatDB.use$interaction$pathway_name)



#### Preprocessing the expression data for cell-cell communication analysis
# To infer the cell state-specific communications, we identify over-expressed ligands
# or receptors in one cell group and then identify over-expressed
# ligand-receptor interactions if either ligand or receptor is over-expressed.

# subset the expression data of signaling genes for saving computation cost
cellchat1 <- subsetData(cellchat1) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)

cellchat2 <- subsetData(cellchat2) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
#execution.time = Sys.time() - ptm
#print(as.numeric(execution.time, units = "secs"))


#### Part II: Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
#For certain functions, each worker needs access to certain global variables.
#If these are larger than the default limit, you will see this error. 
#To get around this, you can set options(future.globals.maxSize = X),where X is the maximum allowed size in bytes.
#So to set it to 1GB, you would run options(future.globals.maxSize = 1000 * 1024^2). 
#Note that this will increase your RAM usage so set this number mindfully.
options(future.globals.maxSize = 1050 * 1024^2)
cellchat1@images$coordinates <- as.matrix(cellchat1@images$coordinates)
cellchat2@images$coordinates <- as.matrix(cellchat2@images$coordinates)

cellchat1 <- computeCommunProb(cellchat1, type = "truncatedMean", trim = 0.1, 
                               distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                               contact.dependent = TRUE, contact.range = 100)

cellchat2 <- computeCommunProb(cellchat2, type = "truncatedMean", trim = 0.1, 
                               distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                               contact.dependent = TRUE, contact.range = 100)

#Users can filter out the cell-cell communication if there are only few cells in certain cell groups. 
#By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)

## Infer the cell-cell communication at a signaling pathway level
#CellChat computes the communication probability on signaling pathway level by summarizing the communication
#probabilities of all ligands-receptors interactions associated with each signaling pathway.

#NB: The inferred intercellular communication network of each ligand-receptor pair
#and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat2 <- computeCommunProbPathway(cellchat2)

##### Calculate the aggregated cell-cell communication network
#We can calculate the aggregated cell-cell communication network by
#counting the number of links or summarizing the communication probability.
cellchat1 <- aggregateNet(cellchat1)
cellchat2 <- aggregateNet(cellchat2)

## netAnalysis_computeCentrality
cellchat1 <- netAnalysis_computeCentrality(cellchat1, slot.name = "netP")
cellchat2 <- netAnalysis_computeCentrality(cellchat2, slot.name = "netP")

#We can also visualize the aggregated cell-cell communication network. 
#For example, showing the number of interactions or the total interaction strength (weights)
#between any two cell groups using circle plot or heatmap plot.
### merge data
object.list <- list(ctl = cellchat1, treat = cellchat2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


# Users can now export the merged CellChat object and the list of the two separate objects for later use
save(object.list, file = "cellchat_object.list_liv_con_trt.RData")
save(cellchat, file = "cellchat_merged_liv_con_trt.RData")


## Part I: Identify altered interactions and cell populations
# Whether the cell-cell communication is enhanced or not?
# Compare the total number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
graph2tif(x = NULL, file='4-1-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,   
          width = 4, height = 3.5, bg = "transparent")


# Compare the number of interactions and interaction strength among different cell populations
# (A) Circle plot showing differential number of interactions or interaction strength
#    among different cell populations across two datasets
# red(or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
graph2tif(x = NULL, file='4-2-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,   
          width = 7, height = 5, bg = "transparent")

# (B) Heatmap showing differential number of interactions or interaction strength
#    among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
graph2tif(x = NULL, file='4-3-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,   
          width = 8, height = 5, bg = "transparent")


# (C) Circle plot showing the number of interactions or interaction strength
#     among different cell populations across multiple datasets
# The above differential network analysis only works for pairwise datasets.
# If there are more datasets for comparison, CellChat can directly show results
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
graph2tif(x = NULL, file='4-4-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,   
          width = 8, height = 5, bg = "transparent")


# (D) Circle plot To simplify the complicated network and gain insights by showing any two/3 cell types
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# between any two cell types using circle plot.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)


###### Compare the major sources and targets in a 2D space
## Identify cell populations with significant changes in
## sending or receiving signals between different datasets
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
graph2tif(x = NULL, file='4-5-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")


### (B) Identify the signaling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "HSCs")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macrophages")
patchwork::wrap_plots(plots = list(gg1,gg2))
graph2tif(x = NULL, file='4-6-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 4, bg = "transparent")



### Identify altered signaling with distinct interaction strength
# (A) Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2

graph2tif(x = NULL, file='4-7-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,   
          width = 6, height = 2.5, bg = "transparent")

### (B) Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)
# combining all the identified signaling pathways from different datasets 
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

graph2tif(x = NULL, file='4-8-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")


#### Part III: Identify the up-gulated and down-regulated signaling ligand-receptor pairs
# Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
graph2tif(x = NULL, file='4-9-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,
          width = 7, height = 5, bg = "transparent") # I did not use this figure


#### the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs 
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in treatment", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in treatment", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
graph2tif(x = NULL, file='4-10-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,
          width = 8, height = 5, bg = "transparent")


# Chord diagram for pathways
pathways.show <- c("ANGPTL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
graph2tif(x = NULL, file='4-11-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,
          width = 14, height = 7, bg = "transparent")  



## Chord figures
levels(object.list[[1]]@idents) 
# [1] "Hepatocytes1"    "Hepatocytes2"    "Hepatocytes3"    "Hepatocytes4"    "Hepatocytes5"   
# [6] "Hepatocytes6"    "HSCs"            "Macrophages"     "Hepatocytes7"    "Hepatocytes8"   
# [11] "B Cells"         "Erythroid cells" "Immune Cells"
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 8, targets.use = 7, lab.cex = 0.5,
                       title.name = paste0("Signaling from macrophage - ", names(object.list)[i]))
}

graph2tif(x = NULL, file='4-12-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,
          width = 14, height = 7, bg = "transparent")  




## 12. Metabolic activity analysis
##     c. Flux balance analysis (FBA)
# FBA - input data preparation
library(Seurat)
library(Matrix)
update.packages("Matrix")
compassFBA <- readRDS("CompassAnal.rds")
DimPlot(compassFBA, reduction = "umap.harmony", group.by = "celltype.stim", label = TRUE)

# export matrix: change 'slot' as 'data'
counts_matrix <- GetAssayData(compassFBA, assay = "SCT", layer = "data")
class(counts_matrix)
counts_matrix <- as(counts_matrix, "dgCMatrix")

# export .mtx format
writeMM(counts_matrix, file = "expression.mtx")

# export row（gene）and column（cells/spot）as tsv files
write.table(rownames(counts_matrix), file = "genes.tsv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(colnames(counts_matrix), file = "sample_names.tsv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# extract and save metadata
metadata_df <- compassFBA@meta.data
head(metadata_df)

write.csv(metadata_df, file = "cell_metadata.csv", quote = TRUE)
