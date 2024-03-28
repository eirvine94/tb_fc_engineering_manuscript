# QC 220302

# Authors: Josh Peters and Edward Irvine

# load libraries
library(tidyverse)
library(Seurat)
library(glue)
library(devtools)
library(MAST)
library(dplyr)
#install_github('immunogenomics/presto')

# load  counts
bctbl <- readRDS("/Volumes/broad_hptmp/jpeters/ei/241_barcodetable.rds")
use <- names(sort(colSums(bctbl), decreasing = TRUE))[3:23]
counts <- bctbl[, use]
counts <- t(counts)
dim(counts)
bc_241 <- counts
colnames(bc_241) <- paste0(colnames(bc_241), "-1")
lmo1 <- CreateSeuratObject(bc_241, assay = "LMO")
lmo1 <- NormalizeData(lmo1, assay = "LMO", normalization.method = "CLR")
lmo1 <- HTODemux(lmo1, assay = "LMO", positive.quantile = 0.99, init = 22, seed = 1, verbose = TRUE, nstarts = 1000)
table(lmo1$hash.ID)
 
bctbl <- readRDS(glue("/Volumes/broad_hptmp/jpeters/ei/242_barcodetable.rds"))
use <- names(sort(colSums(bctbl), decreasing = TRUE))[3:23]
counts <- bctbl[, use]
counts <- t(counts)
bc_242 <- counts
rm(counts)
rm(bctbl)
colnames(bc_242) <- paste0(colnames(bc_242), "-2")
lmo2 <- CreateSeuratObject(bc_242, assay = "LMO")
lmo2 <- NormalizeData(lmo2, assay = "LMO", normalization.method = "CLR")
lmo2 <- HTODemux(lmo2, assay = "LMO", positive.quantile = 0.99, init = 22, seed = 1, verbose = TRUE, nstarts = 1000)
table(lmo2$hash.ID)
 
counts <- Read10X_h5("~/Downloads/filtered_feature_bc_matrix (2).h5")
keep_barcodes <- c(intersect(colnames(bc_241), colnames(counts)), intersect(colnames(bc_242), colnames(counts)))
gex <- counts[, keep_barcodes]
bc1 <- bc_241[, intersect(colnames(bc_241), colnames(counts))]
bc2 <- bc_242[, intersect(colnames(bc_242), colnames(counts))]
bc <- cbind(bc1, bc2)
 
object <- CreateSeuratObject(gex, assay = "RNA", names.delim = "-", names.field = 2)
object$lane <- object$orig.ident
object[["LMO"]] <- CreateAssayObject(counts = bc)
 
o1 <- object[, object$lane == 1]
o1 <- NormalizeData(o1, assay = "LMO", normalization.method = "CLR")
o1 <- HTODemux(o1, assay = "LMO", positive.quantile = 0.99, init = 22, seed = 1, verbose = TRUE, nstarts = 1000)
o1 <- MULTIseqDemux(o1, assay = "LMO", maxiter = 10, verbose = TRUE, quantile = 0.8)
table(o1$hash.ID)
 
o2 <- object[, object$lane == 2]
o2 <- NormalizeData(o2, assay = "LMO", normalization.method = "CLR")
o2 <- HTODemux(o2, assay = "LMO", positive.quantile = 0.99, init = 22, seed = 1, verbose = TRUE, nstarts = 1000)
table(o2$hash.ID)
 
object$sample_id <- "Unassigned"
object$sample_id[Cells(o1)] <- as.character(o1$hash.ID)
object$sample_id[Cells(o2)] <- as.character(o2$hash.ID)
table(object$sample_id)
object$sample_id_v2 <- "Unassigned"
cells1 <- Cells(lmo1)[Cells(lmo1) %in% Cells(object)]
object$sample_id_v2[cells1] <- as.character(lmo1$hash.ID[cells1])
object$sample_id_v2[Cells(lmo2)[Cells(lmo2) %in% Cells(object)]] <- as.character(lmo2$hash.ID[Cells(lmo2) %in% Cells(object)])
table(object$sample_id_v2)
 
fo <- object[, !object$sample_id %in% c("Negative", "Doublet")]
fo[["percent_mito"]] <- PercentageFeatureSet(fo, pattern = "^MT-")
VlnPlot(fo, features = c("nFeature_RNA", "nCount_RNA", "nCount_LMO", "nFeature_LMO", "percent_mito"), ncol = 5)
VlnPlot(fo, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, log = TRUE, split.by = "lane")
fo <- subset(fo, subset = nFeature_RNA > 500 & percent_mito < 25)
fo$cell_barcode <- Cells(fo)
meta <- read_csv("data/ei_metadata_new.csv")
meta <- merge(meta, fo[[]] %>% select(sample_id, cell_barcode), by.x = "barcode_id", by.y = "sample_id")
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell_barcode
fo <- AddMetaData(fo, meta)
saveRDS(fo, "data/filtered_object_v1.rds")
 
fo <- readRDS("data/filtered_object_v1.rds")
fo <- NormalizeData(fo)
fo <- FindVariableFeatures(fo)
fo <- ScaleData(fo)
fo <- RunPCA(fo, npcs = 30)
ElbowPlot(fo, ndims = 30)
fo <- FindNeighbors(fo, dims = 1:20)
fo <- FindClusters(fo, resolution = 0.5, algorithm = 3)
fo <- RunUMAP(fo, dims = 1:20)
fo[[]]

# Add clean sample names metadata
fo$clean_sample_name <- fo$sample_name
fo$clean_sample_name[fo$clean_sample_name == "No mAb"] <- "no Ab"
fo$clean_sample_name[fo$clean_sample_name == "IgG1 24c5"] <- "24c5 IgG1"
fo$clean_sample_name[fo$clean_sample_name == "SEH 24c5"] <- "24c5 SEHFST LS"
fo$clean_sample_name[fo$clean_sample_name == "SAE 24c5"] <- "24c5 SAEAKA"
Idents(object = fo) <- "clean_sample_name"


keep_24c5 <- c("24c5 IgG1", "24c5 SEHFST LS", "no Ab", "Uninfected")
fo_24c5 <- subset(x = fo, idents = keep_24c5)


saveRDS(fo_24c5, "data/filtered_object_v1.rds")


