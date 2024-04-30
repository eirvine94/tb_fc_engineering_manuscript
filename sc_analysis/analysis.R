# description -------------------------------------------------------------

#' Analysis of 3' 10X Genomics GEX on Mtb-infected whole blood
#' Author: Josh Peters

# load libraries ----------------------------------------------------------
library(tidyverse)
library(DropletUtils)
library(glue)
library(Seurat)
library(Azimuth)
library(presto)
library(clustree)
library(rstatix)
library(ggpubr)

# identify cells ----------------------------------------------------------

l1c <- Read10X_h5("data/publication/24h_lane1_raw.h5")
l1cf <- Read10X_h5("data/publication/24h_lane1_filtered.h5")
l1c <- l1c[, colSums(l1c) > 0]

l1ed <- emptyDrops(l1c, lower = 100)
l1_iscell <- l1ed$FDR <= 0.05
l1_cells <- as.data.frame(l1ed) %>% rownames_to_column("barcode") %>% filter(FDR <= 0.05) %>% pull(barcode)
table(l1_cells %in% colnames(l1cf))

l2c <- Read10X_h5("data/publication/24h_lane2_raw.h5")
l2cf <- Read10X_h5("data/publication/24h_lane2_filtered.h5")
l2c <- l2c[, colSums(l2c) > 0]

l2ed <- emptyDrops(l2c, lower = 100)
l2_iscell <- l2ed$FDR <= 0.05
l2_cells <- as.data.frame(l2ed) %>% rownames_to_column("barcode") %>% filter(FDR <= 0.05) %>% pull(barcode)
table(l2_cells %in% colnames(l2cf))
colnames(l2cf) <- gsub("-1", "-2", colnames(l2cf))
colnames(l2c) <- gsub("-1", "-2", colnames(l2c))
l2_cells <- gsub("-1", "-2", l2_cells)

counts <- cbind(l1c[, l1_cells], l2c[, l2_cells])
dim(counts)

object <- CreateSeuratObject(counts, assay = "RNA", names.delim = "-", names.field = 2)
object$lane <- object$orig.ident
table(object$lane)
object$in_10x_filtered <- Cells(object) %in% c(colnames(l1cf), colnames(l2cf))
object$qc_group <- paste0(object$lane, " ", object$in_10x_filtered)

object[["percent_mito"]] <- PercentageFeatureSet(object, pattern = "^MT-")
VlnPlot(object, features = c("percent_mito"), ncol = 1, log = FALSE, group.by = "qc_group")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, log = TRUE, group.by = "qc_group")
VlnPlot(object, features = c("CSF3R"), group.by = "qc_group")

DefaultAssay(object) <- "RNA"
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 30)
object <- FindNeighbors(object, dims = 1:30)
object <- RunUMAP(object, dims = 1:30)
DimPlot(object, group.by = c("lane", "in_10x_filtered"), shuffle = TRUE)

object <- FindClusters(object, resolution = 0.4, algorithm = 3)
DimPlot(object, group.by = c("RNA_snn_res.0.4", "in_10x_filtered"), label = TRUE, shuffle = TRUE) + NoLegend()
markers <- wilcoxauc(object, group_by = "RNA_snn_res.0.4")
markers <- markers %>% filter(padj <= 1e-5) %>% group_by(group) %>% 
  top_n(10, auc) %>% arrange(group, desc(auc))

VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, log = TRUE, group.by = "RNA_snn_res.0.4")
VlnPlot(object, features = c("CSF3R", "NAMPT"), ncol = 2, group.by = "RNA_snn_res.0.4", slot = "data")
FeaturePlot(object, "CSF3R")

object$use_cell <- object$in_10x_filtered
colSums(GetAssayData(object, "data")[c("CSF3R", "NAMPT"), ]) > 4
object$use_cell[colSums(GetAssayData(object, "data")[c("CSF3R", "NAMPT"), ]) > 3] <- TRUE
VlnPlot(object, features = c("CSF3R", "NAMPT"), ncol = 2, group.by = "RNA_snn_res.0.4", slot = "data", split.by = "use_cell")
DimPlot(object, group.by = c("lane", "in_10x_filtered", "use_cell"), shuffle = TRUE)

saveRDS(object, "data/publication/cellqc_object.rds")
use_cells <- object$use_cells
saveRDS(use_cells, "data/publication/use_cells.rds")

object <- readRDS("data/publication/cellqc_object.rds")
use_cells <- readRDS("data/publication/use_cells.rds")

# load data ---------------------------------------------------------------

bctbl <- readRDS("data/publication/241_barcodetable.rds")
use <- names(sort(colSums(bctbl), decreasing = TRUE))[3:23]
counts <- bctbl[, use]
counts <- t(counts)
bc_241 <- counts
colnames(bc_241) <- paste0(colnames(bc_241), "-1")
# lmo1 <- CreateSeuratObject(bc_241, assay = "LMO")
# lmo1 <- NormalizeData(lmo1, assay = "LMO", normalization.method = "CLR")
# lmo1 <- HTODemux(lmo1, assay = "LMO", positive.quantile = 0.99, init = 22, seed = 1, verbose = TRUE, nstarts = 1000)
# lmo1$hash.ID
# recode_bcs <- rownames(bc_241)
# names(recode_bcs) <- 1:21
# d1$Best <- recode(.x = d1$Best, !!!recode_bcs)
# table(d1$Best)
# rownames(d1)
# df <- merge(lmo1[[]], as.data.frame(d1), by = 0)
# table(df$hash.ID == df$Best)
# df$singlet <- df$LMO_maxID == df$Best & df$Doublet == FALSE
# table(df$LMO_maxID == df$Best)

bctbl <- readRDS(glue("data/publication/242_barcodetable.rds"))
use <- names(sort(colSums(bctbl), decreasing = TRUE))[3:23]
counts <- bctbl[, use]
counts <- t(counts)
bc_242 <- counts
colnames(bc_242) <- paste0(colnames(bc_242), "-2")
# lmo2 <- CreateSeuratObject(bc_242, assay = "LMO")
# lmo2 <- NormalizeData(lmo2, assay = "LMO", normalization.method = "CLR")
# lmo2 <- HTODemux(lmo2, assay = "LMO", positive.quantile = 0.99, init = 22, seed = 1, verbose = TRUE, nstarts = 1000)

# load counts -------------------------------------------------------------

# # ca <- Read10X_h5("data/publication/24agg.h5")
# # dim(ca)
# # head(colnames(ca))
# # tail(colnames(ca))
# # 
# # c1f_10x <- Read10X_h5("data/publication/241.h5")
# # dim(c1f_10x)
# c1 <- Read10X_h5("data/publication/241_raw.h5")
# c1 <- c1[, colSums(c1) > 0]
# # dim(c1)
# # c1[1:5, 1:5]
# 
# c1ed <- emptyDrops(c1, lower = 100)
# c1_iscell <- c1ed$FDR <= 0.05
# # sum(c1_iscell, na.rm = TRUE)
# # table(Limited = c1ed$Limited, Significant = c1_iscell)
# c1_cells <- as.data.frame(c1ed) %>% rownames_to_column("barcode") %>% filter(FDR <= 0.05) %>% pull(barcode)
# # c1_cells
# # c1_cells_flipped <- gsub("-1", "-2", c1_cells)
# # table(c1_cells %in% colnames(c1f_10x))
# # table(c1_cells_flipped %in% colnames(ca))
# #table(c1_cells %in% colnames(ca))
# c1f <- c1[, c1_cells]
# # dim(c1f)
# 
# # use_cells <- rownames(e.out)[e.out$FDR <= 0.05 & !is.na(e.out$FDR)]
# # use_cells <- intersect(use_cells, colnames(bc_241))
# # d1 <- hashedDrops(x = bc_241[, use_cells], confident.nmads = 3, pseudo.count = 1, doublet.mixture = FALSE)
# # table(d1$Doublet)
# # table(d1$Confident)
# # table(d1$Best, d1$Confident)
# # c1_metrics <- e.out
# 
# c2 <- Read10X_h5("data/publication/242_raw.h5")
# c2 <- c2[, colSums(c2) > 0]
# # dim(c2)
# # c2[1:5, 1:5]
# 
# c2ed <- emptyDrops(c2, lower = 100)
# c2_iscell <- c2ed$FDR <= 0.05
# # sum(c2_iscell, na.rm = TRUE)
# # table(Limited = c2ed$Limited, Significant = c2_iscell)
# c2_cells <- as.data.frame(c2ed) %>% rownames_to_column("barcode") %>% filter(FDR <= 0.05) %>% pull(barcode)
# #c2_cells <- gsub("-1", "-2", c2_cells)
# #table(c2_cells %in% colnames(ca))
# c2f <- c2[, c2_cells]
# colnames(c2f) <- gsub("-1", "-2", colnames(c2f))
# 
# counts <- cbind(c1f, c2f)
# # dim(counts)

keep_barcodes <- c(intersect(colnames(bc_241), use_cells), intersect(colnames(bc_242), use_cells))
length(keep_barcodes)
gex <- counts[, keep_barcodes]
bc1 <- bc_241[, intersect(colnames(bc_241), use_cells)]
bc2 <- bc_242[, intersect(colnames(bc_242), use_cells)]
bc <- cbind(bc1, bc2)

# create object -----------------------------------------------------------

object <- CreateSeuratObject(gex, assay = "RNA", names.delim = "-", names.field = 2)
object$lane <- object$orig.ident
#table(object$lane)

object[["LMO"]] <- CreateAssayObject(counts = bc)
rm(bc, bc1, bc2, counts, gex, bctbl, c1f, c2f, c1, c2)

# demultiplex samples -----------------------------------------------------

o1 <- object[, object$lane == 1]
o1 <- NormalizeData(o1, assay = "LMO", normalization.method = "CLR")
o1 <- HTODemux(o1, assay = "LMO", positive.quantile = 0.99, init = 22, seed = 1, verbose = TRUE, nstarts = 1000)
#o1 <- MULTIseqDemux(o1, assay = "LMO", maxiter = 10, verbose = TRUE, quantile = 0.8)
# table(o1$hash.ID)

o2 <- object[, object$lane == 2]
o2 <- NormalizeData(o2, assay = "LMO", normalization.method = "CLR")
o2 <- HTODemux(o2, assay = "LMO", positive.quantile = 0.99, init = 22, seed = 1, verbose = TRUE, nstarts = 1000)
# table(o2$hash.ID)

lmo_meta <- rbind(o1[[]] %>% select(LMO_maxID, LMO_secondID, LMO_margin, LMO_classification, LMO_classification.global, hash.ID),
                  o2[[]] %>% select(LMO_maxID, LMO_secondID, LMO_margin, LMO_classification, LMO_classification.global, hash.ID))
object <- AddMetaData(object, lmo_meta)  
# View(object[[]])
table(object$LMO_classification.global)

test_cells <- intersect(colnames(bc_241), use_cells)
d1 <- hashedDrops(x = bc_241[, test_cells], confident.nmads = 3, pseudo.count = 1, doublet.mixture = TRUE)
# table(d1$Doublet)
# table(d1$Confident)
# table(d1$Best, d1$Confident)
recode_bcs <- rownames(bc_241)
names(recode_bcs) <- 1:21
d1$Best <- recode(.x = d1$Best, !!!recode_bcs)
d1 <- as.data.frame(d1) %>% rownames_to_column("cell_barcode")
rownames(d1) <- d1$cell_barcode

test_cells <- intersect(colnames(bc_242), use_cells)
d2 <- hashedDrops(x = bc_242[, test_cells], confident.nmads = 3, pseudo.count = 1, doublet.mixture = TRUE)
# table(d2$Doublet)
# table(d2$Confident)
# table(d2$Best, d2$Confident)
recode_bcs <- rownames(bc_242)
names(recode_bcs) <- 1:21
d2$Best <- recode(.x = d2$Best, !!!recode_bcs)
d2 <- as.data.frame(d2) %>% rownames_to_column("cell_barcode")

d <- rbind(d1, d2)
object <- AddMetaData(object, metadata = d %>% select(-cell_barcode, -Total))
# object[[]]
# View(object[[]])
object$cell_barcode <- Cells(object)
table(object$Confident)

saveRDS(object, "data/publication/raw_object.rds")
rm(o1, o2, bctbl, bc_241, bc_242, d, d1, d2, l1ed, l2ed, lmo_meta)
gc()

# filter object -----------------------------------------------------------

# View(object[[]])
object <- NormalizeData(object, assay = "LMO", normalization.method = "CLR")

object$calls <- "Unassigned"

object$Best[is.na(object$Best)] <- "Unassigned"
object$Confident[is.na(object$Confident)] <- "Unassigned"
object$hash.ID <- as.character(object$hash.ID)

object$calls[object$hash.ID == object$Best] <- as.character(object$hash.ID[object$hash.ID == object$Best])
object$calls[!object$hash.ID %in% c("Negative", "Doublet") & object$hash.ID != object$Best] <- object$hash.ID[!object$hash.ID %in% c("Negative", "Doublet") & object$hash.ID != object$Best]
object$calls[object$Doublet == TRUE & object$LMO_classification.global == "Doublet"] <- "Doublet"
object$calls[object$hash.ID == "Negative"] <- "Negative"
object$calls[object$Confident == TRUE] <- as.character(object$Best[object$Confident == TRUE])
table(object$calls)

saveRDS(object, "data/publication/raw_object.rds")

fo <- object[, !object$calls %in% c("Negative", "Doublet", "Unassigned")]
fo[["percent_mito"]] <- PercentageFeatureSet(fo, pattern = "^MT-")
VlnPlot(fo, features = c("nFeature_RNA", "nCount_RNA", "nCount_LMO", "nFeature_LMO", "percent_mito"), ncol = 5)
VlnPlot(fo, features = c("percent_mito"), ncol = 1, log = FALSE)
VlnPlot(fo, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, log = TRUE, split.by = "lane")
# VlnPlot(fo, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, log = FALSE, split.by = "lane")

# add metadata ------------------------------------------------------------

meta <- read_csv("data/publication/ei_metadata.csv")
meta <- merge(meta, fo[[]] %>% select(calls, cell_barcode), by.x = "barcode_id", by.y = "calls")
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell_barcode
fo <- AddMetaData(fo, meta)
saveRDS(fo, "data/publication/base_object.rds")

# check for batch effects -------------------------------------------------

DefaultAssay(fo) <- "RNA"
fo <- NormalizeData(fo)
fo <- FindVariableFeatures(fo)
fo <- ScaleData(fo)
fo <- RunPCA(fo, npcs = 30)
ElbowPlot(fo, ndims = 30)
fo <- FindNeighbors(fo, dims = 1:30)
fo <- RunUMAP(fo, dims = 1:30)
DimPlot(fo, group.by = c("lane", "donor", "treated", "infected"), shuffle = TRUE)

# compare donors ----------------------------------------------------------

donor_markers <- wilcoxauc(fo, group_by = "donor")
donor_markers <- donor_markers %>% filter(padj <= 0.05 & logFC >= log(1.1))
genes <- read_tsv("data/publication/mart_export.txt")
xy_genes <- genes$`Gene name`[genes$`Chromosome/scaffold name` %in% c("X", "Y")]

# reprocess ---------------------------------------------------------------

fo <- NormalizeData(fo)
fo <- FindVariableFeatures(fo, nfeatures = 3000)
VariableFeatures(fo) <- setdiff(VariableFeatures(fo), xy_genes)
fo <- ScaleData(fo, features = VariableFeatures(fo))
fo <- RunPCA(fo, npcs = 30)
ElbowPlot(fo, ndims = 30)
fo <- JackStraw(fo, dims = 30)
fo <- ScoreJackStraw(fo, dims = 1:30)
JackStrawPlot(fo, dims = 1:30)
fo <- FindNeighbors(fo, dims = 1:30)
fo <- RunUMAP(fo, dims = 1:30)
DimPlot(fo, group.by = c("lane", "donor", "treated", "infected"), 
        shuffle = TRUE, raster = TRUE)

# check clusters ----------------------------------------------------------

fo <- FindClusters(fo, resolution = 0.8, algorithm = 3)
DimPlot(fo, group.by = "RNA_snn_res.0.8", label = TRUE, shuffle = TRUE) + NoLegend()
markers <- wilcoxauc(fo, group_by = "RNA_snn_res.0.8")
markers <- markers %>% filter(padj <= 1e-5) %>% group_by(group) %>% 
  top_n(10, auc) %>% arrange(group, desc(auc))
VlnPlot(fo, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, log = TRUE, pt.size = 0.1)
VlnPlot(fo, features = c("percent_mito"), ncol = 1, log = FALSE, pt.size = 0)
quantile(fo$nFeature_RNA[fo$RNA_snn_res.0.8 == 5], seq(0, 1, 0.1))
quantile(fo$nFeature_RNA[fo$RNA_snn_res.0.8 == 3], seq(0, 1, 0.1))
quantile(fo$nFeature_RNA[fo$RNA_snn_res.0.8 == 15], seq(0, 1, 0.1))
summary(fo$nFeature_RNA[fo$RNA_snn_res.0.8 == 3])
summary(fo$nFeature_RNA[fo$RNA_snn_res.0.8 == 15])

# remove cluster 7
fo <- subset(fo, subset = nFeature_RNA >= 300)
fo <- fo[, fo$RNA_snn_res.0.8 != 6]
fo

# reprocess ---------------------------------------------------------------

fo <- NormalizeData(fo)
fo <- FindVariableFeatures(fo, nfeatures = 3000)
VariableFeatures(fo) <- setdiff(VariableFeatures(fo), xy_genes)
VariableFeatures(fo) <- setdiff(VariableFeatures(fo), rownames(fo)[grep("^MT-", rownames(fo))])
fo <- ScaleData(fo, features = VariableFeatures(fo))
fo <- RunPCA(fo, npcs = 50)
ElbowPlot(fo, ndims = 50)
#fo <- JackStraw(fo, dims = 50)
#fo <- ScoreJackStraw(fo, dims = 1:50)
#JackStrawPlot(fo, dims = 1:50)
fo <- FindNeighbors(fo, dims = 1:30)
fo <- RunUMAP(fo, dims = 1:30)
DimPlot(fo, group.by = c("lane", "donor", "treated", "infected"), 
        shuffle = TRUE, raster = TRUE)

# run azimuth -------------------------------------------------------------

reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")
query <- fo

query <- SCTransform(
  object = query,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

refdata <- lapply(X = c("celltype.l1", "celltype.l2", "celltype.l3"), function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- c("celltype.l1", "celltype.l2", "celltype.l3")
if (TRUE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )
}
query <- TransferData(
  reference = reference$map,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference$map,
  query = query,
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference$map[["refDR"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

query <- NNTransform(
  object = query,
  meta.data = reference$map[[]]
)

query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference$map[["refUMAP"]],
  reduction.key = 'UMAP_'
)

query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)

saveRDS(query, "data/publication/azimuth_object.rds")
query <- readRDS("data/publication/azimuth_object.rds")

fo$predicted_celltype_l1 <- query$predicted.celltype.l1
fo$predicted_celltype_l2 <- query$predicted.celltype.l2
fo$predicted_celltype_l3 <- query$predicted.celltype.l3

saveRDS(fo, "data/publication/processed_object.rds")
rm(query)
gc()

# cluster and define optimal ----------------------------------------------

fo <- readRDS("data/publication/processed_object.rds")
fo <- FindClusters(fo, resolution = seq(0.2, 1.2, 0.1), algorithm = 3)

# check agreement with azimuth predictions
resos <- colnames(fo[[]])[grep("RNA_snn_res.", colnames(fo[[]]))]
ari <- map_dfr(resos, function(x) {
  l1 <- mclust::adjustedRandIndex(fo[["predicted_celltype_l1", drop = TRUE]], fo[[x, drop = TRUE]])
  l2 <- mclust::adjustedRandIndex(fo[["predicted_celltype_l2", drop = TRUE]], fo[[x, drop = TRUE]])
  l3 <- mclust::adjustedRandIndex(fo[["predicted_celltype_l3", drop = TRUE]], fo[[x, drop = TRUE]])
  return(data.frame(l1 = l1, l2 = l2, l3 = l3, res = x))
})
ari$res <- parse_number(gsub("RNA_snn_res.", "", ari$res))
ari <- ari %>% gather(l1:l3, key = "level", value = "ari")
ari <- ari %>% group_by(res) %>% mutate(sum =sum(ari))
ggplot(ari, aes(x = as.factor(res), y = ari, color = level)) +
  geom_point()

# 0.9 resolution

#clustree(fo[[]], prefix = "RNA_snn_res.")
fo <- FindClusters(fo, resolution = seq(0.9), algorithm = 3, n.start = 30, n.iter = 10)
fo$clusters <- fo$RNA_snn_res.0.9
Idents(fo) <- "clusters"
fo <- BuildClusterTree(fo, reorder = TRUE, reorder.numeric = TRUE)
fo$clusters <- Idents(fo)
m <- wilcoxauc(fo)
fm <- m %>% group_by(feature) %>% top_n(1, avgExpr) %>% filter(padj <= 0.05)
topfm <- fm %>% group_by(group) %>% top_n(20, auc)
table(fm$group)
DimPlot(fo, label = TRUE, repel = TRUE, group.by = "clusters")

markers <- wilcoxauc(fo, group_by = "clusters")
markers <- wilcoxauc(fo, group_by = "clusters", groups_use = c(16, 15))
markers <- markers %>% filter(padj <= 1e-5) %>% group_by(group) %>% 
  top_n(10, auc) %>% arrange(group, desc(auc))
VlnPlot(fo, features = c("percent_mito"), group.by = "clusters", ncol = 1, log = FALSE, pt.size = 0)

# call cells --------------------------------------------------------------

toptypesbycluster <- as.data.frame(table(fo$clusters, fo$predicted_celltype_l2)) %>% group_by(Var1) %>% top_n(2, Freq) %>% arrange(Var1, desc(Freq))
topclustersbytype <- as.data.frame(table(fo$predicted_celltype_l2, fo$clusters)) %>% group_by(Var1) %>% top_n(2, Freq) %>% arrange(Var1, desc(Freq))

Idents(fo) <- "clusters"
annotations <- c("Plasmablast", "Neutrophil", "CD14 Mono", "CD16 Mono", "cDC1",
                 "Neutrophil", "ILC", "pDC", "B naive", "B intermediate", 
                 "MAIT", "Mixed", "CD4 TCM", "CD8 TEM", "MAIT Apop",
                 "CD4 TCM Apop", "Treg", "CD8 Naive", "CD4 TCM", "CD4 Naive",
                 "NK Proliferating", "NK CD56Bright", "NK", "CD8 TEM", "gdT")
names(annotations) <- levels(fo)
fo <- RenameIdents(fo, annotations)
fo[["cell_type_l2"]] <- Idents(fo)
fo[["cell_type_l2_num"]] <- paste0(fo$cell_type_l2, " ", fo$clusters)
DimPlot(fo, group.by = "cell_type_l2", label = TRUE, repel = FALSE, ncol = 1) + NoLegend()
DimPlot(fo, group.by = "cell_type_l2_num", label = TRUE, repel = FALSE, ncol = 1, 
        cells = fo$cell_barcode[fo$cell_type_l2 != "Nondescript"]) + NoLegend()

saveRDS(fo, "data/publication/annotated_object.rds")

# plot cell labels --------------------------------------------------------

p <- DimPlot(fo, group.by = "cell_type_l2_num", reduction = "umap", pt.size = 0.25,
             raster = FALSE, label = T, repel = TRUE, shuffle = TRUE, label.size = 2,
             cells = fo$cell_barcode[fo$cell_type_l2 != "Nondescript"]) +
  RemoveAxes() +
  GeneralTheme(14) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  RemoveBackgrounds() +
  colorspace::scale_color_discrete_qualitative("Dark 3") +
  labs(x = "UMAP1", y = "UMAP2", title = "")
plot(p)
SavePlot(plot = p, filename = "celltype_umap", w = 4, h = 4, s = 1, save.data = F)

# proportion stats --------------------------------------------------------

fo$batch <- paste0(fo$donor, fo$lane, " ", fo$sample_name)
sample_count <- fo[[]] %>%
  dplyr::select(batch) %>%
  group_by(batch) %>%
  mutate(batch_frequency = sum(n())) %>%
  distinct()
head(sample_count)

count_table <- fo[[]] %>%
  dplyr::select(batch, cell_type_l2) %>%
  group_by(batch, cell_type_l2, .drop = FALSE) %>%
  dplyr::count(name = "count", .drop = FALSE) %>%
  ungroup()
head(count_table)
count_table <- complete(count_table, batch, cell_type_l2, fill = list(count = 0))
table(count_table$cell_type_l2)
table(count_table$batch)

count_table <- merge(count_table, sample_count, by = "batch", sort = FALSE)
count_table <- count_table %>% mutate(other = batch_frequency - count)
count_table$FA <- count_table$count/count_table$batch_frequency
head(count_table)
count_table <- merge(count_table, fo[[]] %>% select(batch, sample_name) %>% distinct(), by = "batch")
count_table
table(count_table$batch)

comps <- combn(unique(count_table$sample_name), 2)
comps <- lapply(seq_len(ncol(comps)), function(i) comps[,i])
comps <- comps[c(2, 4, 6, 13, 15, 20)]

res <- count_table %>%
  dplyr::group_by(cell_type_l2) %>%
  wilcox_test(FA ~ sample_name, comparisons = comps, p.adjust.method = "none",
              detailed = TRUE) %>%
  adjust_pvalue(p.col = "p", method = "BH") %>%
  add_significance("p.adj")
neutrophil <- res %>% filter(cell_type_l2 == "Neutrophil") %>% adjust_pvalue(p.col = "p", method = "BH")

saveRDS(res, "data/publication/celltypel2_abundance_stats.rds")
write_csv(res, "data/publication/celltypel2_abundance_stats.csv")

# numbered stats ----------------------------------------------------------

fo$batch <- paste0(fo$donor, fo$lane, " ", fo$sample_name)
sample_count <- fo[[]] %>%
  dplyr::select(batch) %>%
  group_by(batch) %>%
  mutate(batch_frequency = sum(n())) %>%
  distinct()
head(sample_count)

count_table <- fo[[]] %>%
  dplyr::select(batch, cell_type_l2_num) %>%
  group_by(batch, cell_type_l2_num, .drop = FALSE) %>%
  dplyr::count(name = "count", .drop = FALSE) %>%
  ungroup()
head(count_table)
count_table <- complete(count_table, batch, cell_type_l2_num, fill = list(count = 0))
table(count_table$cell_type_l2_num)
table(count_table$batch)

count_table <- merge(count_table, sample_count, by = "batch", sort = FALSE)
count_table <- count_table %>% mutate(other = batch_frequency - count)
count_table$FA <- count_table$count/count_table$batch_frequency
count_table <- merge(count_table, fo[[]] %>% select(batch, sample_name) %>% distinct(), by = "batch")
table(count_table$batch)

table(fo$cell_type_l2_num)
ggplot(count_table %>% filter(cell_type_l2_num == "Neutrophil 6"), aes(x = sample_name, y = FA)) +
  geom_boxplot() +
  geom_point()

# plot neutrophils part of myeloid --------------------------------------------------------

fo$batch <- paste0(fo$donor, fo$lane, " ", fo$sample_name)
table(fo$cell_type_l2)
fo$is_myeloid <- ifelse(fo$cell_type_l2 %in% c("Neutrophil", "CD14 Mono", "CD16 Mono", "mregDC", "pDC") , TRUE, FALSE)
sample_count <- fo[[]] %>%
  filter(is_myeloid == TRUE) %>%
  dplyr::select(batch) %>%
  group_by(batch) %>%
  mutate(batch_frequency = sum(n())) %>%
  distinct()
head(sample_count)

count_table <- fo[[]] %>%
  filter(is_myeloid == TRUE) %>%
  dplyr::select(batch, cell_type_l2) %>%
  group_by(batch, cell_type_l2, .drop = FALSE) %>%
  dplyr::count(name = "count", .drop = FALSE) %>%
  ungroup()
head(count_table)
count_table <- complete(count_table, batch, cell_type_l2, fill = list(count = 0))
table(count_table$cell_type_l2)
table(count_table$batch)

count_table <- merge(count_table, sample_count, by = "batch", sort = FALSE)
count_table <- count_table %>% mutate(other = batch_frequency - count)
count_table$FA <- count_table$count/count_table$batch_frequency
head(count_table)

count_table <- merge(count_table, fo[[]] %>% select(batch, sample_name, donor) %>% distinct(), by = "batch")
count_table
table(count_table$batch)
cttp <- count_table %>% filter(cell_type_l2 == "Neutrophil")
unique(cttp$sample_name)
cttp <- cttp %>% filter(sample_name %in% c("IgG1 24c5", "SEH 24c5", "No mAb", "Uninfected"))
cttp$sample_name <- factor(cttp$sample_name, 
                           levels = c("Uninfected", "No mAb", "IgG1 24c5", "SEH 24c5"))
comps <- combn(unique(cttp$sample_name), 2)
comps <- lapply(seq_len(ncol(comps)), function(i) comps[,i])
cttp$log10FA <- log10(cttp$FA)
cttp$log10FA[cttp$FA == 0] <- 0

res <- cttp %>%
  filter(is.finite(FA)) %>%
  wilcox_test(FA ~ sample_name, comparisons = comps, p.adjust.method = "BH",
              detailed = TRUE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position()
p <- ggboxplot(cttp, size = 0.5, width = 0.5, add.params = list(size = 1, alpha = 1),
               x = "sample_name", y = "FA", fill = "sample_name",
               color = "black", add = "jitter", outlier.shape = NA) +
  stat_pvalue_manual(res, label = "{p.adj.signif}", hide.ns = TRUE,
                     tip.length = 0, vjust = 0) +
  labs(x = "Treatment", y = "Fractional abundance", title = "Neutrophils") +
  #scale_y_continuous(labels = scales::comma_format(accuracy = 0.01)) +
                #breaks = c(0.1, 1, 10, 100, 1000),
                #expand = expansion(mult = c(0.05, 0.1))) +
  ggthemes::scale_fill_ptol(name = "") +
  NoLegend() +
  GeneralTheme(14) +
  SpaceAxisTitles() +
  RemoveBackgrounds(outline = TRUE) +
  theme(strip.background = element_rect(color = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major.y = element_line(size = 0.25, color = "gray80", linetype = "dotted"),
        strip.text = element_text(size = 10))
plot(p)
SavePlot(plot = p, filename = "neutrophil_abundance", w = 4.5, h = 4, s = 1.2, save.data = FALSE)

# plot densities ----------------------------------------------------------

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

unique(fo$sample_name)
sample_name <- unique(fo$sample_name)[2]
sub_obj <- fo[, fo$cell_barcode[fo$cell_type_l2 != "Nondescript"]]
cells_to_highlight <- sub_obj$cell_barcode[sub_obj$sample_name == sample_name]
umap_df <- as.data.frame(sub_obj@reductions$umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$UMAP_1, umap_df_sub$UMAP_2, n = 64)
a <- ggplot() +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "black", size = 2) +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "gray95", size = 1) +
  ggrastr::geom_point_rast(data = umap_df_sub, mapping = aes(x = UMAP_1, y = UMAP_2, color = density), size = 1) +
  GeneralTheme(18) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density", breaks = c(min(umap_df_sub$density), max(umap_df_sub$density)), labels = c("Min", "Max")) +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2", title = sample_name)

sample_name <- unique(fo$sample_name)[5]
cells_to_highlight <- sub_obj$cell_barcode[sub_obj$sample_name == sample_name]
umap_df <- as.data.frame(sub_obj@reductions$umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$UMAP_1, umap_df_sub$UMAP_2, n = 64)
b <- ggplot() +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "black", size = 2) +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "gray95", size = 1) +
  ggrastr::geom_point_rast(data = umap_df_sub, mapping = aes(x = UMAP_1, y = UMAP_2, color = density), size = 1) +
  GeneralTheme(18) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density", breaks = c(min(umap_df_sub$density), max(umap_df_sub$density)), labels = c("Min", "Max")) +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2", title = sample_name)

sample_name <- unique(fo$sample_name)[4]
cells_to_highlight <- sub_obj$cell_barcode[sub_obj$sample_name == sample_name]
umap_df <- as.data.frame(sub_obj@reductions$umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$UMAP_1, umap_df_sub$UMAP_2, n = 64)
c <- ggplot() +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "black", size = 2) +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "gray95", size = 1) +
  ggrastr::geom_point_rast(data = umap_df_sub, mapping = aes(x = UMAP_1, y = UMAP_2, color = density), size = 1) +
  GeneralTheme(18) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density", breaks = c(min(umap_df_sub$density), max(umap_df_sub$density)), labels = c("Min", "Max")) +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2", title = sample_name)

grid <- a + b + c
SavePlot(plot = grid, filename = "density_umaps", w = 12, h = 4, s = 1.2, save.data = FALSE)

# summarize de ------------------------------------------------------------

types <- unique(fo$cell_type_l2)
degenes <- map_dfr(.x = types, .f = ~ {
  m1 <- wilcoxauc(fo[, fo$cell_type_l2 == .x], group_by = "sample_name", groups_use = c("IgG1 24c5", "SEH 24c5"))
  m1$comp <- "IgG1 24c5 / SEH 24 c5"
  m2 <- wilcoxauc(fo[, fo$cell_type_l2 == .x], group_by = "sample_name", groups_use = c("No mAb", "SEH 24c5"))
  m2$comp <- "No mAb / SEH 24c5"
  m <- rbind(m1, m2)
  m$cell_type_l2 <- .x
  return(m)
})
dgf <- degenes %>% dplyr::filter(padj <= 0.1 & logFC >= 0 & !feature %in% xy_genes)
dgfs <- dgf %>% group_by(cell_type_l2) %>% summarize(n = n())
write_csv(dgf, "data/publication/antibody_degenes.csv")

a <- ggplot(dgfs, aes(y = fct_reorder(cell_type_l2, n), x = n)) +
  geom_col(color = "black", fill = "gray40") +
  #geom_vline(xintercept = 100, linetype = "dotted", color = "black", size = 1) +
  labs(x = "Number of DE genes", y = "") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 175)) +
  GeneralTheme(14) +
  RemoveBackgrounds(outline = TRUE) +
  SpaceAxisTitles()
plot(a)
SavePlot(a, filename = "degenes_bycelltypel2", root = "", h = 4, w = 4, s = 1, save.data = FALSE)

types <- unique(fo$cell_type_l2_num)
degenes <- map_dfr(.x = types, .f = ~ {
  m1 <- wilcoxauc(fo[, fo$cell_type_l2_num == .x], group_by = "sample_name", groups_use = c("IgG1 24c5", "SEH 24c5"))
  m1$comp <- "IgG1 24c5 / SEH 24 c5"
  m2 <- wilcoxauc(fo[, fo$cell_type_l2_num == .x], group_by = "sample_name", groups_use = c("No mAb", "SEH 24c5"))
  m2$comp <- "No mAb / SEH 24c5"
  m <- rbind(m1, m2)
  m$cell_type_l2_num <- .x
  return(m)
})
degenes <- merge(degenes, fo[[]] %>% select(cell_type_l2, cell_type_l2_num) %>% distinct(), by = "cell_type_l2_num")
head(degenes)
dgf <- degenes %>% filter(pval <= 0.1 & logFC >= 0 & !feature %in% xy_genes)
dgfs_num <- dgf %>% group_by(cell_type_l2_num) %>% summarize(n = n())
write_csv(dgf, "data/publication/antibody_degenes.csv")

unique(fo$sample_name)
table(fo$cell_type_l2)
markers <- wilcoxauc(fo[, fo$cell_type_l2 == "CD14 Mono"], group_by = "sample_name", groups_use = c("IgG1 24c5", "SEH 24c5"))
markers <- markers %>% filter(padj <= 0.1)

# functions ---------------------------------------------------------------

#' RemoveAxes
#'
#' Modified from Seurat::NoAxes()
#'
#' @return
#' @export
RemoveAxes <- function (..., keep.text = FALSE, keep.ticks = FALSE)
{
  blank <- element_blank()
  no_axes_theme <- theme(axis.line.x = blank,
                         axis.line.y = blank,
                         validate = TRUE, ...)
  if (!keep.text) {
    no_axes_theme <- no_axes_theme + theme(
      axis.text.x = blank,
      axis.text.y = blank,
      validate = TRUE, ...)
  }
  
  if (!keep.ticks) {
    no_axes_theme <- no_axes_theme + theme(
      axis.ticks.x = blank,
      axis.ticks.y = blank,
      validate = TRUE, ...)
  }
  
  return(no_axes_theme)
}

#' General plotting theme
#'
#' @param base_size
#' @param ...
#'
#' @return
#' @export
GeneralTheme <- function(base_size, ...) {
  theme_classic(base_size = base_size, ...) +
    ggeasy::easy_all_text_color("black") +
    theme(
      axis.line = element_blank(),
      plot.title = element_text(size =  base_size, color = "black", face = "bold", margin = margin(0,0,4,0)),
      plot.subtitle = element_text(size = base_size - 2, color = "black", margin = margin(0,0,4,0)),
      panel.background = element_rect(fill = "transparent", color = NA, size = 1),
      plot.background = element_rect(fill = "transparent", color = NA, size = 0),
      panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      plot.caption = element_text(hjust = 0, color = "gray40", margin = margin(12)),
      legend.title = element_text(size = base_size, face = "plain"),
      legend.text = element_text(size = base_size - 2),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = "transparent"),
      legend.position = "right",
      legend.justification = "top",
      legend.key.size = unit(1, "line"),
      validate = TRUE
    )
}

#' Remove backgrounds
#'
#' @param outline Keep plot outline
#' @param ...
#'
#' @return
#' @export
RemoveBackgrounds <- function(outline = FALSE, ...)
{
  if (outline) {
    no_bg_theme <- theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1),
                         plot.background = element_rect(fill = "transparent", color = NA, size = 0),
                         legend.background = element_rect(fill = "transparent", color = NA, size = 0),
                         legend.box.background = element_rect(fill = "transparent", color = NA, size = 0),
                         panel.border = element_rect(fill = "transparent", color = NA, size = 0),
                         panel.grid = element_blank(),
                         axis.line = element_blank(),
                         validate = TRUE, ...)
  } else {
    no_bg_theme <- theme(panel.background = element_rect(fill = "transparent", color = NA, size = 0),
                         plot.background = element_rect(fill = "transparent", color = NA, size = 0),
                         legend.background = element_rect(fill = "transparent", color = NA, size = 0),
                         legend.box.background = element_rect(fill = "transparent", color = NA, size = 0),
                         panel.border = element_rect(fill = "transparent", color = NA, size = 0),
                         panel.grid = element_blank(),
                         validate = TRUE, ...)
  }
  
  return(no_bg_theme)
}

#' Save png, pdf and plot object
#'
#' @param plot `ggplot2` object to save
#' @param filename Filename
#' @param root Directories to add to /plots or /data
#' @param h Height
#' @param w Width
#' @param s Scale
#'
#' @return NULL
#' @export
SavePlot <- function(
  plot,
  filename,
  root = "publication",
  h = 6,
  w = 6,
  s = 1,
  save.data = TRUE
) {
  ggplot2::ggsave(plot = plot, filename = glue::glue("plots/{root}/{filename}.png"),
                  scale = s, width = w, height = h, units = "in", dpi = 300)
  ggplot2::ggsave(plot = plot, filename = glue::glue("plots/{root}/{filename}.pdf"),
                  scale = s, width = w, height = h, units = "in", dpi = 300)
  if (save.data) {
    saveRDS(plot, file = glue::glue("data/{root}/plots/{filename}.rds"))
  }
  usethis::ui_done("Saved")
}


#' Space axis titles away from plot area
#'
#' @param scale Increase spacing
#' @param ...
#'
#' @return
#' @export
SpaceAxisTitles <- function(scale = 1, ...) {
  theme (
    axis.title.x = element_text(face = "plain", margin = margin(12*scale, 0, 0, 0)),
    axis.title.y = element_text(face = "plain", margin = margin(0, 12*scale, 0, 0)),
    validate = TRUE
  )
}

