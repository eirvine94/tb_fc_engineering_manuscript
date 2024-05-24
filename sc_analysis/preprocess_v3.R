#
# Author: Edward B. Irvine
#
# Description:
# Single cell RNA-seq analysis of 24c5 WBA data.
# Uses filtered and annotated object sent over by Josh Peters.
#
# Created: 3/31/22
#
# Modified: 3/31/22
#

##################
##### Housekeeping ----------------------------------------------------------------------------------------------------------------------------------
##################

# load libraries
library(tidyverse)
library(Seurat)
library(glue)
library(devtools)
library(MAST)
library(dplyr)
library(ggrepel)

# Read in annotated object
ao <- readRDS("input_data/annotated_object_v3.rds")



###
### Organize data
###
# Add clean sample names metadata
ao$clean_sample_name <- ao$sample_name
ao$clean_sample_name[ao$clean_sample_name == "No mAb"] <- "no Ab"
ao$clean_sample_name[ao$clean_sample_name == "IgG1 24c5"] <- "24c5 IgG1"
ao$clean_sample_name[ao$clean_sample_name == "SEH 24c5"] <- "24c5 SEHFST LS"
ao$clean_sample_name[ao$clean_sample_name == "SAE 24c5"] <- "24c5 SAEAKA"
ao$clean_sample_name[ao$clean_sample_name == "IgG1 A194"] <- "A194 IgG1"
ao$clean_sample_name[ao$clean_sample_name == "SEH A194"] <- "A194 SEHFST LS"
Idents(object = ao) <- "clean_sample_name"

# Keep only 24c5 conditions
keep_24c5 <- c("24c5 IgG1", "24c5 SEHFST LS", "no Ab", "Uninfected")
ao_24c5 <- subset(x = ao, idents = keep_24c5)

# Remove Mixed cells
Idents(object = ao_24c5) <- "cell_type_l2"
keep_mye <- c("Neutrophil", "CD14 Mono", "CD16 Mono", "pDC", "cDC1")
keep_tnk <- c("B intermediate", "B naive", "ILC", "NK Proliferating", "NK CD56Bright", "NK", "NK Proliferating", "CD4 Naive",
              "Treg", "CD4 TCM", "CD4 TCM Apop", "CD8 TEM", "gdT", "MAIT", "MAIT Apop", "CD8 Naive", "Plasmablast")
keep_descript <- c(keep_mye, keep_tnk)
ao_24c5 <- subset(x = ao_24c5, idents = keep_descript)

sum(table(ao_24c5@meta.data$cell_type_l2))


# Make myeloid and lymphoid data subsets
ao_24c5_mye <- subset(x = ao_24c5, idents = keep_mye)
saveRDS(ao_24c5_mye, "input_data/ao_24c5_mye_v3.rds")
ao_24c5_tnk <- subset(x = ao_24c5, idents = keep_tnk)
saveRDS(ao_24c5_tnk, "input_data/ao_24c5_tnk_v3.rds")









##############################
##### Dimensionality reduction ------------------------------------------------------------------------------------------------------------------------------------------------------------
##############################

###
### 24c5: Myeloid + Lymphoid
###
tiff("plots/combined/combined_umap_cell_type_l2_v3.tiff", units = "in", width = 8, height = 8, res = 300)
DimPlot(ao_24c5, reduction = "umap", group.by = "cell_type_l2", label = FALSE) &
  ggtitle("Cell Clusters") &
  xlab("UMAP 1") &
  ylab("UMAP 2") &
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        #legend.text = element_text(size = 16),
        aspect.ratio = 8/8)
dev.off()

pdf("plots/combined/combined_umap_cell_type_l2_v3.pdf", width = 8, height = 8)
DimPlot(ao_24c5, reduction = "umap", group.by = "cell_type_l2", label = TRUE) &
  ggtitle("Cell Clusters") &
  xlab("UMAP 1") &
  ylab("UMAP 2") &
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        #legend.text = element_text(size = 16),
        aspect.ratio = 8/8)
dev.off()








###################
##### Density plots ------------------------------------------------------------------------------------------------------------------------------------
###################

### Density function
get_density <- function(x, y, ...) {
  # Get density of points in 2 dimensions.
  # @param x A numeric vector.
  # @param y A numeric vector.
  # @param n Create a square n by n grid to compute density.
  # @return The density within each square.
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


###
### 24c5: Myeloid + Lymphoid
###

### 24c5 SEHFST LS
treatment <- "24c5 SEHFST LS"
cells_to_highlight <- ao_24c5$cell_barcode[ao_24c5$clean_sample_name == treatment]
umap_df <- as.data.frame(ao_24c5@reductions$umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$UMAP_1, umap_df_sub$UMAP_2, n = 100)
pdf("plots/24c5/mye/24c5_density_SEH_v3.pdf", width = 4, height = 4)
i <- ggplot() +
  geom_point(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "black", size = 2) +
  geom_point(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "gray95", size = 1) +
  geom_point(data = umap_df_sub, mapping = aes(x = UMAP_1, y = UMAP_2, color = density), size = 0.75) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density", breaks = c(min(umap_df_sub$density), max(umap_df_sub$density)), labels = c("Min", "Max")) +
  labs(x = "UMAP 1", y = "UMAP 2", title = treatment) +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        aspect.ratio = 8/8)
i
dev.off()

### 24c5 IgG1
treatment <- "24c5 IgG1"
cells_to_highlight <- ao_24c5$cell_barcode[ao_24c5$clean_sample_name == treatment]
umap_df <- as.data.frame(ao_24c5@reductions$umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$UMAP_1, umap_df_sub$UMAP_2, n = 100)
pdf("plots/24c5/mye/24c5_density_IgG1_v3.pdf", width = 4, height = 4)
i <- ggplot() +
  geom_point(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "black", size = 2) +
  geom_point(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "gray95", size = 1) +
  geom_point(data = umap_df_sub, mapping = aes(x = UMAP_1, y = UMAP_2, color = density), size = 0.75) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density", breaks = c(min(umap_df_sub$density), max(umap_df_sub$density)), labels = c("Min", "Max")) +
  labs(x = "UMAP 1", y = "UMAP 2", title = treatment) +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        aspect.ratio = 8/8)
i
dev.off()

### no Ab
treatment <- "no Ab"
cells_to_highlight <- ao_24c5$cell_barcode[ao_24c5$clean_sample_name == treatment]
umap_df <- as.data.frame(ao_24c5@reductions$umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$UMAP_1, umap_df_sub$UMAP_2, n = 100)
pdf("plots/24c5/mye/24c5_density_noAb_v3.pdf", width = 4, height = 4)
i <- ggplot() +
  geom_point(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "black", size = 2) +
  geom_point(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "gray95", size = 1) +
  geom_point(data = umap_df_sub, mapping = aes(x = UMAP_1, y = UMAP_2, color = density), size = 0.75) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density", breaks = c(min(umap_df_sub$density), max(umap_df_sub$density)), labels = c("Min", "Max")) +
  labs(x = "UMAP 1", y = "UMAP 2", title = treatment) +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        aspect.ratio = 8/8)
i
dev.off()

### Uninfected
treatment <- "Uninfected"
cells_to_highlight <- ao_24c5$cell_barcode[ao_24c5$clean_sample_name == treatment]
umap_df <- as.data.frame(ao_24c5@reductions$umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$UMAP_1, umap_df_sub$UMAP_2, n = 100)
pdf("plots/24c5/mye/24c5_density_uninfected_v3.pdf", width = 4, height = 4)
i <- ggplot() +
  geom_point(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "black", size = 2) +
  geom_point(data = umap_df, mapping = aes(x = UMAP_1, y = UMAP_2), color = "gray95", size = 1) +
  geom_point(data = umap_df_sub, mapping = aes(x = UMAP_1, y = UMAP_2, color = density), size = 0.75) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density", breaks = c(min(umap_df_sub$density), max(umap_df_sub$density)), labels = c("Min", "Max")) +
  labs(x = "UMAP 1", y = "UMAP 2", title = treatment) +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        aspect.ratio = 8/8)
i
dev.off()




