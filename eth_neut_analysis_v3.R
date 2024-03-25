#####
#
# Author: Edward B. Irvine
#
# Description:
# This script analyses the 24c5 whole-blood single-cell RNA-sequencing dataset. Initially, it filters the dataset to focus on neutrophils.
# In the data preparation section, the script identifies differentially expressed genes between sample conditions and extracts significant genes for further analysis, organizing them into subsets for comparison. 
# It then arranges these data subsets for visualization, creating a quad plot to illustrate the differential expression between conditions, specifically comparing log fold changes between different antibody treatments.
#
# Created: 3/31/22
#
# Modified: 3/25/24
#
#####



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
library(ggrastr)
library(presto)
library(msigdbr)
library(fgsea)
library(ggrepel)

# Read in filtered aoect
ao_24c5_mye <- readRDS("input_data/ao_24c5_mye_v3.rds")

# Make neutrophil only object
ao_24c5_neut <- ao_24c5_mye
Idents(object = ao_24c5_neut) <- "cell_type_l2"
ao_24c5_neut <- subset(x = ao_24c5_neut, idents = "Neutrophil")
Idents(object = ao_24c5_neut) <- "clean_sample_name"

# Get number of cells per cluster and per sample of origin
table(ao_24c5_neut@meta.data$cell_type_l2)










###############
##### Prep data ----------------------------------------------------------------------------------------------------------------------------------
###############

# Create ranked gene lists
Idents(object = ao_24c5_neut) <- "clean_sample_name"
gsea_24c5_SEH_IgG1_input <- FindMarkers(object = ao_24c5_neut, ident.1 = "24c5 SEHFST LS", ident.2 = "24c5 IgG1", min.pct = 0, logfc.threshold = 0)
#write.csv(gsea_24c5_SEH_IgG1_input, "gsea/24c5/fgsea/input/gsea_24c5_SEH_IgG1_input_v3.csv")

Idents(object = ao_24c5_neut) <- "clean_sample_name"
gsea_24c5_SEH_noAb_input <- FindMarkers(object = ao_24c5_neut, ident.1 = "24c5 SEHFST LS", ident.2 = "no Ab", min.pct = 0, logfc.threshold = 0)
#write.csv(gsea_24c5_SEH_noAb_input, "gsea/24c5/fgsea/input/gsea_24c5_SEH_noAb_input_v3.csv")

# Pull out genes of interest (p < 0.1)
markers_24c5_neut_SEHvNoAb_quad <- FindMarkers(object = ao_24c5_neut, ident.1 = "24c5 SEHFST LS", ident.2 = "no Ab", min.pct = 0.10, logfc.threshold = 0.25)
markers_24c5_neut_SEHvNoAb_quad <- markers_24c5_neut_SEHvNoAb_quad[markers_24c5_neut_SEHvNoAb_quad$p_val < 0.1, ]
markers_24c5_neut_SEHvIgG1_quad <- FindMarkers(object = ao_24c5_neut, ident.1 = "24c5 SEHFST LS", ident.2 = "24c5 IgG1", min.pct = 0.10, logfc.threshold = 0.25)
markers_24c5_neut_SEHvIgG1_quad <- markers_24c5_neut_SEHvIgG1_quad[markers_24c5_neut_SEHvIgG1_quad$p_val < 0.1, ]

# Make data subsets
noAb_subset <- rownames(markers_24c5_neut_SEHvNoAb_quad)
IgG1_subset <- rownames(markers_24c5_neut_SEHvIgG1_quad)
keep <- c(unique(c(noAb_subset, IgG1_subset)))

SEH_noAb_subset <- gsea_24c5_SEH_noAb_input[rownames(gsea_24c5_SEH_noAb_input) %in% keep, ]
SEH_noAb_subset$gene_name <- rownames(SEH_noAb_subset)
colnames(SEH_noAb_subset)[1:5] <- c("p_val_noAb", "avg_log2FC_noAb", "pct.1_noAb", "pct.2_noAb", "p_val_adj_noAb")
SEH_noAb_subset <- SEH_noAb_subset %>%
  arrange(gene_name)

SEH_IgG1_subset <- gsea_24c5_SEH_IgG1_input[rownames(gsea_24c5_SEH_IgG1_input) %in% keep, ]
SEH_IgG1_subset$gene_name <- rownames(SEH_IgG1_subset)
colnames(SEH_IgG1_subset)[1:5] <- c("p_val_IgG1", "avg_log2FC_IgG1", "pct.1_IgG1", "pct.2_IgG1", "p_val_adj_IgG1")
SEH_IgG1_subset <- SEH_IgG1_subset %>%
  arrange(gene_name)

SEH_noAb_IgG1_subset_df <- data.frame(cbind(SEH_noAb_subset, SEH_IgG1_subset))










###############
##### Quad plot -----------------------------------------------------------------------------------------------------------------------------------
###############

# Organize data for plot
SEH_noAb_IgG1_subset_df$sign <- SEH_noAb_IgG1_subset_df$avg_log2FC_noAb*SEH_noAb_IgG1_subset_df$avg_log2FC_IgG1
SEH_noAb_IgG1_subset_df$label <- SEH_noAb_IgG1_subset_df$gene_name
SEH_noAb_IgG1_subset_df$label[SEH_noAb_IgG1_subset_df$sign < 0] <- "" 
SEH_noAb_IgG1_subset_df <- SEH_noAb_IgG1_subset_df[SEH_noAb_IgG1_subset_df$pct.1_noAb > 0.25, ]
#write.csv(SEH_noAb_IgG1_subset_df, "220407_SEH_noAb_IgG1_subset_df_v3.csv")

# Generate plot
font_size <- 15
png("plots/SEH_quad_plot_v3.png", units = "in", width = 5.6, height = 5.6, res = 300)
quad_plot <- SEH_noAb_IgG1_subset_df %>%
  ggplot(aes(x = avg_log2FC_noAb, y = avg_log2FC_IgG1, label = label)) +
  annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "#941100", alpha = 0.65) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0 , fill= "#005493", alpha = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.25, color = "black", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.25, color = "black", alpha = 0.85) +
  geom_point(size = 1.75, alpha = 0.65) +
  geom_label_repel(size = 4, max.overlaps = 32) +
  scale_x_continuous(limits = c(-0.85, 0.85)) +
  scale_y_continuous(limits = c(-2.9, 2.9)) +
  xlab("Log2 fold change (SEHFST LS / no Ab)") +
  ylab("Log2 fold change (SEHFST LS / IgG1)") +
  theme_linedraw() +
  theme(axis.text = element_text(size = font_size), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        axis.title = element_text(size = font_size),
        legend.text=element_text(size = font_size),
        aspect.ratio = 5/5)
quad_plot
dev.off()


