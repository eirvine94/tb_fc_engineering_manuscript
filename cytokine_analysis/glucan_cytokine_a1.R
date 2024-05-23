#
# Author: Edward B. Irvine
#
# Description:
# Luminex cytokine analysis from supernatants of 24c5 WBA performed on 8/1/18 and 9/19/18.
#
# Created: 17 December 2021
#
# Modified: 23 May 2024
#

###################
###### Housekeeping -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###################

# Load required libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(factoextra)
library(ggforce)
library(ggrepel)
library(RColorBrewer)
library(colorspace)
library(gplots)

# Read in data
dat <- read.csv("211215_24c5_WBA_cytokines.csv")










#######################
###### Pre-process data -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######################

colnames(dat)[4:ncol(dat)] <- c("IL-1β", "IL-10", "IL-6", "GM-CSF", "IL-5", "IFNγ", "TNFα", "IL-2", "IL-4", "IL-8")

# Drop unneeded samples
dat <- dat %>% 
  filter(Outcome != "none")

# Make plate 1 (8/1/18) subset
dat_plate1 <- dat %>%
  filter(Plate == 1)

# Make plate 2 (9/19/18) subset
dat_plate2 <- dat %>%
  filter(Plate == 2)


# Take replicate average
dat_plate1$Outcome <- as.factor(dat_plate1$Outcome)
dat_plate1 <- aggregate(. ~ Sample, dat_plate1, mean)
dat_plate1$Outcome[dat_plate1$Outcome == 1] <- "Non-restrictive"
dat_plate1$Outcome[dat_plate1$Outcome == 2] <- "Restrictive"
rownames(dat_plate1) <- dat_plate1$Sample
dat_plate1["Sample"] <- list(NULL)

dat_plate2$Outcome <- as.factor(dat_plate2$Outcome)
dat_plate2 <- aggregate(. ~ Sample, dat_plate2, mean)
dat_plate2$Outcome[dat_plate2$Outcome == 1] <- "Non-restrictive"
dat_plate2$Outcome[dat_plate2$Outcome == 2] <- "Restrictive"
rownames(dat_plate2) <- dat_plate2$Sample
dat_plate2["Sample"] <- list(NULL)

# Create color palatte
restrictive.col <- "#007F00"
nonrestrictive.col <- "grey70"

dat_plate1$Color <- rep(nonrestrictive.col, nrow(dat_plate1))
dat_plate1$Color[dat_plate1$Outcome == "Restrictive"] <- restrictive.col
dat_plate1$Color <- as.factor(dat_plate1$Color)
#write.csv(dat_plate1, "211215 24c5 WBA Cytokine Luminex Summarized P1 Data.csv")

dat_plate2$Color <- rep(nonrestrictive.col, nrow(dat_plate2))
dat_plate2$Color[dat_plate2$Outcome == "Restrictive"] <- restrictive.col
dat_plate2$Color <- as.factor(dat_plate2$Color)
#write.csv(dat_plate2, "211215 24c5 WBA Cytokine Luminex Summarized P2 Data.csv")










##########
###### PCA -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########

###
### Plate 1
###

# Run PCA
pca_plate1 <- prcomp(as.matrix(dat_plate1[3:12]), center = TRUE, scale = TRUE)

# PCA score plot
tiff("pca_plate1_score.tiff", units = "in", width = 6, height = 6, res = 300)
fviz_pca_ind(pca_plate1,
             habillage = dat_plate1$Outcome,
             repel = TRUE,
             geom = c("point"),
             pointsize = 4,
             alpha = 0.8,
             invisible = "quali",
             pch = 19) + 
  xlim(-5.35, 5.35) +
  ylim(-5.35, 5.35) +
  theme_linedraw() +
  scale_color_manual(values = c(nonrestrictive.col, restrictive.col)) +
  geom_text_repel(aes(label = rownames(dat_plate1), color = dat_plate1$Outcome), size = 4.5) +
  theme(plot.title = element_blank(),
        axis.text = element_text(size = 22), 
        axis.title = element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        legend.title = element_blank(),
        legend.text = element_text(size=22),
        legend.position = "top",
        aspect.ratio = 6/6) 
dev.off()

# PCA loading plot
tiff("pca_plate1_loading.tiff", units = "in", width = 6, height = 6, res = 300)
fviz_pca_var(pca_plate1,
             geom = c("arrow"),
             repel = TRUE,
             col.var = "grey70") +
  theme_linedraw() +
  xlim(-1, 1) +
  ylim(-1, 1) +
  geom_text_repel(aes(label = colnames(dat_plate1)[3:12]), size = 5.5, color = "black") +
  theme(plot.title = element_blank(),
        axis.text = element_text(size = 22), 
        axis.title = element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        legend.title = element_text(size=22),
        legend.text = element_text(size=22),
        aspect.ratio = 6/6) 
dev.off()





###
### Plate 2
###

# Run PCA
pca_plate2 <- prcomp(as.matrix(dat_plate2[3:12]), center = TRUE, scale = TRUE)

# PCA score plot
pca_score_plot <- fviz_pca_ind(pca_plate2,
                               habillage = dat_plate2$Outcome,
                               repel = TRUE,
                               geom = c("point"),
                               pointsize = 4,
                               #alpha = 0.8,
                               invisible = "quali",
                               pch = 19) + 
  xlim(-4.5, 4.5) +
  ylim(-4.5, 4.5) +
  theme_linedraw() +
  scale_color_manual(values = c(nonrestrictive.col, restrictive.col)) +
  geom_text_repel(aes(label = rownames(dat_plate2), color = dat_plate2$Outcome), size = 4.5) +
  theme(plot.title = element_blank(),
        axis.text = element_text(size = 22), 
        axis.title = element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        legend.title = element_blank(),
        legend.text = element_text(size=22),
        legend.position = "top",
        aspect.ratio = 6/6) 
ggsave("pca_plate2_score.eps", plot = pca_score_plot, device = "eps", width = 6, height = 6, units = 'in')


# PCA loading plot
pca_loading_plot <- fviz_pca_var(pca_plate2, 
                                 geom = c("arrow"),
                                 repel = TRUE,
                                 col.var = "grey70") +
  theme_linedraw() +
  xlim(-1, 1) +
  ylim(-1, 1) +
  geom_text_repel(aes(label = colnames(dat_plate2)[3:12]), size = 5.5, color = "black") +
  theme(plot.title = element_blank(),
        axis.text = element_text(size = 22), 
        axis.title = element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        legend.title = element_text(size=22),
        legend.text = element_text(size=22),
        aspect.ratio = 6/6) 
ggsave("pca_plate2_loading.eps", plot = pca_loading_plot, device = "eps", width = 6, height = 6, units = 'in')




##############
###### Heatmap -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##############

#Make heatmap
MWB <- colorRampPalette(c("dodgerblue", "white", "darkmagenta"))

###
### Plate 1
###

# Z-score data
plate1_norm <- dat_plate1[ , 3:(ncol(dat_plate1)-1)]
plate1_norm <- scale(plate1_norm, center = TRUE, scale = TRUE)

# Make heatmep
tiff(file = "plate1_cytokine_heatmap.tiff", width = 7.7, height = 5.12, units = "in", res = 300)
heatmap.2(x = t(plate1_norm), 
          dendrogram = "both", 
          col = MWB, 
          density.info = "none",  
          colsep = 1:nrow(plate1_norm), 
          sepcolor = "lightgrey", 
          rowsep = 1:ncol(plate1_norm), 
          sepwidth = c(0.0010, 0.0010), 
          ColSideColors = as.character(dat_plate1$Color),
          margins = c(7, 8), 
          cexRow = 1.7, 
          cexCol = 1.7, 
          trace = "none",
          key = FALSE, 
          lhei = c(1,6),
          lwid = c(0.4,4),
          srtCol = 45)
dev.off()

# Make legend
tiff(file = "plate1_cytokine_heatmapLegend.tiff", width = 7.7, height = 3, units = "in", res = 300)
heatmap.2(x = t(plate1_norm), 
          dendrogram = "both", 
          col = MWB, 
          density.info = "none",  
          colsep = 1:nrow(plate1_norm), 
          sepcolor = "lightgrey", 
          rowsep = 1:ncol(plate1_norm), 
          sepwidth = c(0.0010, 0.0010), 
          ColSideColors = as.character(dat_plate1$Color),
          #margins = c(7, 8), 
          cexRow = 0.1, 
          cexCol = 0.1, 
          trace = "none",
          key = TRUE, 
          keysize = 3,
          #lhei = c(1,6),
          #lwid = c(0.4,4),
          srtCol = 45)
dev.off()





###
### Plate 2
###

# Z-score data
plate2_norm <- dat_plate2[ , 3:(ncol(dat_plate2)-1)]
plate2_norm <- scale(plate2_norm, center = TRUE, scale = TRUE)

# Make heatmep
pdf(file = "plate2_cytokine_heatmap.pdf", width = 7.7, height = 5.12)
heat_plot <- heatmap.2(x = t(plate2_norm), 
                       dendrogram = "both", 
                       col = MWB, 
                       density.info = "none",  
                       colsep = 1:nrow(plate2_norm), 
                       sepcolor = "lightgrey", 
                       rowsep = 1:ncol(plate2_norm), 
                       sepwidth = c(0.0010, 0.0010), 
                       ColSideColors = as.character(dat_plate2$Color),
                       margins = c(7, 8), 
                       cexRow = 1.7, 
                       cexCol = 1.7, 
                       trace = "none",
                       key = FALSE, 
                       lhei = c(1,6),
                       lwid = c(0.4,4),
                       srtCol = 45)
dev.off()


# Make legend
pdf(file = "plate2_cytokine_heatmapLegend.pdf", width = 7.7, height = 3)
legend_plot <- heatmap.2(x = t(plate2_norm), 
                         dendrogram = "both", 
                         col = MWB, 
                         density.info = "none",  
                         colsep = 1:nrow(plate2_norm), 
                         sepcolor = "lightgrey", 
                         rowsep = 1:ncol(plate2_norm), 
                         sepwidth = c(0.0010, 0.0010), 
                         ColSideColors = as.character(dat_plate2$Color),
                         #margins = c(7, 8), 
                         cexRow = 0.1, 
                         cexCol = 0.1, 
                         trace = "none",
                         key = TRUE, 
                         keysize = 3,
                         #lhei = c(1,6),
                         #lwid = c(0.4,4),
                         srtCol = 45)
dev.off()







