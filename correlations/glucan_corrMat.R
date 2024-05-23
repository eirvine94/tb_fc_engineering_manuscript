#
# Author: Edward B. Irvine
#
# Description: 
# Correlation matrix of down-selected 24c5 Fc-variant panel.
#
# Created: 29 November 2021
#
# Modified: 23 May 2024
#

##################
##### Housekeeping -------------------------------------------------------------------------------------------------------------------------------------------------------
##################

# Load required libraries
library(corrplot)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Read in data
dat <- read.csv("24c5 Functions and Anti-microbial.csv")










######################
##### Pre-process data ---------------------------------------------------------------------------------------------------------------------------------------------------
######################

# Clean-up column names
colnames(dat) <- c("Variant", "NK cell (CD107a)", "NK cell (MIP-1β)", "NK cell (IFNγ)", "Monocyte Phago",
                   "Neutrophil Phago", "Complement", "MΦ Restriction", "Whole-blood assay")

# Drop missing data
dat <- na.omit(dat)

# Generate correlation matrix
corr_mat <- cor(dat[1:nrow(dat), 2:ncol(dat)], method = "spearman")

corr_p_mat <- cor_pmat(dat[1:nrow(dat), 2:ncol(dat)], method = "spearman", exact = FALSE)
corr_p_mat








########################
##### Correlation matrix -------------------------------------------------------------------------------------------------------------------------------------------------
########################

# Set color palette
pal <- colorRampPalette(c("#007F00", "white", "#F73772"))(200)


# Plot
#diag(corr_mat) = NA
pdf("glucan_corrMat_v2.pdf", width = 6, height = 6)
corrplot(corr_mat, order = 'AOE', type="upper", method = "ellipse", p.mat = corr_p_mat, tl.pos="lt", sig.level = c(0.0001, 0.001, 0.01, 0.05), insig = 'label_sig', pch.cex = 2, col = pal, tl.col = 'black', na.label = "  ")
corrplot(corr_mat, order = 'AOE', type="lower", method = "number", add=T, tl.pos="lt", cl.pos="n", col = pal, tl.col = 'black', na.label = "  ")
dev.off()






#############################
##### Individual correlations --------------------------------------------------------------------------------------------------------------------------------------------
#############################

colnames(dat) <- c("Variant", "NK cell (CD107a)", "NK cell (MIP-1β)", "NK cell (IFNγ)", "Monocyte Phago",
                   "Neutrophil Phago", "Complement", "MΦ Restriction", "Whole blood assay")

# Generate WBA vs. ADNP plot
my_text <- paste("rho = -0.532; p-value = 0.044", sep = "")
plot <- ggscatter(data = dat,
                  x = "Neutrophil Phago",
                  y = "Whole blood assay",
                  size = 3,
                  add = "reg.line", 
                  add.params = list(color = "#007F00", fill = "lightgray"),
                  conf.int = TRUE)
plot2 <- plot + 
  ggtitle(my_text) +
  xlab("Neutrophil phagocytosis") +
  ylab("Whole-blood assay") +
  theme_linedraw() +
  theme(plot.title = element_text( size = 18),
        axis.text = element_text(size = 14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        axis.title = element_text(size = 18),
        legend.position = "none",
        aspect.ratio = 5/6.5) 
plot2
ggsave("WBA_ADNP_corr_v2.pdf", plot = plot2, device = "pdf", width = 6.5, height = 5, units = 'in')



# Generate WBA vs. ELISA plot
wba_elisa_dat <- read.csv("230926_glucan_wba_elisa_corr.csv")
colnames(wba_elisa_dat) <- c("Variant", "Whole blood assay", "ELISA")

wba_elisa_dat <- na.omit(wba_elisa_dat)
cor(x = wba_elisa_dat$ELISA, y = wba_elisa_dat$`Whole blood assay`, method = "spearman")
my_text <- paste("rho = 0.100; p-value = 0.724", sep = "")
plot <- ggscatter(data = wba_elisa_dat,
                  x = "ELISA",
                  y = "Whole blood assay",
                  size = 3,
                  add = "reg.line", 
                  add.params = list(color = "#007F00", fill = "lightgray"),
                  conf.int = TRUE)
plot3 <- plot + 
  ggtitle(my_text) +
  xlab("α-glucan ELISA") +
  ylab("Whole-blood assay") +
  theme_linedraw() +
  theme(plot.title = element_text( size = 18),
        axis.text = element_text(size = 14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        axis.title = element_text(size = 18),
        legend.position = "none",
        aspect.ratio = 5/6.5) 
plot3
ggsave("WBA_ELISA_corr.pdf", plot = plot3, device = "pdf", width = 6.5, height = 5, units = 'in')




