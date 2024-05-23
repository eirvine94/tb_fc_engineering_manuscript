#
# Author: Edward B. Irvine
#
# Description: Polar plots of down-selected 24c5 Fc-variant panel.
#
# Created: 29 November 2021
#
# Modified: 23 May 2024
#

##################
##### Housekeeping -------------------------------------------------------------------------------------------------------------------------------------------------------
##################

# Load required libraries
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# Read in data
dat <- read.csv("glucan_data_maxNorm.csv")










######################
##### Define functions -------------------------------------------------------------------------------------------------------------------------------------------------
######################

# Polar plot
polar_plot <- function(dat, variant, pal) {
  
  # Name plot file
  eps_name <- paste(variant, "_polarPlot.eps", sep = "")

  # Generate plot
  plot <- dat %>%
    ggplot(aes(x = rownames(dat), y = UQ(as.name(variant)), fill = rownames(dat))) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    theme_minimal() +
    coord_polar() +
    ggtitle(paste(variant)) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(limits = c(0,1)) +
    theme(panel.grid.major = element_line(colour = "black", linetype = "dotted", size = 0.25),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 22, hjust = 0.5, vjust = -2))
  
  # Save plot
  ggsave(filename = eps_name, plot = plot, device = "eps", width = 5, height = 5, units = "in", dpi = 300)
  return(plot)
}

# Legend plot
legend_plot <- function(dat, pal) {
  
  # Make legend data frame
  legend_frame <- data.frame(c(1, 1, 1, 1, 1, 1))
  colnames(legend_frame) <- "legend"
  rownames(legend_frame) <- rownames(dat)
  
  # Generate plot
  plot <- legend_frame %>%
    ggplot(aes(x = rownames(legend_frame), y = legend, fill = rownames(legend_frame))) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    theme_minimal() +
    coord_polar() +
    scale_fill_manual(values = pal) +
    scale_y_continuous(limits = c(0,1)) +
    theme(panel.grid.major = element_line(colour = "white", linetype = "dotted"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 22, hjust = 0.5, vjust = -5))
  
  # Save plot
  ggsave(filename = "legend_polarPlot.eps", plot = plot, device = "eps", width = 5, height = 5, units = "in", dpi = 300)
  return(plot)
}










######################
##### Pre-process data ---------------------------------------------------------------------------------------------------------------------------------------------------
######################

# Clean-up column names
colnames(dat) <- c("Variant", "NK (CD107a)", "NK (MIP-1β)", "NK (IFNγ)", "Monocyte Phago", "Neutrophil Phago", "Complement")

# Keep only down-selected variants
keep <- c("IgG1", "IgG2", "IgG4", "LALA", "N297Q", "YTE", "E380A", "MLNS",
          "IgG3 RH", "SAEAKA", "HFST", "SEHFST LS", "I332E", "SDIE", "SDIEAL")
dat_selected <- dat %>% 
  slice(match(keep, Variant))

# Transpose data frame
dat_selected <- t(dat_selected)
colnames(dat_selected) <- dat_selected["Variant", ]
dat_selected <- dat_selected[2:nrow(dat_selected), ]
features <- rownames(dat_selected)
dat_selected <- data.frame(apply(dat_selected, 2, function(x) as.numeric(as.character(x))))
rownames(dat_selected) <- features










#################
##### Polar plots ------------------------------------------------------------------------------------------------------------------------------------------------------
#################

# Set color palette
pal <- c("#CCEBC5", "#FBB4AE", "#B3CDE3", "#DECBE4", "#FFFFCC", "#FED9A6")

# Make polar plot for each alpha-glucan Fc-variant
plot_list = list()
for (var in colnames(dat_selected)) {
  plot <- polar_plot(dat = dat_selected, variant = var, pal = pal)
  plot_list[[var]] = plot
}

# Plot legend
legend <- legend_plot(dat = dat_selected, pal = pal)




