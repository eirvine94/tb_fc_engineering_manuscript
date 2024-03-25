#####
#
# Author: Edward B. Irvine
#
# Description:
# This script performed Gene Ontology analysis in CD14 monocystes. 
# Genes upregulated in 24c5 SEHFST LS over both 24c5 IgG1 and the no Ab condition were used for this enrichment analysis.
# Biological process data associated with different conditions from a gene ontology (GO) analysis are imported. 
# It then converts the adjusted p-values from these datasets into their negative logarithms for better visualization. 
# The main function, lollipop_plot, is designed to create lollipop plots, specifically for displaying the significance and effect size (as odds ratios) of different biological processes based on GO terms. 
# This function is applied to create lollipop plots for each set of data, representing upregulated and downregulated biological processes, with customizations including color, dot size, and font settings for clarity and emphasis. 
# These plots are then saved as TIFF files, providing a clear, graphical representation of the enrichment analysis results.
#
# Created: 4/5/22
#
# Modified: 3/25/24
#
#####

##################
##### Housekeeping ----------------------------------------------------------------------------------------------------------------------------------
##################

# Load required packages
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

# Import GO Biological Processes data
red_quad_bp <- read.csv("input_data/enrichr_cd14mono_redQuad_gobp_v3.csv")
blue_quad_bp <- read.csv("input_data/enrichr_cd14mono_blueQuad_gobp_v3.csv")

# Create -log10 p-value feature
red_quad_bp$y <- -log10(red_quad_bp$Adjusted.P.value)
blue_quad_bp$y <- -log10(blue_quad_bp$Adjusted.P.value)










###############################
##### Lollipop enrichment plots ---------------------------------------------------------------------------------------------------------------------
###############################

# Function to generate lollipop enrichment plots
lollipop_plot <- function(dat, gene_set, title, font_size, dot_size, positive) {
  
  # Subset on desired gene set
  #dat <- dat[dat$Source == gene_set, ]
  
  print(dat)
  
  dat$Odds.Ratio <- round(dat$Odds.Ratio, 2)
  if (gene_set %in% c("gobp", "gomf", "gocc")) {
    dat$Term <- str_sub(dat$Term, 1, nchar(dat$Term) - 13)
  }
  
  if (positive == TRUE) {
    direction <- "up"
    col <- "#941100"
  } else {
    direction <- "down"
    col <- "#005493"
  }
  
  # Generate plot
  pop_plot <- dat %>%
    ggdotchart(x = "Term", 
               y = 'y', 
               label = "Odds.Ratio",
               color = col,
               #palette = pal,
               sorting = "descending", 
               add = "segments", 
               add.params = list(color = col, size = 1, alpha = 0.65),
               rotate = TRUE,
               dot.size = dot_size,
               font.label = list(color = "white", size = 9.5, vjust = 0.5, face = "bold"),
               xlab = FALSE) +
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    geom_hline(yintercept = -log10(0.05), color = "black", size = 0.75, linetype = "dashed", alpha = 0.65) +
    theme_linedraw() +
    ggtitle(title) +
    ylab("-log10(padj)") +
    #scale_y_continuous(limits = c(0, 2.4)) +
    theme(axis.text.x = element_text(size = font_size), 
          axis.text.y = element_text(size = font_size), 
          axis.title.x = element_text(size = font_size),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
          plot.title = element_text(size = font_size + 2, face = "bold"),
          legend.position = "none",
          aspect.ratio = 5.75/4.5) 
  pop_plot
  
  # Save plot
  tiff_name <- paste("plots/24c5/mye/enrichr/", gene_set, "_", direction, "_lollipop_cd14mono_v3.tiff", sep = "")
  ggsave(filename = tiff_name, plot = pop_plot, width = 8.7, height = 5.5)
  
  return(pop_plot)
}





# Generate plots
gobp_lollipop_up <- lollipop_plot(dat = red_quad_bp, gene_set = "gobp", title = "GO | Biological Processes", font_size = 10, dot_size = 12, positive = TRUE)
gobp_lollipop_down <- lollipop_plot(dat = blue_quad_bp, gene_set = "gobp", title = "GO | Biological Processes", font_size = 10, dot_size = 12, positive = FALSE)


