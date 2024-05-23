#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Edward B. Irvine

Description:
Performs hierarchical clustering on 24c5 Fc-variant panel functional data.

Created: 14 April 2021
    
Modifed: 23 May 2024

"""


##################
##### Housekeeping ----------------------------------------------------------------------------------------------------------------------------------------
##################

# Import required libraries
import pandas as pd
from scipy.cluster import hierarchy
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# Import data
data = pd.read_csv("24c5_hclust_input_dat.csv")










######################
##### Pre-process data ------------------------------------------------------------------------------------------------------------------------------------
######################

### Clean data
data = data.rename(index = data["Unnamed: 0"])
data = data.iloc[:, 1:-1]


### Z-score data
# Create the object
scaler = StandardScaler()

# Calculate the mean and standard deviation
scaler.fit(data)

# Transform the values
data_scaled = pd.DataFrame(scaler.transform(data))

# Restore row and column names
data_scaled = data_scaled.rename(columns = pd.Series(data.columns), index = pd.Series(data.index))

# Drop rows w/ missing values
data_scaled = data_scaled.dropna(axis = 0)

# Drop additional rows
data_scaled = data_scaled.drop(["LPLIL LS", "SDIESA LS", "SDIEGA LS"])








#############################
##### Hierarchical clustering ------------------------------------------------------------------------------------------------------------------------------
#############################

###
### Complete linkage
###

# Hierarchically cluster data
hclust_data_complete = hierarchy.linkage(data_scaled, method = "complete")
hclust_data_complete

# Plot dendrogram
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.figure(figsize = (4, 12))
dendrogram_complete = hierarchy.dendrogram(hclust_data_complete, 
                                      labels = data_scaled.index, 
                                      orientation = "left", 
                                      leaf_font_size = 12, 
                                      color_threshold = 2.8)
#plt.title("Cluster Dendrogram", fontsize = 15, fontweight = "bold")
plt.xlabel('Euclidean Distance', fontsize = 15)
plt.tick_params(axis='x', which='major', labelsize=15)
plt.tight_layout()
plt.savefig('24c5_hclust_dendrogram.eps', dpi = 300, format = "eps")







