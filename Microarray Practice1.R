# importing libraries
library(GEOquery)
library(tidyverse)
library(affy)

# Get supplemantary files
getGEOSuppFiles("GSE148537")

# Untar files
untar("GSE148537/GSE148537_RAW.tar")

# read in .cel files
raw.data <- ReadAffy(celfile.path = "Data/")
raw.data

# now we use RMA normalization method to our data
normalize.data <- rma(raw.data)

# get expression estimates
normalized.expr <- exprs(normalize.data)

# Change the generated matrix to dataframe
normalized.expr <- as.data.frame(normalized.expr)

# map probe IDs to gene symboles
gse <- getGEO("GSE148537", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbole mapping
feature.data <- gse$GSE148537_series_matrix.txt.gz@featureData@data
# generate a subset for gene symboles
feature.data <- feature.data[, c(1,11)]

# merge feature.data to normalized.expr dataframes
normalized.expr <- normalized.expr %>%
  rownames_to_column(var = "ID") %>%
  inner_join(., feature.data, by= "ID")

# now we have our final expession values for visualization