library(readr)
library(matrixStats)
library(dplyr)

# Import file 
log_tpm <- read.csv("log_tpm_full.csv", row.names = 1)

# Exclude condition with knockout genes
log_tpm <- log_tpm[, -grep("_del", colnames(log_tpm))]

# Drop columns with isoform
genes_with_isoforms <- grep("_\\d+$", rownames(log_tpm), value = TRUE)

log_tpm <- subset(log_tpm, !(rownames(log_tpm) %in% genes_with_isoforms))


### CODICE BRILLI ###
#load the table containing the regulatory network from regulondb in the form of an edge list
regulator <- read.table(file="tableDataReg.csv",
                  header=TRUE, 
                  sep=",")

# let's remove some of the lines, if they don't refer to a protein regulator:
regulator <- regulator[which(regulator[, 2] != "ppGpp"), ]

#if they are flagged with WEAK confidence
w <- which(trimws(regulator[,7])=="W")
if(length(w)>0){
  regulator <- regulator[-w,]
}

#or if the confidence is actually unkown
w <- which(trimws(regulator[,7])=="?")
if(length(w)>0){
  regulator <- regulator[-w,]
}
#this leaves us with
nrow(regulator)

# Loading files from ECOcyc
map_bnum <- read.delim("mapbnum.txt", header = TRUE)
map_bnum <- map_bnum[c("Gene.Name", "Accession.1")]

# Map between my dataset and the file of ecocyc
log_tpm$gene_number <- rownames(log_tpm)
log_tpm <- merge(log_tpm, map_bnum, by.x = "gene_number", by.y = "Accession.1", all.x = TRUE)

# Rearrange the dataset
log_tpm <- log_tpm[, c("Gene.Name", setdiff(names(log_tpm), "Gene.Name"))]

# Removing unmapped genes bnumber
log_tpm <- subset(log_tpm, !is.na(Gene.Name))

log_tpm_mean <- apply(log_tpm[,3:ncol(log_tpm)], 1, mean)
hist(log_tpm_mean, breaks = 50)

log_tpm_median <- apply(log_tpm[,3:ncol(log_tpm)], 1, median)
hist(log_tpm_median, breaks = 50)

log_tpm_max <- apply(log_tpm[,3:ncol(log_tpm)], 1, max)
hist(log_tpm_max, breaks = 50)

# In this project, we want to introduce a proxy for a transcription factor activity based on the expression level of genes
# regulated by that transcription factor.

# For a positive regulator, its activity can be taken as being a function (linear? sigmoidal?) of
# the expression level of its targets (mean? median? max?).
# Unsupervised learning features = expr level, outcome activity 
