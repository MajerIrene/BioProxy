library(readr)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

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

# Histograms of Summary Statistics
log_tpm_mean <- data.frame(value = apply(log_tpm[,3:ncol(log_tpm)], 1, mean))
mean_hist <- ggplot(log_tpm_mean, aes(x = value)) +
             geom_histogram(binwidth = 0.5, fill = "skyblue", 
                           color = "black", bins = 100) +
             labs(x = "Mean log-TPM", y = "Frequency") +
             theme_minimal() +
             theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())
shapiro.test(log_tpm_mean$value) #not gaussian

log_tpm_median <- data.frame(value = apply(log_tpm[,3:ncol(log_tpm)], 1, median))
median_hist <- ggplot(log_tpm_median, aes(x = value)) +
               geom_histogram(binwidth = 0.5, fill = "lightgreen", 
                              color = "black", bins = 100) +
               labs(x = "Median log-TPM", y = "Frequency") +
               theme_minimal() +
               theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank())
shapiro.test(log_tpm_median$value) #not gaussian

log_tpm_max <- data.frame(value = apply(log_tpm[,3:ncol(log_tpm)], 1, max))
max_hist <- ggplot(log_tpm_max, aes(x = value)) +
               geom_histogram(binwidth = 0.5, fill = "lavender", 
                              color = "black", bins = 100) +
               labs(x = "Max log-TPM", y = "Frequency") +
               theme_minimal() +
               theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank())
shapiro.test(log_tpm_max$value) #not gaussian

log_tpm_min <- data.frame(value = apply(log_tpm[,3:ncol(log_tpm)], 1, min))
min_hist <- ggplot(log_tpm_min, aes(x = value)) +
            geom_histogram(binwidth = 0.5, fill = "lightpink", 
                           color = "black", bins = 100) +
            labs(x = "Min log-TPM", y = "Frequency") +
            theme_minimal() +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())
shapiro.test(log_tpm_min$value) #not gaussian

grid.arrange(mean_hist, median_hist, max_hist, min_hist, nrow = 2, ncol = 2,
             top = textGrob("Histograms of Summary Statistics", 
                            gp=gpar(fontsize=16)))

# Histogram of how many genes are regulated by each gene
positive_reg <- regulator[regulator$X6.function == "+",]
negative_reg <- regulator[regulator$X6.function == "-",]
unique_regulators <- unique(regulator$X3.RegulatorGeneName)

pos_counts <- c()
neg_counts <- c()

for(reg in unique_regulators){
  pos_counts <- c(pos_counts, 
                  count(positive_reg[positive_reg$X3.RegulatorGeneName == reg,]))
  neg_counts <- c(neg_counts, 
                  count(negative_reg[negative_reg$X3.RegulatorGeneName == reg,]))
}

pos_counts <- unlist(pos_counts)
names(pos_counts) <- unique_regulators

neg_counts <- unlist(neg_counts)
names(neg_counts) <- unique_regulators

pos_counts <- data.frame(value = pos_counts)
shapiro.test(pos_counts$value) #not gaussian

neg_counts <- data.frame(value = neg_counts)
shapiro.test(neg_counts$value) #not gaussian

total_counts <- data.frame(value = pos_counts$value + neg_counts$value)
shapiro.test(total_counts$value) #not gaussian

pos_counts_hist <- ggplot(pos_counts, aes(x = value)) +
                   geom_histogram(binwidth = 10, fill = "skyblue", 
                                  color = "black", bins = 20) +
                   labs(x = "Positive Regulations Count", y = "Frequency") +
                   theme_minimal() +
                   theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank())

neg_counts_hist <- ggplot(neg_counts, aes(x = value)) +
                   geom_histogram(binwidth = 10, fill = "lightgreen", 
                                  color = "black", bins = 20) +
                   labs(x = "Negative Regulations Count", y = "Frequency") +
                   theme_minimal() +
                   theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank())

total_counts_hist <- ggplot(total_counts, aes(x = value)) +
                   geom_histogram(binwidth = 10, fill = "lightpink", 
                                  color = "black", bins = 20) +
                   labs(x = "Total Regulations Count", y = "Frequency") +
                   theme_minimal() +
                   theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank())

grid.arrange(pos_counts_hist, neg_counts_hist, total_counts_hist, nrow = 1)

# crp positively controls 310 genes, lets do a subset
crp_target <- positive_reg[positive_reg$X3.RegulatorGeneName == "crp",]
crp_target_exp <- log_tpm[log_tpm$Gene.Name %in% crp_target$X5.regulatedName,]
crp_target_exp$mean_exp <- apply(crp_target_exp[,3:ncol(crp_target_exp)], 1, mean)

# caccapupu <- apply(crp_target_exp[,3:(ncol(crp_target_exp)-1)], 2, mean)
# plot(caccapupu, log_tpm[log_tpm$Gene.Name == "crp", 3:ncol(log_tpm)])

# Scaling with RobustScaler
# log_tpm_norm <- log_tpm
# log_tpm_norm[,3:ncol(log_tpm_norm)] <- 
#   apply(log_tpm_norm[,3:ncol(log_tpm_norm)], 1, 
#         function(x){(x-median(x))/IQR(x)})
