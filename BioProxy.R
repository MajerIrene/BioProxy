library(readr)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(caret)
library(glmnet)
library(sigmoid)
library(vip)

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

#or if the confidence is actually unknown
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

#removing duplicate genes (it also has all expression values equal 0 so very bad)
log_tpm <- subset(log_tpm, !(log_tpm$Gene.Name == "insI2"))

#setting rownames and dropping the first 2 columns
rownames(log_tpm) <- log_tpm$Gene.Name
log_tpm <- log_tpm[,3:ncol(log_tpm)]

#transpose log_tpm
log_tpm <- t(log_tpm)

# Histograms of Summary Statistics
log_tpm_mean <- data.frame(value = apply(log_tpm, 2, mean))
mean_hist <- ggplot(log_tpm_mean, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", 
                 color = "black", bins = 100) +
  labs(x = "Mean log-TPM", y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
shapiro.test(log_tpm_mean$value) #not gaussian

log_tpm_median <- data.frame(value = apply(log_tpm, 2, median))
median_hist <- ggplot(log_tpm_median, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "lightgreen", 
                 color = "black", bins = 100) +
  labs(x = "Median log-TPM", y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
shapiro.test(log_tpm_median$value) #not gaussian

log_tpm_max <- data.frame(value = apply(log_tpm, 2, max))
max_hist <- ggplot(log_tpm_max, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "lavender", 
                 color = "black", bins = 100) +
  labs(x = "Max log-TPM", y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
shapiro.test(log_tpm_max$value) #not gaussian

log_tpm_min <- data.frame(value = apply(log_tpm, 2, min))
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

# crp positively controls 310 genes
crp_exp <- log_tpm[,colnames(log_tpm) == "crp"]
crp_target <- positive_reg[positive_reg$X3.RegulatorGeneName == "crp",]
crp_target_exp <- log_tpm[,colnames(log_tpm) %in% crp_target$X5.regulatedName]
#crp_target_exp$mean_exp <- apply(crp_target_exp[,3:ncol(crp_target_exp)], 1, mean)

crp_target_mean_conditions <- apply(crp_target_exp, 1, mean)

training <- data.frame(crp_exp = unlist(crp_exp), target = unlist(crp_target_mean_conditions))

fit_control <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

crp_lm <- train(crp_exp ~ target, 
                data = training,
                method = "lm",
                trControl = fit_control)

scatterplot <- ggplot(training, aes(x = target, y = crp_exp)) +
  geom_point(size = 2) + 
  geom_abline(intercept = crp_lm$finalModel$coefficients[[1]],
              slope = crp_lm$finalModel$coefficients[[2]],
              color = "red",
              linewidth = 1)  

crp_res <- data.frame(fitted.values = crp_lm$finalModel$fitted.values,
                      residuals = crp_lm$finalModel$residuals)

residuals <-  ggplot(crp_res, aes(x = fitted.values, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red")

#try to fit with 300 features
crp_exp <- data.frame(crp = crp_exp)
crp_full <- data.frame(cbind(crp_exp, crp_target_exp))

crp_lm2 <- train(crp ~., 
                 data = crp_full,
                 method = "lm",
                 preProcess = c("center", "scale", "pca"),
                 trControl = fit_control)

crp_lm2$finalModel$coefficients

crp_res2 <- data.frame(fitted.values = crp_lm2$finalModel$fitted.values,
                      residuals = crp_lm2$finalModel$residuals)

residuals <-  ggplot(crp_res2, aes(x = fitted.values, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red")

crp_prova <- data.frame(fitted.values = crp_lm2$finalModel$fitted.values,
                        real_values = crp_lm2$finalModel$model$.outcome)

prova <-  ggplot(crp_prova, aes(x = fitted.values, y = real_values)) +
          geom_point() 

prova2 <- ggplot(crp_full, aes(x = gpsA, y = crp)) +
          geom_point() +
          geom_abline(intercept = 12.2460, slope = -0.2683, color = "red")

acab <- lm(crp ~ gpsA + fucR + malY + raiA + frlR, data = crp_full)

plot(crp_lm2$finalModel)


#to analize the pca contribution of each gene in the final model with 64 PC
pca_results <- data.frame(crp_lm2$preProcess$rotation)

#pca in another way, this also consider up to 188 PC
crp_preProc <- preProcess(crp_full, method = c("center", "scale", "pca"))
crp_pca <- predict(crp_preProc, crp_full)

plot(crp_pca[,1], crp_pca[,2])
plot(crp_pca[,2], crp_pca[,3])
plot(crp_pca[,3], crp_pca[,4])

#lasso regression (L1 regularization)
crp_lasso <- train(crp ~., 
                   data = crp_full,
                   method = "lasso",
                   preProcess = c("center", "scale", "pca"),
                   trControl = fit_control)

crp_lasso$finalModel$beta.pure

print(crp_lasso)

plot(crp_lasso$finalModel, xvar = "penalty")

### PCR 
cv_model_pcr <- train(crp ~., 
                      data = crp_full, 
                      method = "pcr",
                      trControl = fit_control,
                      preProcess = c("zv", "center", "scale"), 
                      tuneLength = 160)

print(cv_model_pcr)

cv_model_pcr$bestTune

cv_model_pcr$results %>%
        dplyr::filter(ncomp == pull(cv_model_pcr$bestTune))
        ggplot(cv_model_pcr)

prova_vi <- vi(cv_model_pcr, num_features = 50, method = "model")

formula1 <- paste0("crp  ~ ", paste(prova_vi$Variable[1:70], collapse = " + "))
formula2 <- paste0("crp  ~ ", paste(prova_vi$Variable[1:100], collapse = " + "))
formula3 <- paste0("crp  ~ ", paste(prova_vi$Variable[1:50], collapse = " + "))

crp_proxy1 <- train(as.formula(formula1), 
                   data = crp_full,
                   method = "lm",
                   trControl = fit_control)

crp_proxy2 <- train(as.formula(formula2), 
                   data = crp_full,
                   method = "lm",
                   trControl = fit_control)

crp_proxy3 <- train(as.formula(formula3), 
                    data = crp_full,
                    method = "lm",
                    trControl = fit_control)

anova(crp_proxy3$finalModel, crp_proxy2$finalModel, crp_proxy1$finalModel)
