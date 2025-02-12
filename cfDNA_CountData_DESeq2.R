################################################################################################
# Install necessary packages
install.packages("rlang")
################################################################################################

################################################################################################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
################################################################################################

################################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("SARTools")
BiocManager::install("BiocStyle")  # For Bioconductor-style reports
BiocManager::install("pheatmap")  # For heatmaps

#require(SARTools)
################################################################################################

################################################################################################
install.packages("cli")
install.packages("ggplot2")
install.packages("ggforce")
install.packages("rmarkdown")
################################################################################################

################################################################################################
devtools::install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")
################################################################################################

################################################################################################
# Load Required Libraries
library(cli)
library(DESeq2)
library(org.Hs.eg.db)
library(devtools)
library(SARTools)
library(ggplot2)
library(ggforce)
library(pheatmap)
library(rmarkdown)
library(BiocStyle)
library(dplyr) 

################################################################################################
# Set working directory and define file paths
workDir <- getwd()
projectName <- "DESeq2 Data Analysis - cfDNA count Data"
author <- "Sakuntha Gunarathna"
targetFile <- "target.txt"
rawDir <- workDir
inputFile <- "DESeq2 Input.csv"

# Define DESeq2 and Analysis Parameters
varInt <- "group"  # Variable of interest (e.g., Healthy vs Cancer)
alpha <- 0.05
typeTrans <- "VST" 
cpmCutoff <- 0.1  
gene.selection <- "pairwise"
featuresToRemove <- c()  # No features to remove
batch <- NULL  # No batch effect assumed
pAdjustMethod <- "BH"
condRef <- "Healthy"  # Reference group

# Define Colors for Visualization
colors <- c("#0067a5", "#be0032","#e68fac", "#f3c300", 
            "#875692", "#f38400","#a1caf1", "#c2b280",
            "#848482", "#008856","#f99379", "#604e97")

################################################################################################
# Load Count Data (Only Once)
readcounts <- read.delim(file = inputFile, header = TRUE, sep = ",", comment.char = "#", 
                         stringsAsFactors = FALSE, row.names = 1)

# Convert to Matrix for DESeq2
counts <- as.matrix(readcounts)

################################################################################################
# Write Individual Sample Data (if needed)
lapply(names(readcounts), function(x) {
  write.table(cbind(rownames(readcounts), readcounts[[x]]), 
              file = paste("output", x, sep = ""), 
              row.names = FALSE, quote = FALSE, 
              col.names = FALSE, sep = "\t")
})

################################################################################################
# Check Parameters for DESeq2 Analysis
checkParameters.DESeq2(projectName=projectName, author=author, targetFile=targetFile,
                       rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                       condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod,
                       colors=colors, fitType="local", cooksCutoff=TRUE, 
                       independentFiltering=TRUE, typeTrans="rlog", locfunc="median")

################################################################################################
# Load Target File
target_data <- read.delim(targetFile, header = TRUE, stringsAsFactors = TRUE)
print(unique(target_data[[varInt]]))  # Verify groups in the target file

################################################################################################
# Create DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = target_data, design = ~ group)

# Run Differential Expression Analysis
dds <- DESeq(dds)

# Extract Results
res <- results(dds, alpha = alpha)

# Extract Normalized Counts
normalized_counts <- counts(dds, normalized = TRUE)

# Save Normalized Counts
write.csv(normalized_counts, file = "DESeq2_normalized_counts.csv")

########################################################################################################################################
########################################################################################################################################

# Summary of the analysis (boxplots, dispersions, export table, histograms, MA plot)
# Extract results
res <- results(dds, alpha=alpha)

# Summary of differential expression
summary(res)

# Save full results
summaryResults <- as.data.frame(res)
write.csv(summaryResults, file="DESeq2_summary_results.csv")


################################################################################################
# Generate Key DESeq2 Visualizations and Save as PNG & PDF

# Dispersion Plot
png("DispersionPlot.png", width=6, height=6, units="in", res=300)
plotDispEsts(dds)
dev.off()
pdf("DispersionPlot.pdf", width=6, height=6)
plotDispEsts(dds)
dev.off()

# MA Plot
png("MA_Plot.png", width=6, height=6, units="in", res=300)
plotMA(res, ylim=c(-2,2))
dev.off()
pdf("MA_Plot.pdf", width=6, height=6)
plotMA(res, ylim=c(-2,2))
dev.off()

# PCA Plot (Small datasets (≤30 samples) Reduces noise, stabilizes low-count features)
vsd <- rlog(dds, blind=FALSE)
png("PCA_Plot.png", width=6, height=6, units="in", res=300)
plotPCA(vsd, intgroup="group")
dev.off()
pdf("PCA_Plot.pdf", width=6, height=6)
plotPCA(vsd, intgroup="group")
dev.off()

# PCA Plot (Large datasets (≥30 samples)Scales well, avoids excessive shrinkage)
#vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
#png("PCA_Plot.png", width=6, height=6, units="in", res=300)
#plotPCA(vsd, intgroup="group")
#dev.off()
#pdf("PCA_Plot.pdf", width=6, height=6)
#plotPCA(vsd, intgroup="group")
#dev.off()

# Sample-to-Sample Distance Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)


while(dev.cur() > 1) dev.off()
# Generate PNG heatmap
png("Sample_Distance_Heatmap.png", width=6, height=6, units="in", res=300)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)
dev.off()


while(dev.cur() > 1) dev.off()
# Generate PDF heatmap
pdf("Sample_Distance_Heatmap.pdf", width=6, height=6)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)
dev.off()


# Boxplot of Normalized Counts
boxplot_data <- log2(assay(vsd) + 1)
boxplot_data[is.na(boxplot_data)] <- 0  # Replace NaNs with 0

png("Boxplot_Normalized_Counts.png", width=6, height=6, units="in", res=300)
boxplot(boxplot_data, main="Boxplot of Normalized cfDNA Counts", las=2)
dev.off()

pdf("Boxplot_Normalized_Counts.pdf", width=6, height=6)
boxplot(boxplot_data, main="Boxplot of Normalized cfDNA Counts", las=2)
dev.off()


# Histogram of p-values
png("Histogram_pvalues.png", width=6, height=6, units="in", res=300)
hist(res$pvalue, breaks=50, col="skyblue", main="Histogram of p-values", xlab="p-value")
dev.off()
pdf("Histogram_pvalues.pdf", width=6, height=6)
hist(res$pvalue, breaks=50, col="skyblue", main="Histogram of p-values", xlab="p-value")
dev.off()

################################################################################################
#Save Final Summary to a single CSV file

# Read the input CSV file
input_data <- read.csv(inputFile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
normalized_counts <- read.csv("DESeq2_normalized_counts.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
summary_results <- read.csv("DESeq2_summary_results.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Rename the first column to "GeneID" for consistency
colnames(input_data)[1] <- "GeneID"
colnames(normalized_counts)[1] <- "GeneID"
colnames(summary_results)[1] <- "GeneID"

# Merge the datasets based on the 'GeneID' column
merged_data <- input_data %>%
  inner_join(normalized_counts, by = "GeneID") %>%
  inner_join(summary_results, by = "GeneID")

# Remove duplicate suffixes like ".x" and ".y"
colnames(merged_data) <- gsub("\\.x$|\\.y$", "", colnames(merged_data))

# Save the merged dataset to a new CSV file
outputFile <- "Final_DESeq2_Data_Summary.csv"
write.csv(merged_data, file = outputFile, row.names = FALSE)

# Print first few rows to verify merging
print(head(merged_data))


# Save R Session
save.image(file=paste0(projectName, ".RData"))

