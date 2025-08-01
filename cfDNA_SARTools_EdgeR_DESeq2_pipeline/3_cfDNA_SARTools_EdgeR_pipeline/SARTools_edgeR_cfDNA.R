################################################################################################
# Install Required Packages
################################################################################################

install.packages("rlang")        # Needed for tidy evaluation
install.packages("cli")          # Command-line interface enhancements
install.packages("ggplot2")      # For plotting
install.packages("ggforce")      # For adding ellipses to MDS plots

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")        # Although this is DESeq2, SARTools supports both DESeq2 and edgeR
BiocManager::install("org.Hs.eg.db")  # Human genome annotation package

devtools::install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")

################################################################################################
# Load Required Libraries
################################################################################################

library(cli)
library(DESeq2)
library(org.Hs.eg.db)
library(devtools)
require(SARTools)
library(SARTools)
library(ggplot2)
library(ggforce)

################################################################################################
# Load Count Data and Targets
################################################################################################

# This file should be the output of the Python script:
# "Count_Matrix.txt" (generated from deepTools multiBigwigSummary results)
readcounts <- read.delim(file = "Count_Matrix.txt", 
                         header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)

# Set Gene IDs as row names
rownames(readcounts) <- readcounts$Gene.ID

# Drop the Gene ID column (since it's now used as rownames)
readcounts <- readcounts[,c(2:ncol(readcounts))]

# Export counts per sample into individual text files (required for SARTools)
lapply(names(readcounts), function(x) {
  write.table(cbind(rownames(readcounts), readcounts[[x]]), 
              file = paste("output_", x, ".txt", sep = ""), 
              row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
})

################################################################################################
# Project Parameters
################################################################################################

workDir <- getwd()                        
projectName <- "SARTools_cfDNA_Analysis"  
author <- "Sakuntha Devaka Gunarathna"    
targetFile <- "target.txt"                
rawdir <- getwd()                         
varInt <- "group"                         
alpha <- 0.05                             
typeTrans <- "VST"                        
cpmCutoff <- 0.1                          
gene.selection <- "pairwise"              
normalizationMethod <- "TMM"              

colors <- c("#0067a5", "#be0032", "#e68fac", "#f3c300", 
            "#875692", "#f38400", "#a1caf1", "#c2b280",
            "#848482", "#008856", "#f99379", "#604e97")

featuresToRemove <- c()        
batch <- NULL                  
pAdjustMethod <- "BH"          
condRef <- "Healthy"           

################################################################################################
# Check Input Parameters
################################################################################################

checkParameters.edgeR(projectName=projectName, author=author, targetFile=targetFile,
                      rawDir=rawdir, featuresToRemove=featuresToRemove, varInt=varInt,
                      condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff, gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod, colors=colors)

################################################################################################
# Load Target File (sample metadata)
################################################################################################

target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)
unique(target[[varInt]])  

################################################################################################
# Load Counts
################################################################################################

counts <- loadCountData(target=target, rawDir=rawdir, featuresToRemove=featuresToRemove)

################################################################################################
# QC and Descriptive Plots
################################################################################################

majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

################################################################################################
# Run edgeR Differential Analysis
################################################################################################

out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

################################################################################################
# MDS + Clustering
################################################################################################

exploreCounts(object=out.edgeR$dge, group=target[,varInt], gene.selection=gene.selection, col=colors)

################################################################################################
# Custom MDS Plot (With Ellipses and Labels)
################################################################################################

createMDSPlotWithEllipses <- function(dge, group, col, label_font_size = 1) {
  mds <- plotMDS(dge, plot=FALSE)
  mds_df <- data.frame(x=mds$x, y=mds$y, group=group, sample=rownames(dge$samples))

  p <- ggplot(mds_df, aes(x=x, y=y, color=group, label=sample)) +
    geom_point(size=2) +
    geom_mark_ellipse(aes(group=group, fill=group), alpha=0.2) +
    geom_text(aes(label=sample), hjust=0.5, vjust=-2, size=label_font_size) +
    scale_color_manual(values=col) +
    scale_fill_manual(values=col) +
    theme_minimal() +
    labs(title="Multi-Dimensional Scaling Plot", 
         x="Leading logFC dimension 1", 
         y="Leading logFC dimension 2") +
    theme(legend.position="right") +
    theme(plot.margin = margin(20, 20, 20, 20))

  return(list(plot=p, mds_df=mds_df))
}

result <- createMDSPlotWithEllipses(out.edgeR$dge, target[,varInt], colors)
p <- result$plot
mds_df <- result$mds_df

write.csv(mds_df, file="MDS_Data.csv", row.names=FALSE)
ggsave(filename="MDS_Plot_with_Ellipses.png", plot=p, width=6, height=6, units="in")
ggsave(filename="MDS_Plot_with_Ellipses.pdf", plot=p, width=6, height=6, units="in")

################################################################################################
# MDS Plot Without Labels
################################################################################################

createMDSPlotNoLabels <- function(dge, group, col) {
  mds <- plotMDS(dge, plot=FALSE)
  mds_df <- data.frame(x=mds$x, y=mds$y, group=group, sample=rownames(dge$samples))

  p <- ggplot(mds_df, aes(x=x, y=y, color=group)) +
    geom_point(size=3) +
    geom_mark_ellipse(aes(group=group, fill=group), alpha=0.1) +
    scale_color_manual(values=col) +
    scale_fill_manual(values=col) +
    theme_minimal() +
    labs(title="MDS Plot (No Labels)", 
         x="Leading logFC dimension 1", 
         y="Leading logFC dimension 2") +
    theme(legend.position="right") +
    theme(plot.margin = margin(20, 20, 20, 20))

  return(list(plot=p, mds_df=mds_df))
}

result2 <- createMDSPlotNoLabels(out.edgeR$dge, target[,varInt], colors)
p2 <- result2$plot
ggsave(filename="MDS_Plot_NoLabels.png", plot=p2, width=6, height=6, units="in")
ggsave(filename="MDS_Plot_NoLabels.pdf", plot=p2, width=6, height=6, units="in")

################################################################################################
# Summarize Results and Generate Report
################################################################################################

summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=counts, alpha=alpha, col=colors)

writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults, 
                  cpmCutoff=cpmCutoff, majSequences=majSequences, workDir=workDir, 
                  projectName=projectName, author=author, targetFile=targetFile, rawDir=rawdir, 
                  featuresToRemove=featuresToRemove, varInt=varInt, condRef=condRef, batch=batch, 
                  alpha=alpha, pAdjustMethod=pAdjustMethod, colors=colors,
                  gene.selection=gene.selection, normalizationMethod=normalizationMethod)
