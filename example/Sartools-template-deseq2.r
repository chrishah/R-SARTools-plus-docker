################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### November 28th, 2019
### designed to be executed with SARTools 1.7.3
################################################################################
# Only the first time ever to install
#install.packages("devtools")
#library(devtools) 
#devtools::install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")

library(SARTools)
library(dplyr)

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session

workDir <- "/home/sbi6dap/mydata/RNAseq-Analysis"      # *** working directory for the R session

projectName <- "projectName"                         # *** name of the project
author <- "Your name"                                # *** author of the statistical analysis/report

targetFile <- "metadata.txt"                           # *** path to the design/target file
rawDir <- "FeatureCounts/rmdup"                                      # *** path to the directory containing raw counts files

varInt <- "Type"                                    # *** factor of interest
condRef <- "Control"                                      # *** reference biological condition
batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example

featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default), "local" or "mean"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

colors <- c("#f3c300", "#875692", "#f38400",         # vector of colors of each biological condition on the plots
            "#a1caf1", "#be0032", "#c2b280",
            "#848482", "#008856", "#e68fac",
            "#0067a5", "#f99379", "#604e97")

forceCairoGraph <- FALSE

IDtoGene <- read.delim("/mnt/clusters/sponsa/data/classdata/Bioinformatics/REFS/Homo_sapiens.GRCh38.100.map.txt", header=T)

################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)

if (forceCairoGraph) options(bitmapType="cairo")

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# save image of the R session
#save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)

# Add coulmmn to tables with gene name
tables <- list.files(path="tables", pattern="*.txt", full.names=TRUE, recursive=FALSE)
lapply(tables, function(x) {
  t <-read.table(x, header=TRUE)
  t.cols <- select(t, Id, log2FoldChange, pvalue, padj)
  t.annot <- merge(t.cols, IDtoGene, by='Id')
  write.table(t.annot, file=paste(x, ".annot.csv", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
})



####################################################################
####################################################################
###### Making Heatmaps #############################################
####################################################################
####################################################################



#################
## Input files ##
#################

# Files must be in the format with a header (samples) and gene names (rows)
# Make sure to use the normalised data!
infile = data.frame(read.table("tables/DeletionvsControl.complete.txt", sep="\t", header = TRUE, row.names= 1))
infile <- infile[,15:28]
infile<-infile[!is.na(infile$p.adj),]
infile<-infile[order(infile$p.adj),]
x<-infile[1:100,]

# Metadata column names must match sample names!
x.meta = data.frame(read.table(targetFile, header = TRUE, row.names= 1))

# Ensure ordering is correct
y <- x[ , order(names(x))]
y.meta <- x.meta[order(row.names(x.meta)), ]

# Convert the data into gene-scaled format 
y.scaled = t(apply(y, 1, function(x) {
  q10 = quantile(x, 0.1)
  q90 = quantile(x, 0.9)
  x[x < q10] = q10
  x[x > q90] = q90
  scale(x)
}))
colnames(y.scaled) = colnames(y)


## NOTE: You now have two datasets: 
##         -- y        (raw normalised count data)
##         -- y.scaled (thes same data, but each gene is on it's own scale)
##
## Choose which to plot yourself!

###################
## Basic Heatmap ##
###################
library(RColorBrewer)
library(gplots)

# This function does a basic heatmap based on the vaules alone, and can cluster together samples and genes

# Choosing colours. You can ignore this and delete the "col = my_palette" line in the function to use the default instead
my_palette <- colorRampPalette(brewer.pal(9,"Reds"))

# Here I am logging the data, and adding 1 to ensure not logging a zero value
heatmap.2(data.matrix(y.scaled),
          Rowv=TRUE,                    ## Should the script cluster on rows?
          Colv=TRUE,                    ## Should the script cluster on columns?
          col = my_palette,             ## Colour choice
          cexRow=0.5,
          margins=c(5,7),
          xlab = "samples", ylab = "genes", main = "My Significant Genes",     ## Labels!
          trace="none",
          density.info=c("none")
)

#####################
## Complex Heatmap ##
#####################

library(ComplexHeatmap)

# This function can use the metadata file to split
# your samples up, or add different annotations

# Set the factor to split by:
factor <- as.vector(y.meta$Type)

# Plot the heatmap
Heatmap(data.matrix(y.scaled),
        cluster_columns = TRUE,
        #col = colorRampPalette(c("white", "black", "red"))(50),
        top_annotation = columnAnnotation(abc = factor),
        row_names_gp = gpar(fontsize = 8),
        column_split = factor,
        cluster_column_slices = FALSE,
        #show_row_names=T,
        #show_column_names=T,
        #border = TRUE
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson"
)
