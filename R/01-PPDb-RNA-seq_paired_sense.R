##########################################################################
# RNA-seq analysis: PPD-b stimulated vs unstimulated peripheral blood    #
#                           paired-end reads.                            #
#                                                                        #
#           --- R workflow for analyses of known sense genes ---         #
##########################################################################

# Based on the workflow created by Nalpas, N.C. (2014)
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474

# Author of current version (4.0.0): Correia, C.N.
# DOI badge of current version:
# Last updated on 01/08/2017

# R version 3.4.0 (2017-04-21) -- "You Stupid Darkness"
# Bioconductor version 3.3 (BiocInstaller 1.26.0)

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/PPD-B-sti_vs_uns/edgeR_sense")
getwd()

# Define variables for working and file directories
workDir <- getwd()
workDir
fileDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/PPD-B-sti_vs_uns/edgeR_sense/PPDbCounts"

# Load previously saved data
load("PPDb-RNA-seq_paired_sense.RData")

############################################
# 02 Load and/or install required packages #
############################################

library(Biobase)
library(edgeR)
library(AnnotationFuncs)
library(org.Bt.eg.db)
library(MASS)
library(devtools)
library(plyr)
library(tidyverse)
library(stringr)
library(magrittr)
library(biobroom)
library(ggplot2)
library(ggjoy)
library(ggrepel)
library(VennDiagram)
library(extrafont)
library(grDevices)
library(grid)
library(gridExtra)
library(svglite)
library(tools)

# Uncomment functions below to install packages in case you don't have them

# Bioconductor packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("AnnotationFuncs")
#biocLite("edgeR")
#biocLite("Biobase")
#biocLite("geneLenDataBase")
#biocLite("org.Bt.eg.db")


# CRAN packages
#install.packages("devtools")
#install.packages("plyr")
#install.packages("lubridate")
#install.packages("tidyverse")
#install.packages("stringr")
#install.packages("magrittr")
#install.packages("ggplot2")
#install.packages("MASS")
#install.packages("extrafont")
#install.packages("gdata")
#install.packages("ggrepel")
#install.packages("ggjoy")
#install.packages("gridExtra")
#install.packages("svglite")
#install.packages("VennDiagram")

#########################################
# 03 Import featureCounts sense counts  #
#########################################

# Create a vector of file names
files <- list.files(path        = fileDir,
                    pattern     = "^A6",
                    all.files   = TRUE,
                    full.names  = FALSE,
                    recursive   = FALSE,
                    ignore.case = FALSE)

files

# Create a dataframe with raw counts for all samples
rawCounts <- readDGE(path         = fileDir,
                     files        = files,
                     header       = TRUE,
                     comment.char = "#",
                     columns      = c(1, 7))
names(rawCounts)
head(rawCounts$samples)
head(rawCounts$counts)


#################################
# 04 Clean column and row names #
#################################

# Define pattern to be replaced
str_pattern <- c("_sense-counts", "-")

# Correct column names
colnames(rawCounts$counts) %<>%
  str_replace(str_pattern[1], "") %>% str_replace(str_pattern[2], "_")

# Correct row names
rownames(rawCounts$counts) %<>%
  str_replace("BGD.*,", "") %>%
  str_replace(",miRBase:.*", "") %>%
  str_replace("GeneID:", "")

rownames(rawCounts$samples) %<>%
  str_replace(str_pattern[1], "") %>% str_replace(str_pattern[2], "_")


head(rawCounts$counts)
head(rawCounts$samples)

#############################
# 05 Add sample information #
#############################

# Treatment group (infection)
rawCounts$samples$group <- rownames(rawCounts$samples)
rawCounts$samples$group %<>%
  str_replace("^A.+_W_1_\\w", "pre_infection") %>%
  str_replace("^A.+", "post_infection") %>%
  factor(levels = c("pre_infection", "post_infection"))

# Animal ID
rawCounts$samples$animal <- rownames(rawCounts$samples)
rawCounts$samples$animal %<>%
  str_replace("_W.+_\\w", "") %>%
  factor()

# Time points
rawCounts$samples$time.point <- rownames(rawCounts$samples)
rawCounts$samples$time.point %<>%
  str_replace("A\\d\\d\\d\\d_", "") %>%
  str_replace("_(P|U)", "") %>%
  factor(levels = c("W_1", "W1", "W2", "W6", "W10", "W12"))

# PPD-b stimulation
rawCounts$samples$ppdb_stimulation <- rownames(rawCounts$samples)
rawCounts$samples$ppdb_stimulation %<>%
  str_replace("A\\d\\d\\d\\d_.+_", "") %>%
  str_replace("P", "PPDb_stimulated") %>%
  str_replace("U", "Non_stimulated") %>%
  factor(levels = c("Non_stimulated", "PPDb_stimulated"))


head(rawCounts$samples)

#################################
# 06 Get gene names and symbols #
#################################

# Create annotation table with counts information
annotCounts <- as.data.frame(rawCounts$counts)
columns(org.Bt.eg.db)

# Retrieve gene names using NCBI Ref Seq gene identifiers
annotCounts$gene_name <- mapIds(org.Bt.eg.db,
                                keys      = rownames(annotCounts),
                                column    = "GENENAME",
                                keytype   = "ENTREZID",
                                multiVals = "first")

# Retrieve gene symbols using NCBI Ref Seq gene identifiers
annotCounts$gene_symbol <- mapIds(org.Bt.eg.db,
                                  keys      = rownames(annotCounts),
                                  column    = "SYMBOL",
                                  keytype   = "ENTREZID",
                                  multiVals = "first")

head(annotCounts)
dim(annotCounts)

# Ouptut data
write.table(x         = annotCounts,
            file      = "PPDb-RNA-seq_sense_annot-rawcounts.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

#####################
# 07 Create DGElist #
#####################

# Assign required information to variables
gene_annotation <- dplyr::select(annotCounts,
                                 gene_name,
                                 gene_symbol)

raw_counts <- as.data.frame(rawCounts$counts)

group <- rawCounts$samples$group

# Use newly assigned variables to create DGElist
PPDb_dgelist <- DGEList(counts       = raw_counts,
                        group        = group,
                        genes        = gene_annotation,
                        lib.size     = NULL,
                        norm.factors = NULL,
                        remove.zeros = FALSE)

names(PPDb_dgelist)
dim(PPDb_dgelist)
head(PPDb_dgelist$counts)
head(PPDb_dgelist$samples)
head(PPDb_dgelist$genes)

# Include addtional experimental information into DGElist
identical(rownames(rawCounts$samples), rownames(PPDb_dgelist$samples))

PPDb_dgelist$samples$animal <- rawCounts$samples$animal
PPDb_dgelist$samples$time.point <- rawCounts$samples$time.point
PPDb_dgelist$samples$ppdb_stimulation <- rawCounts$samples$ppdb_stimulation

head(PPDb_dgelist$samples)


################################################
# 08 Density plot: raw gene counts per library #
################################################

# Tidy DGElist and plot data
PPDb_dgelist %>%
  tidy() %>%
  ggplot() +
      geom_density(aes(x     = log10(count + 1),
                       group = sample)) +
      theme_bw() +
      ylab("Density of raw gene counts per sample") +
      xlab(expression(paste(log[10], "(counts + 1)"))) -> density_raw


density_raw

# Export high quality image
ggsave("PPDb-density_plot_raw_counts.png",
       plot      = density_raw,
       device    = "png",
       limitsize = FALSE,
       dpi       = 600,
       path      = workDir)

###########################################
# 09 Remove zero and lowly expressed tags #
###########################################

# Filter non expressed tags (all genes that have zero counts in all samples)
PPDb_no_zeros <- PPDb_dgelist[rowSums(PPDb_dgelist$counts) > 0, ]
dim(PPDb_no_zeros$counts)
head(PPDb_no_zeros$counts)
colnames(PPDb_no_zeros$counts)

# Filter lowly expressed tags, retaining only tags with
# more than 1 count per million in 10 or more libraries
# (10 libraries correspond to 10 biological replicates and represent
# one PPDb-stimulated [P] or Non-stimulated [U] group at any given time point)
PPDb_filt <- PPDb_no_zeros[rowSums(cpm(PPDb_no_zeros) > 1) >= 10, ]
dim(PPDb_filt$counts)
head(PPDb_filt$counts)

# Output file of filtered counts
write.table(x         = PPDb_filt$counts,
            file      = "PPDb_filt_counts.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)


##############################
# 10 Recompute library sizes #
##############################

PPDb_filt$samples$lib.size <- colSums(PPDb_filt$counts)
head(PPDb_filt$samples)
head(PPDb_dgelist$samples)


###########################################################################
# 11 Calculate normalisation factors using Trimmed Mean of M-values (TMM) #
###########################################################################

# With edgeR, counts are not transformed in any way after normalisation
PPDb_norm <- calcNormFactors(PPDb_filt, method = "TMM")
head(PPDb_norm$samples)

##########################################################
# 12 Tidy DGElist for exploratory data analysis plotting #
##########################################################

PPDb_norm %>%
  tidy() %>%
  dplyr::mutate(ppdb_stimulation = sample, labels = sample) -> tidy_PPDbNorm

# Add PPDb_stimulation info
tidy_PPDbNorm$ppdb_stimulation %<>%
  stringr::str_replace("A\\d\\d\\d\\d_.+_", "") %>%
  stringr::str_replace("P", "PPDb_stimulated") %>%
  stringr::str_replace("U", "Non_stimulated") %>%
  factor(levels = c("Non_stimulated", "PPDb_stimulated"))

# Add plotting labels
tidy_PPDbNorm$labels %<>%
  stringr::str_replace("A", "") %>%
  stringr::str_replace("_(P|U)", "") %>%
  stringr::str_replace("W_1", "W-1") %>%
  factor(levels = c("6511_W-1", "6511_W1", "6511_W2", "6511_W10",
                    "6514_W-1", "6514_W1", "6514_W2", "6514_W10",
                    "6520_W-1", "6520_W1", "6520_W2", "6520_W10",
                    "6522_W-1", "6522_W1", "6522_W2", "6522_W10",
                    "6526_W-1", "6526_W1", "6526_W2", "6526_W10",
                    "6635_W-1", "6635_W1", "6635_W2", "6635_W10",
                    "6636_W-1", "6636_W1", "6636_W2", "6636_W10",
                    "6637_W-1", "6637_W1", "6637_W2", "6637_W10",
                    "6644_W-1", "6644_W1", "6644_W2", "6644_W10",
                    "6698_W-1", "6698_W1", "6698_W2", "6698_W10"))

###########################################################
# 13 Joyplot: density of filtered gene counts per library #
###########################################################

ggplot(tidy_PPDbNorm, aes(x = log10(count + 1),
                          y = labels)) +
    scale_y_discrete(limits = rev(levels(tidy_PPDbNorm$labels))) +
    geom_joy(aes(fill = ppdb_stimulation), alpha = 0.5) +
    scale_fill_manual("Treatment", values = c("#91bfdb", "#fc8d59")) +
    theme_bw() +
    facet_grid(. ~ ppdb_stimulation) +
    ylab("Density of gene counts (CPM > 1) per sample") +
    xlab(expression(paste(log[10], "(counts + 1)"))) -> density_norm


density_norm

# Export high quality image
ggsave("PPDb-density_plot_filt_counts.svg",
       plot      = density_norm,
       device    = "svg",
       limitsize = FALSE,
       dpi       = 600,
       path      = workDir)

#########################################
# 14 MDS #
#########################################



# Plot MDS of all samples
test<-plotMDS.DGEList(x = PPDb_norm, method = "bcv", plot = FALSE)
ggplot(as.data.frame(test), aes(x = x, y = y)) +
  geom_point(aes(colour = names(test$x)))


png(filename = "MDS_all_samples.png",
    width    = 1366,
    height   = 768,
    units    = "px")

BiocGenerics::plotPCA(t(PPDb_norm))

dev.off()

# Plot MDS of Week -1 time point
png(filename = "MDS_W-1.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x = PPDb_norm[, grep(pattern = "_W.1_", x = colnames(PPDb_norm))])

dev.off()

# Plot MDS of Week +1 time point
png(filename = "MDS_W1.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x = PPDb_norm[, grep(pattern = "_W1_", x = colnames(PPDb_norm))])

dev.off()

# Plot MDS of Week +2 time point
png(filename = "MDS_W2.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x = PPDb_norm[, grep(pattern = "_W2_", x = colnames(PPDb_norm))])

dev.off()

# Plot MDS of Week +10 time point
png(filename = "MDS_W10.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x = PPDb_norm[, grep(pattern = "_W10_", x = colnames(PPDb_norm))])

dev.off()

###############################
# Define experimental factors # ----
###############################

head(PPDb_norm$samples)

condition <- factor(PPDb_norm$samples$group,
                    levels = c("U", "P"))

time.point <- factor(PPDb_norm$samples$time.point,
                     levels = c("W_1", "W1", "W2", "W10"))

animal <- factor(PPDb_norm$samples$animal)

batch <- factor(PPDb_norm$samples$batch)

cond.time <- factor(paste(PPDb_norm$samples$group,
                          PPDb_norm$samples$time.point,
                          sep="."),
                          levels = c("U.W_1", "U.W1", "U.W2", "U.W10",
                                     "P.W_1", "P.W1", "P.W2", "P.W10"))

##############################################
# Create a design matrix for paired analysis # ----
##############################################

# Create a design matrix with animal as a blocking factor
block_animal <- model.matrix(~animal + cond.time,
                             data = PPDb_norm$samples)


dim(block_animal)
dim(PPDb_norm$samples)
head(block_animal)

#colnames(block_animal) <- gsub(pattern = "(time.point)|(condition)",
#                               replacement = "",
#                               x = colnames(block_animal),
#                               perl = TRUE)
#block_animal


# Output the design matrix info
write.table(x         = block_animal,
            file      = "PPDb_design-matrix-animal-block.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

############################################################################
# Estimate the dispersion parameter for each tag using the Cox-Reid method # ----
#                      (for multi-factor data)                             #
############################################################################

PPDb_disp <- estimateGLMCommonDisp(y       = PPDb_norm,
                                   design  = block_animal,
                                   verbose = TRUE)

PPDb_disp <- estimateGLMTrendedDisp(y      = PPDb_disp,
                                    design = block_animal)

PPDb_disp <- estimateGLMTagwiseDisp(y      = PPDb_disp,
                                    design = block_animal)

names(PPDb_disp)

# Plot the dispersion
png(filename = "BCV_PPDb.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotBCV(PPDb_disp)

dev.off()

# Show the calculated dispersion
PPDb_disp$common.dispersion

# And show its square root, the coefficient of biological variation
sqrt(PPDb_disp$common.dispersion)
sqrt(PPDb_disp$tagwise.dispersion)
# Create a matrix of the tagwise dispersion associated with gene information
Tagwisedisp <- cbind(PPDb_disp$genes, PPDb_disp$tagwise.dispersion)
head(Tagwisedisp)
dim(Tagwisedisp)

# Write into a table the calculated tagwise dispersion
write.matrix(x    = Tagwisedisp,
             file = "PPDb_Tagwise_dispersion.txt",
             sep  = "\t")

##################################################################
# Determine differential expression using negative binomial GLMs # ----
##################################################################

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
PPDb_fit <- glmFit(y = PPDb_disp, design = block_animal)
names(PPDb_fit)
colnames(PPDb_fit$design)

# Test for differential expression between the different time points/treatments,
# using the coefficients from PPDbfit$design
pre1.lrt <- glmLRT(PPDb_fit, coef = "cond.timeP.W_1")
DE.pre1 <- topTags(object        = pre1.lrt,
                   n             = "inf",
                   adjust.method = "BH")
head(DE.pre1$table)

W1.lrt <- glmLRT(PPDb_fit, coef = "cond.timeP.W1")
DE.W1 <- topTags(object        = W1.lrt,
                 n             = "inf",
                 adjust.method = "BH")
head(DE.W1$table)

W2.lrt <- glmLRT(PPDb_fit, coef = "cond.timeP.W2")
DE.W2 <- topTags(object        = W2.lrt,
                 n             = "inf",
                 adjust.method = "BH")
head(DE.W2$table)

W10.lrt <- glmLRT(PPDb_fit, coef = "cond.timeP.W10")
DE.W10 <- topTags(object        = W10.lrt,
                  n             = "inf",
                  adjust.method = "BH")
head(DE.W10$table)

### Previous test with DE call when specifying contrasts (no intercept in design matrix)
#my.contrasts <- makeContrasts(pre1   = P.W_1-U.W_1,
#                              W1     = P.W1-U.W1,
#                              W2     = P.W2-U.W2,
#                              W10    = P.W10-U.W10,
#                             levels = PPDb_fit)
#pre1.lrt <- glmLRT(PPDb_fit,
#                   contrast = my.contrasts[,"pre1"])
#DE.pre1 <- topTags(object        = pre1.lrt,
#                   n             = "inf",
#                   adjust.method = "BH")
#head(DE.pre1$table)

#W1.lrt <- glmLRT(PPDb_fit,
#                 contrast = my.contrasts[,"W1"])
#DE.W1 <- topTags(object        = W1.lrt,
#                 n             = "inf",
#                 adjust.method = "BH")
#head(DE.W1$table)

#W2.lrt <- glmLRT(PPDb_fit,
#                 contrast = my.contrasts[,"W2"])
#DE.W2 <- topTags(object        = W2.lrt,
#                 n             = "inf",
#                 adjust.method = "BH")
#head(DE.W2$table)

#W10.lrt <- glmLRT(PPDb_fit,
#                  contrast = my.contrasts[,"W10"])
#DE.W10 <- topTags(object        = W10.lrt,
#                 n             = "inf",
#                  adjust.method = "BH")
#head(DE.W10$table)

#########################################
# Output lists of DE genes: FDR < 0.05  # ----
# Comparisons P vs U at each time point #
#########################################

FDR_0.05_DE.pre1 <- subset(DE.pre1$table, FDR < 0.05)
write.table(x         = FDR_0.05_DE.pre1,
            file      = "FDR_0-05_pre1.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_DE.W1 <- subset(DE.W1$table, FDR < 0.05)
write.table(x         = FDR_0.05_DE.W1,
            file      = "FDR_0-05_W1.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_DE.W2 <- subset(DE.W2$table, FDR < 0.05)
write.table(x         = FDR_0.05_DE.W2,
            file      = "FDR_0-05_W2.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_DE.W10 <- subset(DE.W10$table, FDR < 0.05)
write.table(x         = FDR_0.05_DE.W10,
            file      = "FDR_0-05_W10.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

###############################################################
# Output lists of DE genes: FDR < 0.05 and absolute logFC > 1 # ----
###############################################################

FDRlogFC_DE.pre1 <- subset(DE.pre1$table, abs(logFC) > 1 & FDR < 0.05)
write.table(x         = FDRlogFC_DE.pre1,
            file      = "FDRlogFC_pre1.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDRlogFC_DE.W1 <- subset(DE.W1$table, abs(logFC) > 1 & FDR < 0.05)
write.table(x         = FDRlogFC_DE.W1,
            file      = "FDRlogFC_W1.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDRlogFC_DE.W2 <- subset(DE.W2$table, abs(logFC) > 1 & FDR < 0.05)
write.table(x         = FDRlogFC_DE.W2,
            file      = "FDRlogFC_W2.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDRlogFC_DE.W10 <- subset(DE.W10$table, abs(logFC) > 1 & FDR < 0.05)
write.table(x         = FDRlogFC_DE.W10,
            file      = "FDRlogFC_W10.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

###################################################
# Output lists of all DE genes at each time point # ----
#  [without selection by FDR or absolute logFC]   #
###################################################

write.table(x         = DE.pre1,
            file      = "DE_pre1.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE.W1,
            file      = "DE_W1.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE.W2,
            file      = "DE_W2.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE.W10,
            file      = "DE_W10.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

#########################################################
# Merge all DE call data from the different time points # ----
#              into a single data frame                 #
#     [without selection by FDR or absolute logFC]      #
#########################################################

Full_DE_PPDb <- merge(x  = DE.pre1$table,
                      y  = DE.W1$table
                      [,(ncol(DE.W1$table) - 4) :
                          ncol(DE.W1$table)],
                      by = "row.names")
head(Full_DE_PPDb)

rownames(Full_DE_PPDb) <- Full_DE_PPDb[, 1] # Correct row names before
Full_DE_PPDb <- Full_DE_PPDb[, -1]          # merging the other data frames.

colnames(Full_DE_PPDb) <- gsub(pattern     = ".x$",
                               replacement = "_pre1",
                               x           = colnames(Full_DE_PPDb),
                               perl        = TRUE)
colnames(Full_DE_PPDb) <- gsub(pattern     = ".y$",
                               replacement = "_W1",
                               x           = colnames(Full_DE_PPDb),
                               perl        = TRUE)
head(Full_DE_PPDb)

Full_DE_PPDb <- merge(x  = Full_DE_PPDb,
                      y  = DE.W2$table
                      [, (ncol(DE.W2$table) - 4) :
                           ncol(DE.W2$table)],
                      by = "row.names")
head(Full_DE_PPDb)

rownames(Full_DE_PPDb) <- Full_DE_PPDb[, 1]
Full_DE_PPDb <- Full_DE_PPDb[, -1]

Full_DE_PPDb <- merge(x  = Full_DE_PPDb,
                      y  = DE.W10$table
                      [, (ncol(DE.W10$table) - 4) :
                           ncol(DE.W10$table)],
                      by = "row.names")
head(Full_DE_PPDb)

rownames(Full_DE_PPDb) <- Full_DE_PPDb[, 1]
Full_DE_PPDb <- Full_DE_PPDb[, -1]

colnames(Full_DE_PPDb) <- gsub(pattern     = ".x$",
                               replacement = "_W2",
                               x           = colnames(Full_DE_PPDb),
                               perl        = TRUE)
colnames(Full_DE_PPDb) <- gsub(pattern     = ".y$",
                               replacement = "_W10",
                               x           = colnames(Full_DE_PPDb),
                               perl        = TRUE)
head(Full_DE_PPDb)

write.table(x         = Full_DE_PPDb,
            file      = "Full_DE_PPDb.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

########################################
# Output lists of DE genes: FDR < 0.01 # ----
########################################

FDR_0.01_DE.pre1 <- subset(DE.pre1$table, FDR < 0.01)
write.table(x         = FDR_0.01_DE.pre1,
            file      = "FDR_0-01_pre1.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.01_DE.W1 <- subset(DE.W1$table, FDR < 0.01)
write.table(x         = FDR_0.01_DE.W1,
            file      = "FDR_0-01_W1.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.01_DE.W2 <- subset(DE.W2$table, FDR < 0.01)
write.table(x         = FDR_0.01_DE.W2,
            file      = "FDR_0-01_W2.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.01_DE.W10 <- subset(DE.W10$table, FDR < 0.01)
write.table(x         = FDR_0.01_DE.W10,
            file      = "FDR_0-01_W10.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

#########################################
# Output lists of DE genes: FDR < 0.001 # ----
#########################################

FDR_0001_DE.pre1 <- subset(DE.pre1$table, FDR < 0.001)
write.table(x         = FDR_0001_DE.pre1,
            file      = "FDR_0001_pre1.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0001_DE.W1 <- subset(DE.W1$table, FDR < 0.001)
write.table(x         = FDR_0001_DE.W1,
            file      = "FDR_0001_W1.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0001_DE.W2 <- subset(DE.W2$table, FDR < 0.001)
write.table(x         = FDR_0001_DE.W2,
            file      = "FDR_0001_W2.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0001_DE.W10 <- subset(DE.W10$table, FDR < 0.001)
write.table(x         = FDR_0001_DE.W10,
            file      = "FDR_0001_W10.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

########################################
# Venn diagram of DE genes: FDR < 0.05 # ----
########################################

# Turn gene IDs of significant DE genes (FDR < 0.05)
# per time point into vectors
Pre1_0.05.vector <- c(rownames(FDR_0.05_DE.pre1))
W1_0.05.vector <- c(rownames(FDR_0.05_DE.W1))
W2_0.05.vector <- c(rownames(FDR_0.05_DE.W2))
W10_0.05.vector <- c(rownames(FDR_0.05_DE.W10))

# Turn log files off from VennDiagram package
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Plot Venn diagram for all time points
venn.plot <- venn.diagram(list(A = Pre1_0.05.vector,
                               B = W10_0.05.vector,
                               C = W1_0.05.vector,
                               D = W2_0.05.vector),
                          filename        = "Venn_DE_FDR_0-05.svg",
                          imagetype       = "svg",
                          col             = "transparent",
                          fill            = c("#ffffcc",
                                              "#225ea8",
                                              "#a1dab4",
                                              "#41b6c4"),
                          alpha           = 0.50,
                          label.col       = "#003333",
                          cex             = 10,
                          fontfamily      = "Raavi",
                          category.names  = c("-1 wk",
                                              "+10 wk",
                                              "+1 wk",
                                              "+2 wk"),
                          cat.col         = "black",
                          cat.cex         = 10,
                          cat.pos         = c(-11, 11, 0, 0),
                          cat.dist        = c(0.21, 0.21, 0.1, 0.1),
                          cat.fontfamily  = "Raavi",
                          rotation.degree = 360,
                          margin          = 0,
                          height          = 85,
                          width           = 85,
                          units           = 'cm',
                          compression     = 'lzw',
                          resolution      = 300)

###############################################################
# Venn diagram of DE genes: FDR < 0.05 and absolute logFC > 1 # ----
###############################################################

# Turn gene IDs of significant DE genes (FDR < 0.05 and absolute logFC > 1)
# per time point into vectors
Pre1.vector <- c(rownames(FDRlogFC_DE.pre1))
W1.vector <- c(rownames(FDRlogFC_DE.W1))
W2.vector <- c(rownames(FDRlogFC_DE.W2))
W10.vector <- c(rownames(FDRlogFC_DE.W10))

# Turn log files off from VennDiagram package
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Plot Venn diagram for all time points
venn.plot <- venn.diagram(list(A = Pre1.vector,
                               B = W10.vector,
                               C = W1.vector,
                               D = W2.vector),
                          filename        = "Venn_DE_FDRlogFC.svg",
                          imagetype       = "svg",
                          col             = "transparent",
                          fill            = c("#ffffcc",
                                              "#225ea8",
                                              "#a1dab4",
                                              "#41b6c4"),
                          alpha           = 0.50,
                          label.col       = "#003333",
                          cex             = 10,
                          fontfamily      = "Raavi",
                          category.names  = c("-1 wk",
                                              "+10 wk",
                                              "+1 wk",
                                              "+2 wk"),
                          cat.col         = "black",
                          cat.cex         = 10,
                          cat.pos         = c(-11, 11, 0, 0),
                          cat.dist        = c(0.21, 0.21, 0.1, 0.1),
                          cat.fontfamily  = "Raavi",
                          rotation.degree = 360,
                          margin          = 0,
                          height          = 85,
                          width           = 85,
                          units           = 'cm',
                          compression     = 'lzw',
                          resolution      = 300)

########################################
# Venn diagram of DE genes: FDR < 0.01 # ----
########################################

# Turn gene IDs of significant DE genes (FDR < 0.01)
# per time point into vectors
Pre1_0.01.vector <- c(rownames(FDR_0.01_DE.pre1))
W1_0.01.vector <- c(rownames(FDR_0.01_DE.W1))
W2_0.01.vector <- c(rownames(FDR_0.01_DE.W2))
W10_0.01.vector <- c(rownames(FDR_0.01_DE.W10))

# Turn log files off from VennDiagram package
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Plot Venn diagram for all time points
venn.plot <- venn.diagram(list(A = Pre1_0.01.vector,
                               B = W10_0.01.vector,
                               C = W1_0.01.vector,
                               D = W2_0.01.vector),
                          filename        = "Venn_DE_FDR_0-01.svg",
                          imagetype       = "svg",
                          col             = "transparent",
                          fill            = c("#ffffcc",
                                              "#225ea8",
                                              "#a1dab4",
                                              "#41b6c4"),
                          alpha           = 0.50,
                          label.col       = "#003333",
                          cex             = 10,
                          fontfamily      = "Raavi",
                          category.names  = c("-1 wk",
                                              "+10 wk",
                                              "+1 wk",
                                              "+2 wk"),
                          cat.col         = "black",
                          cat.cex         = 10,
                          cat.pos         = c(-11, 11, 0, 0),
                          cat.dist        = c(0.21, 0.21, 0.1, 0.1),
                          cat.fontfamily  = "Raavi",
                          rotation.degree = 360,
                          margin          = 0,
                          height          = 85,
                          width           = 85,
                          units           = 'cm',
                          compression     = 'lzw',
                          resolution      = 300)

################################################################
# Volcano plots of DE genes: FDR < 0.05 and absolute logFC > 1 # ----
################################################################

# Pre-infection -1wk time point
dfPre1 <- as.data.frame(DE.pre1,
                        row.names        = rownames(DE.pre1$table),
                        cut.names        = FALSE,
                        col.names        = colnames(DE.pre1$table),
                        fix.empty.names  = FALSE,
                        stringsAsFactors = FALSE)
dfPre1 <- transform(dfPre1,
                    threshold = as.factor(abs(dfPre1$logFC) > 1
                                          & dfPre1$FDR < 0.05))

volcanoPre1 <- ggplot(data       = dfPre1,
                      aes(x      = logFC,
                          y      = -log10(FDR),
                          colour = threshold)) +
    geom_point(alpha = 0.4,
               size  = 1.75) +
    xlim(c(-6, 6)) +
    ylim(c(-1, 20)) +
    theme(legend.position   = "right",
          legend.background = element_rect(colour = "black"),
          text              = element_text(size   = 16,
                                           family = "Raavi")) +
    ggtitle("-1 wk") +
    xlab(expression(paste(log[2], "-fold change"))) +
    ylab(expression(paste(-log[10], " adjusted p-value"))) +
    scale_colour_manual("Passed \ncut-off",
                        labels = c("False", "True"),
                        values = c('#41b6c4', '#225ea8'))
volcanoPre1

ggsave("volcano_pre1.svg",
       plot      = volcanoPre1,
       limitsize = FALSE,
       width     = 20,
       height    = 20,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

# Post-infection +1wk time point
dfW1 <- as.data.frame(DE.W1,
                      row.names        = rownames(DE.W1$table),
                      cut.names        = FALSE,
                      col.names        = colnames(DE.W1$table),
                      fix.empty.names  = FALSE,
                      stringsAsFactors = FALSE)
dfW1 <- transform(dfW1,
                  threshold = as.factor(abs(dfW1$logFC) > 1
                                        & dfW1$FDR < 0.05))

volcanoW1 <- ggplot(data       = dfW1,
                    aes(x      = logFC,
                        y      = -log10(FDR),
                        colour = threshold)) +
    geom_point(alpha = 0.4,
               size  = 1.75) +
    xlim(c(-6, 7)) +
    ylim(c(-1, 19)) +
    theme(legend.position   = "right",
          legend.background = element_rect(colour = "black"),
          text              = element_text(size   = 16,
                                           family = "Raavi")) +
    ggtitle("+1 wk") +
    xlab(expression(paste(log[2], "-fold change"))) +
    ylab(expression(paste(-log[10], " adjusted p-value"))) +
    scale_colour_manual("Passed \ncut-off",
                        labels = c("False", "True"),
                        values = c('#41b6c4', '#225ea8'))
volcanoW1

ggsave("volcano_W1.svg",
       plot      = volcanoW1,
       limitsize = FALSE,
       width     = 18,
       height    = 18,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

# Post-infection +2wk time point
dfW2 <- as.data.frame(DE.W2,
                      row.names        = rownames(DE.W2$table),
                      cut.names        = FALSE,
                      col.names        = colnames(DE.W2$table),
                      fix.empty.names  = FALSE,
                      stringsAsFactors = FALSE)
dfW2 <- transform(dfW2,
                  threshold = as.factor(abs(dfW2$logFC) > 1
                                        & dfW2$FDR < 0.05))

volcanoW2 <- ggplot(data       = dfW2,
                    aes(x      = logFC,
                        y      = -log10(FDR),
                        colour = threshold)) +
    geom_point(alpha = 0.4,
               size  = 1.75) +
    xlim(c(-9, 8)) +
    ylim(c(-1, 27)) +
    theme(legend.position   = "right",
          legend.background = element_rect(colour = "black"),
          text              = element_text(size   = 16,
                                           family = "Raavi")) +
    ggtitle("+2 wk") +
    xlab(expression(paste(log[2], "-fold change"))) +
    ylab(expression(paste(-log[10], " adjusted p-value"))) +
    scale_colour_manual("Passed \ncut-off",
                        labels = c("False", "True"),
                        values = c('#41b6c4', '#225ea8'))
volcanoW2

ggsave("volcano_W2.svg",
       plot      = volcanoW2,
       limitsize = FALSE,
       width     = 18,
       height    = 18,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

# Post-infection +10wk time point
dfW10 <- as.data.frame(DE.W10,
                       row.names        = rownames(DE.W10$table),
                       cut.names        = FALSE,
                       col.names        = colnames(DE.W10$table),
                       fix.empty.names  = FALSE,
                       stringsAsFactors = FALSE)
dfW10 <- transform(dfW10,
                   threshold = as.factor(abs(dfW10$logFC) > 1
                                         & dfW10$FDR < 0.05))

volcanoW10 <- ggplot(data       = dfW10,
                     aes(x      = logFC,
                         y      = -log10(FDR),
                         colour = threshold)) +
    geom_point(alpha = 0.4,
               size  = 1.75) +
    xlim(c(-10, 11)) +
    ylim(c(-1, 65)) +
    theme(legend.position   = "right",
          legend.background = element_rect(colour = "black"),
          text              = element_text(size   = 16,
                                           family = "Raavi")) +
    ggtitle("+10 wk") +
    xlab(expression(paste(log[2], "-fold change"))) +
    ylab(expression(paste(-log[10], " adjusted p-value"))) +
    scale_colour_manual("Passed \ncut-off",
                        labels = c("False", "True"),
                        values = c('#41b6c4', '#225ea8'))
volcanoW10

ggsave("volcano_W10.svg",
       plot      = volcanoW10,
       limitsize = FALSE,
       width     = 20,
       height    = 22,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

#############################################
# Combine all volcano plots into one figure # ----
#############################################

all_volcano <- arrangeGrob(volcanoPre1,
                           volcanoW1,
                           volcanoW2,
                           volcanoW10,
                           ncol = 2)
arrange_all_volcano <- grid.arrange(all_volcano, padding = unit(1, "line"))
ggsave("all_volcano_plots.png",
       plot      = arrange_all_volcano,
       width     = 70,
       height    = 70,
       dpi       = 300,
       path      = workDir,
       units     = "cm",
       limitsize = FALSE)

#########################################
# Volcano plots of DE genes: FDR < 0.05 # ----
#########################################

# Pre-infection -1wk time point
dfPre1_0.05 <- transform(dfPre1,
                         threshold = as.factor(dfPre1$FDR < 0.05))

volcanoPre1_0.05 <- ggplot(data       = dfPre1_0.05,
                           aes(x      = logFC,
                               y      = -log10(FDR),
                               colour = threshold)) +
                        geom_point(alpha = 0.4, size  = 1.75) +
                        geom_text_repel(data = subset(dfPre1_0.05,
                                                      abs(logFC) > 3
                                                      & FDR < 0.001),
                                        aes(label = gene.symbol),
                                        fontface = "bold",
                                        colour = "#404040",
                                        size = 4,
                                        box.padding = unit(0.35, "lines"),
                                        point.padding = unit(0.3, "lines")) +
                        xlim(c(-6, 6)) +
                        ylim(c(-1, 20)) +
                        theme(legend.position   = "right",
                              legend.background = element_rect(colour = "black"),
                              axis.title        = element_text(size   = 16),
                              text              = element_text(size   = 14,
                                                               family = "Raavi")) +
                        ggtitle("-1 wk") +
                        xlab(expression(paste(log[2], "-fold change"))) +
                        ylab(expression(paste(-log[10], " FDR"))) +
                        scale_colour_manual("FDR \n< 0.05",
                                            labels = c("False", "True"),
                                            values = c('#41b6c4', '#225ea8'))
volcanoPre1_0.05

ggsave("labeled_volcano_pre1_0.05.svg",
       plot      = volcanoPre1_0.05,
       limitsize = FALSE,
       width     = 20,
       height    = 20,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

# Post-infection +1wk time point
dfW1_0.05 <- transform(dfW1,
                       threshold = as.factor(dfW1$FDR < 0.05))

volcanoW1_0.05 <- ggplot(data       = dfW1_0.05,
                         aes(x      = logFC,
                             y      = -log10(FDR),
                             colour = threshold)) +
                    geom_point(alpha = 0.4, size  = 1.75) +
                    geom_text_repel(data = subset(dfW1_0.05,
                                                  abs(logFC) > 3
                                                  & FDR < 0.001),
                    aes(label = gene.symbol),
                    fontface = "bold",
                    colour = "#404040",
                    size = 4,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines")) +
                    xlim(c(-6, 7)) +
                    ylim(c(-1, 19)) +
                    theme(legend.position   = "right",
                          legend.background = element_rect(colour = "black"),
                          axis.title        = element_text(size   = 16),
                          text              = element_text(size   = 14,
                                                           family = "Raavi")) +
                    ggtitle("+1 wk") +
                    xlab(expression(paste(log[2], "-fold change"))) +
                    ylab(expression(paste(-log[10], " FDR"))) +
                    scale_colour_manual("FDR \n< 0.05",
                                        labels = c("False", "True"),
                                        values = c('#41b6c4', '#225ea8'))
volcanoW1_0.05

ggsave("labeled_volcano_W1_0-05.svg",
       plot      = volcanoW1_0.05,
       limitsize = FALSE,
       width     = 18,
       height    = 18,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

# Post-infection +2wk time point
dfW2_0.05 <- transform(dfW2,
                       threshold = as.factor(dfW2$FDR < 0.05))

volcanoW2_0.05 <- ggplot(data       = dfW2_0.05,
                         aes(x      = logFC,
                             y      = -log10(FDR),
                             colour = threshold)) +
                    geom_point(alpha = 0.4, size  = 1.75) +
                    geom_text_repel(data = subset(dfW2_0.05,
                                                  abs(logFC) > 4
                                                  & FDR < 0.001),
                    aes(label = gene.symbol),
                    fontface = "bold",
                    colour = "#404040",
                    size = 4,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines")) +
                    xlim(c(-9, 8)) +
                    ylim(c(-1, 27)) +
                    theme(legend.position   = "right",
                          legend.background = element_rect(colour = "black"),
                          axis.title        = element_text(size   = 16),
                          text              = element_text(size   = 14,
                                                           family = "Raavi")) +
                    ggtitle("+2 wk") +
                    xlab(expression(paste(log[2], "-fold change"))) +
                    ylab(expression(paste(-log[10], " FDR"))) +
                    scale_colour_manual("FDR \n< 0.05",
                                        labels = c("False", "True"),
                                        values = c('#41b6c4', '#225ea8'))
volcanoW2_0.05

ggsave("labeled_volcano_W2_0-05.svg",
       plot      = volcanoW2_0.05,
       limitsize = FALSE,
       width     = 18,
       height    = 18,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

# Post-infection +10wk time point
dfW10_0.05 <- transform(dfW10,
                        threshold = as.factor(dfW10$FDR < 0.05))

volcanoW10_0.05 <- ggplot(data       = dfW10_0.05,
                          aes(x      = logFC,
                              y      = -log10(FDR),
                              colour = threshold)) +
                    geom_point(alpha = 0.4, size  = 1.75) +
                    geom_text_repel(data = subset(dfW10_0.05,
                                                  abs(logFC) > 5
                                                  & FDR < 0.001),
                    aes(label = gene.symbol),
                    fontface = "bold",
                    colour = "#404040",
                    size = 4,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines")) +
                    xlim(c(-10, 11)) +
                    ylim(c(-1, 65)) +
                    theme(legend.position   = "right",
                          legend.background = element_rect(colour = "black"),
                          axis.title        = element_text(size   = 16),
                          text              = element_text(size   = 14,
                                                           family = "Raavi")) +
                    ggtitle("+10 wk") +
                    xlab(expression(paste(log[2], "-fold change"))) +
                    ylab(expression(paste(-log[10], " FDR"))) +
                    scale_colour_manual("FDR \n< 0.05",
                                        labels = c("False", "True"),
                                        values = c('#41b6c4', '#225ea8'))
volcanoW10_0.05

ggsave("labeled_volcano_W10_0-05.svg",
       plot      = volcanoW10_0.05,
       limitsize = FALSE,
       width     = 20,
       height    = 22,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

#########################################
# Volcano plots of DE genes: FDR < 0.01 # ----
#########################################

# Pre-infection -1wk time point
dfPre1_0.01 <- transform(dfPre1,
                         threshold = as.factor(dfPre1$FDR < 0.01))

volcanoPre1_0.01 <- ggplot(data       = dfPre1_0.01,
                           aes(x      = logFC,
                               y      = -log10(FDR),
                               colour = threshold)) +
    geom_point(alpha = 0.4,
               size  = 1.75) +
    xlim(c(-6, 6)) +
    ylim(c(-1, 20)) +
    theme(legend.position   = "right",
          legend.background = element_rect(colour = "black"),
          axis.title        = element_text(size   = 16),
          text              = element_text(size   = 14,
                                           family = "Raavi")) +
    ggtitle("-1 wk") +
    xlab(expression(paste(log[2], "-fold change"))) +
    ylab(expression(paste(-log[10], " adjusted p-value"))) +
    scale_colour_manual("FDR \n< 0.01",
                        labels = c("False", "True"),
                        values = c('#41b6c4', '#225ea8'))
volcanoPre1_0.01

ggsave("volcano_pre1_0.01.svg",
       plot      = volcanoPre1_0.01,
       limitsize = FALSE,
       width     = 20,
       height    = 20,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

# Post-infection +1wk time point
dfW1_0.01 <- transform(dfW1,
                       threshold = as.factor(dfW1$FDR < 0.01))

volcanoW1_0.01 <- ggplot(data       = dfW1_0.01,
                         aes(x      = logFC,
                             y      = -log10(FDR),
                             colour = threshold)) +
    geom_point(alpha = 0.4,
               size  = 1.75) +
    xlim(c(-6, 7)) +
    ylim(c(-1, 19)) +
    theme(legend.position   = "right",
          legend.background = element_rect(colour = "black"),
          axis.title        = element_text(size   = 16),
          text              = element_text(size   = 14,
                                           family = "Raavi")) +
    ggtitle("+1 wk") +
    xlab(expression(paste(log[2], "-fold change"))) +
    ylab(expression(paste(-log[10], " adjusted p-value"))) +
    scale_colour_manual("FDR \n< 0.01",
                        labels = c("False", "True"),
                        values = c('#41b6c4', '#225ea8'))
volcanoW1_0.01

ggsave("volcano_W1_0-01.svg",
       plot      = volcanoW1_0.01,
       limitsize = FALSE,
       width     = 18,
       height    = 18,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

# Post-infection +2wk time point
dfW2_0.01 <- transform(dfW2,
                       threshold = as.factor(dfW2$FDR < 0.01))

volcanoW2_0.01 <- ggplot(data       = dfW2_0.01,
                         aes(x      = logFC,
                             y      = -log10(FDR),
                             colour = threshold)) +
    geom_point(alpha = 0.4,
               size  = 1.75) +
    xlim(c(-9, 8)) +
    ylim(c(-1, 27)) +
    theme(legend.position   = "right",
          legend.background = element_rect(colour = "black"),
          axis.title        = element_text(size   = 16),
          text              = element_text(size   = 14,
                                           family = "Raavi")) +
    ggtitle("+2 wk") +
    xlab(expression(paste(log[2], "-fold change"))) +
    ylab(expression(paste(-log[10], " adjusted p-value"))) +
    scale_colour_manual("FDR \n< 0.01",
                        labels = c("False", "True"),
                        values = c('#41b6c4', '#225ea8'))
volcanoW2_0.01

ggsave("volcano_W2_0-01.svg",
       plot      = volcanoW2_0.01,
       limitsize = FALSE,
       width     = 18,
       height    = 18,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

# Post-infection +10wk time point
dfW10_0.01 <- transform(dfW10,
                        threshold = as.factor(dfW10$FDR < 0.01))

volcanoW10_0.01 <- ggplot(data       = dfW10_0.01,
                          aes(x      = logFC,
                              y      = -log10(FDR),
                              colour = threshold)) +
    geom_point(alpha = 0.4,
               size  = 1.75) +
    xlim(c(-10, 11)) +
    ylim(c(-1, 65)) +
    theme(legend.position   = "right",
          legend.background = element_rect(colour = "black"),
          axis.title        = element_text(size   = 16),
          text              = element_text(size   = 14,
                                           family = "Raavi")) +
    ggtitle("+10 wk") +
    xlab(expression(paste(log[2], "-fold change"))) +
    ylab(expression(paste(-log[10], " adjusted p-value"))) +
    scale_colour_manual("FDR \n< 0.01",
                        labels = c("False", "True"),
                        values = c('#41b6c4', '#225ea8'))
volcanoW10_0.01

ggsave("volcano_W10_0-01.svg",
       plot      = volcanoW10_0.01,
       limitsize = FALSE,
       width     = 20,
       height    = 22,
       dpi       = 300,
       path      = workDir,
       units     = "cm")

######################################################
# Gene Ontology (GO) enrichment analysis of DE genes # ----
######################################################

columns(org.Bt.eg.db)

GO.FDRlogFC.pre1 <- FDRlogFC_DE.pre1
GO.FDRlogFC.pre1$GO.terms <- mapIds(org.Bt.eg.db,
                                    keys      = rownames(GO.FDRlogFC.pre1),
                                    column    = "GO",
                                    keytype   = "ENTREZID",
                                    multiVals = "first")

head(GO.FDRlogFC.pre1)


GO.DE.pre1 <- goana(pre1.lrt,
                      geneid = rownames(pre1.lrt),
                      FDR = 0.05,
                      trend = FALSE,
                      species.KEGG = "bta")

KEGG.DE.pre1 <- kegga(pre1.lrt,
                    geneid = rownames(pre1.lrt),
                    FDR = 0.05,
                    trend = FALSE,
                    species.KEGG = "bta")
head(KEGG.DE.pre1)
top <- topKEGG(KEGG.DE.pre1)

# Count NAs in ENSEMBL tags of Full_DE_PPDb
dim(Full_DE_PPDb)
sum(is.na(Full_DE_PPDb$ENSEMBL.tag))

# Transform rownames into first column because dplyr will discard rownames
# later on and we want to keep this information
full.DE <- Full_DE_PPDb
full.DE <- cbind(RefSeqID = as.numeric(rownames(full.DE)), full.DE)
head(full.DE)

# Subset rows with unknown ENSEMBL tag
NA_DE_PPDb <- dplyr::filter(full.DE, is.na(full.DE$ENSEMBL.tag))
dim(NA_DE_PPDb)
head(NA_DE_PPDb)

# Subset rows with ENSEMBL tag
DE_PPDb_ensembl <- dplyr::filter(full.DE, !is.na(full.DE$ENSEMBL.tag))
dim(DE_PPDb_ensembl)
head(DE_PPDb_ensembl)












####################
# Save .RData file # ----
####################

save.image(file = "PPDb-RNA-seq_paired_sense.RData")

#######################
# Save R session info # ----
#######################

sink("PPDb-RNAseq_paired_sense_R-Session-Info.txt")
sessionInfo()
sink()

#######
# END #
#######
