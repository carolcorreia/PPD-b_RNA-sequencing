##########################################################################
# RNA-seq analysis: PPD-b stimulated vs unstimulated peripheral blood    #
#                           paired-end reads.                            #
#                                                                        #
#           --- R workflow for analyses of known sense genes ---         #
#                                 Part 1                                 #
##########################################################################

# Based on the workflow created by Nalpas, N.C. (2014)
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474

# Author of current version (4.0.0): Correia, C.N.
# DOI badge of current version:
# Last updated on 19/09/2017

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
imgDir <- paste0(workDir, "/Figures")

# Load previously saved data
load("PPDb-RNA-seq_paired_sense.RData")

############################################
# 02 Load and/or install required packages #
############################################

library(Biobase)
library(edgeR)
library(AnnotationFuncs)
library(org.Bt.eg.db)
library(devtools)
library(plyr)
library(tidyverse)
library(stringr)
library(magrittr)
library(biobroom)
library(ggjoy)
library(ggrepel)

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
#install.packages("ggrepel")
#install.packages("ggjoy")

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
length(files)

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

# Correct column names in counts
colnames(rawCounts$counts) %<>%
  str_replace("_sense-counts", "") %>%
  str_replace("-", "_")

# Correct row names in counts
rownames(rawCounts$counts) %<>%
  str_replace("BGD.*,", "") %>%
  str_replace(",miRBase:.*", "") %>%
  str_replace("GeneID:", "")

# Correct row names in samples
rownames(rawCounts$samples) %<>%
  str_replace("_sense-counts", "") %>%
  str_replace("-", "_")

# Check data frames
head(rawCounts$counts)
head(rawCounts$samples)

#################################
# 05 Get gene names and symbols #
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
annotCounts %>%
  dplyr::select(gene_symbol, gene_name, everything()) %>%
  rownames_to_column(var = "EntrezID") %>%
  write_csv(file.path(paste0(workDir, "/PPDb-sense_annot-rawcounts.csv")),
            col_names = TRUE)

#########################################
# 06 Add sample information for DGElist #
#########################################

# Treatment group (avoid using underscores)
rawCounts$samples$group <- rownames(rawCounts$samples)
rawCounts$samples$group %<>%
  str_replace("A\\d\\d\\d\\d_.+_", "") %>%
  str_replace("P", "PPDbStimulated") %>%
  str_replace("U", "NonStimulated") %>%
  factor(levels = c("NonStimulated", "PPDbStimulated"))

# Animal ID (cannot start with numbers)
rawCounts$samples$animal <- rownames(rawCounts$samples)
rawCounts$samples$animal %<>%
  str_replace("_W.+_\\w", "") %>%
  factor()

# Time points (avoid using underscores)
rawCounts$samples$time.point <- rownames(rawCounts$samples)
rawCounts$samples$time.point %<>%
  str_replace("A\\d\\d\\d\\d_", "") %>%
  str_replace("_(P|U)", "") %>%
  str_replace("W_1", "Wm1") %>%
  factor(levels = c("Wm1", "W1", "W2", "W10"))

# Check data frame
head(rawCounts$samples)

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

# Export image
ggsave("PPDb-density_plot_raw_counts.png",
       plot      = density_raw,
       device    = "png",
       limitsize = FALSE,
       dpi       = 300,
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

# Ouptut filtered counts
PPDb_filt$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "EntrezID") %>%
  write_csv(file.path(paste0(workDir, "/PPDb-sense_filt_counts.csv")),
            col_names = TRUE)

##############################
# 10 Recompute library sizes #
##############################

PPDb_filt$samples$lib.size <- colSums(PPDb_filt$counts)
head(PPDb_filt$samples)
head(PPDb_dgelist$samples)

###########################################################################
# 11 Calculate normalisation factors using Trimmed Mean of M-values (TMM) #
###########################################################################

# With edgeR, counts are not transformed in any way after
# calculating normalisation factors
PPDb_norm <- calcNormFactors(PPDb_filt, method = "TMM")
head(PPDb_norm$samples)


#######################
# 12 Save .RData file # ----
#######################

save.image(file = "PPDb-RNA-seq_paired_sense.RData")

##########################
# 13 Save R session info # ----
##########################

devtools::session_info()

######################################
# Proceed to Part 2 of this analysis #
######################################

# File: 02-PPDb-RNA-seq_paired_sense.R
