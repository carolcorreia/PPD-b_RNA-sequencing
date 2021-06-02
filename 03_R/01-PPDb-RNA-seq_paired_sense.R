##########################################################################
# RNA-seq analysis: PPD-b stimulated vs unstimulated peripheral blood    #
#                           paired-end reads.                            #
#                                                                        #
#           --- R workflow for analyses of known sense genes ---         #
#                                 Part 1                                 #
##########################################################################

# Author: Carolina N. Correia
# Last updated on 02/06/2021

############################################
# 01 Load and/or install required packages #
############################################

library(tidyverse)
library(magrittr)
library(devtools)
library(here)
library(extrafont)
library(Cairo)
library(Biobase)
library(edgeR)
library(org.Bt.eg.db)
library(biobroom)

# Uncomment functions below to install packages in case you don't have them

# CRAN packages
#install.packages("devtools")
#install.packages("tidyverse")
#install.packages("here")
#install.packages("broom")
#install.packages("extrafont")
#install.packages("Cairo")
#install.packages("remotes")
#remotes::install_github("romainfrancois/nothing")

# Bioconductor and BiocManager
#install.packages("BiocManager")
#BiocManager::install(version = "3.13")

# Bioconductor packages
#BiocManager::install("Biobase")
#BiocManager::install("limma")
#BiocManager::install("edgeR")
#BiocManager::install("geneLenDataBase")
#BiocManager::install("org.Bt.eg.db")
#BiocManager::install("biobroom")

##############################################
# 02 Working directory, fonts, and time zone #
##############################################

# Check working directory
here()

# Define variables for subdirectories
countsDir <- here("PPDbCounts")
imgDir <- here("Figures")
tablesDir <- here("Tables")

# Set time zone
Sys.setenv(TZ = "Europe/Dublin")

# Register fonts with R for the PDF output device
font_import()

# Load fonts
loadfonts()

#########################################
# 03 Import featureCounts sense counts  #
#########################################

# Create a vector of file names
files <- list.files(path        = countsDir,
                    pattern     = "^A6",
                    all.files   = TRUE,
                    full.names  = FALSE,
                    recursive   = FALSE,
                    ignore.case = FALSE)

files
length(files)

# Create a data frame with raw counts for all samples
rawCounts <- readDGE(path         = countsDir,
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
  str_replace(",VGNC:VGNC:.*", "") %>%
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
  write_csv(here(tablesDir, "PPDb-sense_annot-rawcounts.csv"),
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
  str_replace("W_1", "Pre1") %>%
  factor(levels = c("Pre1", "W1", "W2", "W10"))

# Check data frame
rawCounts$samples

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

# Include additional experimental information into DGElist
identical(rownames(rawCounts$samples), rownames(PPDb_dgelist$samples))

PPDb_dgelist$samples$animal <- rawCounts$samples$animal
PPDb_dgelist$samples$time.point <- rawCounts$samples$time.point

#Check samples data frame
head(PPDb_dgelist$samples)

# Check factor levels
levels(PPDb_dgelist$samples$group)
levels(PPDb_dgelist$samples$animal)
levels(PPDb_dgelist$samples$time.point)

################################################
# 08 Density plot: raw gene counts per library #
################################################

# Tidy DGElist and plot data
PPDb_dgelist %>%
  broom::tidy() %>%
  ggplot() +
      geom_density(aes(x     = log10(count + 1),
                       group = sample)) +
      theme_bw(base_size = 16, base_family = "Calibri") +
      ylab("Density of raw gene counts per sample") +
      xlab(expression(paste(log[10], "(counts + 1)"))) -> density_raw


density_raw

# Export image
ggsave("PPDb-density_plot_raw_counts.png",
       path      = imgDir,
       plot      = density_raw,
       device    = "png",
       limitsize = FALSE,
       dpi       = 300)

############################################
# 09 Remove zero and lowly expressed genes #
############################################

# Remove lowly expressed genes
# The filterByExpr function keeps rows that have worthwhile counts in a
# minimum number of samples (two samples in this case because the smallest
# group size is two: PPDbStimulated and NonStimulated)
keep <- filterByExpr(PPDb_dgelist, group = group)
summary(keep)
PPDb_filt <- PPDb_dgelist[keep, , keep.lib.sizes = FALSE]
dim(PPDb_filt$counts)
head(PPDb_filt$counts)

# Ouptut filtered counts
PPDb_filt$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "EntrezID") %>%
  write_csv(here(tablesDir, "PPDb-sense_filt_counts.csv"),
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

##########################
# 12 Save R session info #
##########################

devtools::session_info()

#######################
# 13 Save .RData file #
#######################

# Detach all loaded packages (except base R packages):
require(nothing, quietly = TRUE)

# Save all environment objects:
save.image(file = "PPDb-RNA-seq_paired_sense.RData")

######################################
# Proceed to Part 2 of this analysis #
######################################

# File: 02-PPDb-RNA-seq_paired_sense.R
