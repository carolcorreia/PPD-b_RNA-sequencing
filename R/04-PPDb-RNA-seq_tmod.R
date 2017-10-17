######################################################################
# PPD-b stimulated vs unstimulated peripheral blood RNA-seq Analysis #
#                    --- tmod ---                      #
######################################################################

# Author: Carolina N. Correia
# GitHub Repository DOI:
# Date: October 17th 2017

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/PPD-B-sti_vs_uns/edgeR_sense")

# Define variables for specific directories
workDir <- getwd()
tablesDir <- paste0(workDir, "/Tables")
imgDir <- paste0(workDir, "/Figures")

# Load previously saved data
load("PPDb-RNA-seq-tmod.rda")

############################################
# 02 Load and/or install required packages #
############################################

# Load packages
library(tidyverse)
library(magrittr)
library(tmod)


# Uncomment functions below to install packages in case you don't have them

#install.packages("tmod")

############################################
# 03 Get human orthologous genes from HCOP #
############################################

# Download Human - Cow orthologues (15 column file) using the
# HGNC Comparison of Orthology Predictions (HCOP) bulk downloads tool:
# "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_cow_hcop_fifteen_column.txt.gz"
# Downloaded on 14th October 2017

hcop_data <- read_delim("human_cow_hcop_fifteen_column.txt",
                        delim = "\t",
                        quote = "")

# Check data frame
View(hcop_data)

# Select desired variables
hcop_data %>%
  dplyr::select(cow_entrez_gene, cow_symbol,
                human_entrez_gene, human_symbol) -> HGNC_symbols

# Check data frame
View(HGNC_symbols)

##################################################
# 04 Match cow gene symbols to human orthologues #
##################################################

# Import all genes previously tested for DE with edgeR
# and sort by smallest FDR

# -1 wk
read_csv("Wm1_AllGenes.csv", col_types = cols(EntrezID = "c")) %>%
  dplyr::inner_join(HGNC_symbols,
                   by = c("EntrezID" = "cow_entrez_gene")) %>%
  dplyr::arrange(FDR) -> Wm1_HGNC


# +1 wk
read_csv("W1_AllGenes.csv", col_types = cols(EntrezID = "c")) %>%
  dplyr::inner_join(HGNC_symbols,
                    by = c("EntrezID" = "cow_entrez_gene")) %>%
  dplyr::arrange(FDR) -> W1_HGNC


# +2 wk
read_csv("W2_AllGenes.csv", col_types = cols(EntrezID = "c")) %>%
  dplyr::inner_join(HGNC_symbols,
                    by = c("EntrezID" = "cow_entrez_gene")) %>%
  dplyr::arrange(FDR) -> W2_HGNC


# +10 wk
read_csv("W10_AllGenes.csv", col_types = cols(EntrezID = "c")) %>%
  dplyr::inner_join(HGNC_symbols,
                    by = c("EntrezID" = "cow_entrez_gene")) %>%
  dplyr::arrange(FDR) -> W10_HGNC


# Check data frames
Wm1_HGNC
W1_HGNC
W2_HGNC
W10_HGNC

##################################################
# 05 Apply tmod CERNO test with standard modules #
##################################################

# Apply the CERNO test at each time point with FDR < 0.001

# -1 wk
Wm1_tmod <- tmodCERNOtest(Wm1_HGNC$human_symbol,
                          qval = 0.001,
                          order.by = "pval",
                          mset = "all")

head(Wm1_tmod)
dim(Wm1_tmod)

# +1 wk
W1_tmod <- tmodCERNOtest(W1_HGNC$human_symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = "all")

head(W1_tmod)
dim(W1_tmod)

# +2 wk
W2_tmod <- tmodCERNOtest(W2_HGNC$human_symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = "all")

head(W2_tmod)
dim(W2_tmod)

# +10 wk
W10_tmod <- tmodCERNOtest(W10_HGNC$human_symbol,
                          qval = 0.001,
                          order.by = "pval",
                          mset = "all")

head(W10_tmod)
dim(W10_tmod)

# Output data
tmodlists <- list(Wm1_tmod, W1_tmod, W2_tmod, W10_tmod)
IDcols <- rep(c("Module_id"), each = 4)
tmodfiles <- c(paste0(c("Wm1_tmod", "W1_tmod", "W2_tmod", "W10_tmod"),
                    ".csv"))

tmodlists %<>% map2(.y = IDcols, ~ dplyr::rename(.x, !!.y := "ID"))

pwalk(list(tmodlists, file.path(tablesDir, tmodfiles)),
      write_csv,
      col_names = TRUE)

#########################################
# 06 Compare modules across time points #
#########################################

# Create a named list with the CERNO test results
time_points <- list(Wm1_tmod, W1_tmod, W2_tmod, W10_tmod)
names(time_points) <- c("Wm1", "W1", "W2", "W10")

# Summarise module information and order it by q-values
tmod_summary <- tmodSummary(time_points, clust = "qval")

# Check data frame
head(tmod_summary)
dim(tmod_summary)

# Output data
tmod_summary %>%
  dplyr::rename(Module_id = ID) %>%
  write_csv(path = file.path(paste0(tablesDir, "/tmod_summary.csv")),
            col_names = TRUE)

################################################
# 07 Plot: panel of modules across time points #
################################################

# Create data frame for pie object
# (the universe of 17,789 ortholog genes is the same for all time points,
# so it does not matter which one is used as the input character vector
# of gene symbols)
pie <- tmodDecideTests(Wm1_HGNC$human_symbol,
                       lfc = cbind(Wm1_HGNC$logFC,
                                   W1_HGNC$logFC,
                                   W2_HGNC$logFC,
                                   W10_HGNC$logFC),
                       pval = cbind(Wm1_HGNC$FDR,
                                    W1_HGNC$FDR,
                                    W2_HGNC$FDR,
                                    W10_HGNC$FDR),
                       lfc.thr = 0.5,
                       pval.thr = 0.001)

names(pie) <- c("Wm1", "W1", "W2", "W10")

# Run this chunk together
cairo_pdf(filename = file.path(paste0(imgDir, "/panel_plot.pdf")),
          width    = 10,
          height   = 12,
          family   = "Calibri",
          fallback_resolution = 300)
tmodPanelPlot(time_points,
              pval.thr = 10^-3,
              pie = pie,
              pie.style = "pie",
              pie.colors = c("#008837", "#e5e5e5", "#7b3294"),
              clust = "effect")
dev.off()

###############################################
# 08 Plots: evidence plot of selected modules #
###############################################

# +10 wk
evidencePlot(W10_HGNC$human_symbol, "LI.M16")
evidencePlot(W10_HGNC$human_symbol, "LI.M27.0")

############################
# 09 Import MSigDB modules #
############################

# http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.1/msigdb_v6.1.xml
# Current MSigDB xml file
# 17/10/17

msig <- tmodImportMSigDB("msigdb_v6.1.xml")
msig

# Subset hallmark genes
hall <- msig$MODULES$Category == "H"

# Subset kegg pathways
kegg <- msig$MODULES$Subcategory == "CP:KEGG"

#######################################################
# 10 Apply tmod CERNO test with MSigDB hallmark genes #
#######################################################

# Apply the CERNO test at each time point with FDR < 0.001

# -1 wk
Wm1_msig <- tmodCERNOtest(Wm1_HGNC$human_symbol,
                          qval = 0.001,
                          order.by = "pval",
                          mset = msig[hall])

head(Wm1_msig)
dim(Wm1_msig)

# +1 wk
W1_msig <- tmodCERNOtest(W1_HGNC$human_symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = msig[hall])

head(W1_msig)
dim(W1_msig)

# +2 wk
W2_msig <- tmodCERNOtest(W2_HGNC$human_symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = msig[hall])

head(W2_msig)
dim(W2_msig)

# +10 wk
W10_msig <- tmodCERNOtest(W10_HGNC$human_symbol,
                          qval = 0.001,
                          order.by = "pval",
                          mset = msig[hall])

head(W10_msig)
dim(W10_msig)

# Output data
msiglists <- list(Wm1_msig, W1_msig, W2_msig, W10_msig)
msigIDcols <- rep(c("Module_id_msig"), each = 4)
msigfiles <- c(paste0(c("Wm1_msig", "W1_msig", "W2_msig", "W10_msig"),
                      ".csv"))

msiglists %<>% map2(.y = msigIDcols, ~ dplyr::rename(.x, !!.y := "ID"))

pwalk(list(msiglists, file.path(tablesDir, msigfiles)),
      write_csv,
      col_names = TRUE)

######################################################
# 11 Apply tmod CERNO test with MSigDB KEGG pathways #
######################################################

# Apply the CERNO test at each time point with FDR < 0.001

# -1 wk
Wm1_kegg <- tmodCERNOtest(Wm1_HGNC$human_symbol,
                          qval = 0.001,
                          order.by = "pval",
                          mset = msig[kegg])

head(Wm1_kegg)
dim(Wm1_kegg)

# +1 wk
W1_kegg <- tmodCERNOtest(W1_HGNC$human_symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = msig[kegg])

head(W1_kegg)
dim(W1_kegg)

# +2 wk
W2_kegg <- tmodCERNOtest(W2_HGNC$human_symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = msig[kegg])

head(W2_kegg)
dim(W2_kegg)

# +10 wk
W10_kegg <- tmodCERNOtest(W10_HGNC$human_symbol,
                          qval = 0.001,
                          order.by = "pval",
                          mset = msig[kegg])

head(W10_kegg)
dim(W10_kegg)

# Output data
kegglists <- list(Wm1_kegg, W1_kegg, W2_kegg, W10_kegg)
keggIDcols <- rep(c("Module_id_kegg"), each = 4)
keggfiles <- c(paste0(c("Wm1_kegg", "W1_kegg", "W2_kegg", "W10_kegg"),
                      ".csv"))

kegglists %<>% map2(.y = keggIDcols, ~ dplyr::rename(.x, !!.y := "ID"))

pwalk(list(kegglists, file.path(tablesDir, keggfiles)),
      write_csv,
      col_names = TRUE)



############################
# 12 Save Save .RData file #
############################

save.image("PPDb-RNA-seq-tmod.rda")

##########################
# 13 Save R session info #
##########################

devtools::session_info()

#######
# END #
#######