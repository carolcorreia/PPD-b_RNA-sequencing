######################################################################
# PPD-b stimulated vs unstimulated peripheral blood RNA-seq Analysis #
#                           --- tmod ---                             #
######################################################################

# Author: Carolina N. Correia
# GitHub Repository DOI:
# Date: November 22nd 2017

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
library(Cairo)

# Uncomment functions below to install packages in case you don't have them

#install.packages("tmod")

############################################
# 03 Get human orthologous genes from HCOP #
############################################

# Download Human - Cow orthologues (HomoloGene only) using the
# HGNC Comparison of Orthology Predictions (HCOP) search tool:
# https://www.genenames.org/cgi-bin/hcop
# Downloaded on 14th October 2017

hcop_data <- read.delim2("hcop_homologene.txt")

# Check data frame
View(hcop_data)
colnames(hcop_data)

# Select desired variables
hcop_data %>%
  dplyr::select(`Primary.NCBI.Gene.ID`,
                `Ortholog.NCBI.Gene.ID`,
                `Ortholog.symbol`) %>%
  as.tibble() -> HGNC_symbols

HGNC_symbols$Primary.NCBI.Gene.ID %<>% as.character()
HGNC_symbols$Ortholog.NCBI.Gene.ID %<>% as.character()

# Check data frame
HGNC_symbols

################################################
# 04 Map cow gene symbols to human orthologues #
################################################

# Import all genes previously tested for DE with edgeR
# and sort by smallest FDR

# -1 wk
read_csv("Wm1_AllGenes.csv", col_types = cols(EntrezID = "c")) %>%
  dplyr::inner_join(HGNC_symbols,
                   by = c("EntrezID" = "Primary.NCBI.Gene.ID")) %>%
  dplyr::filter(!`Ortholog.symbol` == "") %>%
  dplyr::filter(duplicated(`Ortholog.symbol`) == FALSE) %>%
  dplyr::arrange(FDR) -> Wm1_HGNC


# +1 wk
read_csv("W1_AllGenes.csv", col_types = cols(EntrezID = "c")) %>%
  dplyr::inner_join(HGNC_symbols,
                    by = c("EntrezID" = "Primary.NCBI.Gene.ID")) %>%
  dplyr::filter(!`Ortholog.symbol` == "") %>%
  dplyr::filter(duplicated(`Ortholog.symbol`) == FALSE) %>%
  dplyr::arrange(FDR) -> W1_HGNC


# +2 wk
read_csv("W2_AllGenes.csv", col_types = cols(EntrezID = "c")) %>%
  dplyr::inner_join(HGNC_symbols,
                    by = c("EntrezID" = "Primary.NCBI.Gene.ID")) %>%
  dplyr::filter(!`Ortholog.symbol` == "") %>%
  dplyr::filter(duplicated(`Ortholog.symbol`) == FALSE) %>%
  dplyr::arrange(FDR) -> W2_HGNC


# +10 wk
read_csv("W10_AllGenes.csv", col_types = cols(EntrezID = "c")) %>%
  dplyr::inner_join(HGNC_symbols,
                    by = c("EntrezID" = "Primary.NCBI.Gene.ID")) %>%
  dplyr::filter(!`Ortholog.symbol` == "") %>%
  dplyr::filter(duplicated(`Ortholog.symbol`) == FALSE) %>%
  dplyr::arrange(FDR) -> W10_HGNC


# Check data frames
Wm1_HGNC
W1_HGNC
W2_HGNC
W10_HGNC

# Output data
hgnclists <- list(Wm1_HGNC, W1_HGNC, W2_HGNC, W10_HGNC)
hgncfiles <- c(paste0(c("Wm1_HGNC", "W1_HGNC", "W2_HGNC", "W10_HGNC"),
                      ".csv"))

pwalk(list(hgnclists, file.path(tablesDir, hgncfiles)),
      write_csv,
      col_names = TRUE)

###############################################
# 05 Apply tmod CERNO test with LI/DC modules #
###############################################

# Apply the CERNO test at each time point with FDR < 0.001

# -1 wk
Wm1_tmod <- tmodCERNOtest(Wm1_HGNC$Ortholog.symbol,
                          qval = 0.001,
                          order.by = "pval",
                          mset = "all")

head(Wm1_tmod)
dim(Wm1_tmod)

# +1 wk
W1_tmod <- tmodCERNOtest(W1_HGNC$Ortholog.symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = "all")

head(W1_tmod)
dim(W1_tmod)

# +2 wk
W2_tmod <- tmodCERNOtest(W2_HGNC$Ortholog.symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = "all")

head(W2_tmod)
dim(W2_tmod)

# +10 wk
W10_tmod <- tmodCERNOtest(W10_HGNC$Ortholog.symbol,
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
# (the universe of 11,510 ortholog genes is the same for all time points,
# so it does not matter which one is used as the input character vector
# of gene symbols)
pie <- tmodDecideTests(Wm1_HGNC$Ortholog.symbol,
                       lfc = cbind(Wm1_HGNC$logFC,
                                   W1_HGNC$logFC,
                                   W2_HGNC$logFC,
                                   W10_HGNC$logFC),
                       pval = cbind(Wm1_HGNC$FDR,
                                    W1_HGNC$FDR,
                                    W2_HGNC$FDR,
                                    W10_HGNC$FDR),
                       lfc.thr = 0.5,
                       pval.thr = 0.001,
                       mset = "all")

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

#####################################
# 08 Plots: ROC of LI/DC chemokines #
#####################################

# -1 wk
cairo_pdf(filename = file.path(paste0(imgDir, "/chemoLI_Wm1.pdf")),
          width    = 5,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
evidencePlot(Wm1_HGNC$Ortholog.symbol,
             legend = "right",
             main = "-1 wk",
             c("LI.M27.0", "LI.M27.1", "LI.M29", "LI.M86.0"),
             mset = "all",
             col = c("#67001f", "#d6604d",
                     "#5aae61", "#053061"))
dev.off()

# +1 wk
cairo_pdf(filename = file.path(paste0(imgDir, "/chemoLI_W1.pdf")),
          width    = 5,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
evidencePlot(W1_HGNC$Ortholog.symbol,
             legend = "right",
             main = "+1 wk",
             c("LI.M27.0", "LI.M27.1", "LI.M29", "LI.M86.0"),
             mset = "all",
             col = c("#67001f", "#d6604d",
                     "#5aae61", "#053061"))
dev.off()

# +2 wk
cairo_pdf(filename = file.path(paste0(imgDir, "/chemoLI_W2.pdf")),
          width    = 5,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
evidencePlot(W2_HGNC$Ortholog.symbol,
             legend = "right",
             main = "+2 wk",
             c("LI.M27.0", "LI.M27.1", "LI.M29", "LI.M86.0"),
             mset = "all",
             col = c("#67001f", "#d6604d",
                     "#5aae61", "#053061"))
dev.off()

# +10 wk
cairo_pdf(filename = file.path(paste0(imgDir, "/chemoLI_W10.pdf")),
          width    = 5,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
evidencePlot(W10_HGNC$Ortholog.symbol,
             legend = "right",
             main = "+10 wk",
             c("LI.M27.0", "LI.M27.1", "LI.M29", "LI.M86.0"),
             mset = "all",
             col = c("#67001f", "#d6604d",
                     "#5aae61", "#053061"))
dev.off()

#####################################
# 09 Plots: ROC of LI/DC interferon #
#####################################

# +2 wk
cairo_pdf(filename = file.path(paste0(imgDir, "/interLI_W2.pdf")),
          width    = 5,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
evidencePlot(W2_HGNC$Ortholog.symbol,
             legend = "right",
             main = "+2 wk",
             c("LI.M127", "DC.M3.4", "DC.M5.12"),
             mset = "all",
             col = c("#67001f", "#5aae61", "#053061"))
dev.off()

# +10 wk
cairo_pdf(filename = file.path(paste0(imgDir, "/interLI_W10.pdf")),
          width    = 5,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
evidencePlot(W10_HGNC$Ortholog.symbol,
             legend = "right",
             main = "+10 wk",
             c("LI.M127", "LI.M111.1", "DC.M3.4", "DC.M5.12"),
             mset = "all",
             col = c("#67001f", "#d6604d",
                     "#5aae61", "#053061"))
dev.off()

#################################
# 10 Plots: ROC of LI/DC others #
#################################

# +1 wk
cairo_pdf(filename = file.path(paste0(imgDir, "/otherLI_W1.pdf")),
          width    = 5,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
evidencePlot(W1_HGNC$Ortholog.symbol,
             legend = "right",
             main = "+1 wk",
             c("DC.M4.14"),
             mset = "all",
             col = c("#67001f"))
dev.off()

# +2 wk
cairo_pdf(filename = file.path(paste0(imgDir, "/otherLI_W2.pdf")),
          width    = 5,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
evidencePlot(W2_HGNC$Ortholog.symbol,
             legend = "right",
             main = "+2 wk",
             c("LI.M35.1", "LI.M115", "DC.M4.14"),
             mset = "all",
             col = c("#67001f", "#5aae61", "#053061"))
dev.off()

# +10 wk
cairo_pdf(filename = file.path(paste0(imgDir, "/otherLI_W10.pdf")),
          width    = 5,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
evidencePlot(W10_HGNC$Ortholog.symbol,
             legend = "right",
             main = "+10 wk",
             c("LI.M24", "LI.M35.1", "LI.M115", "DC.M4.14"),
             mset = "all",
             col = c("#67001f", "#d6604d",
                     "#5aae61", "#053061"))
dev.off()

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
Wm1_msig <- tmodCERNOtest(Wm1_HGNC$Ortholog.symbol,
                          qval = 0.001,
                          order.by = "pval",
                          mset = msig[hall])

head(Wm1_msig)
dim(Wm1_msig)

# +1 wk
W1_msig <- tmodCERNOtest(W1_HGNC$Ortholog.symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = msig[hall])

head(W1_msig)
dim(W1_msig)

# +2 wk
W2_msig <- tmodCERNOtest(W2_HGNC$Ortholog.symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = msig[hall])

head(W2_msig)
dim(W2_msig)

# +10 wk
W10_msig <- tmodCERNOtest(W10_HGNC$Ortholog.symbol,
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

##################################################
# 11 Compare Hallmark modules across time points #
##################################################

# Named msig lists
msiglists <- list(Wm1_msig, W1_msig, W2_msig, W10_msig)
names(msiglists) <- c("Wm1", "W1", "W2", "W10")

# Summarise module information and order it by q-values
hall_summary <- tmodSummary(msiglists, clust = "qval")

# Check data frame
head(hall_summary)
dim(hall_summary)

# Output data
hall_summary %>%
  dplyr::rename(Module_id = ID) %>%
  write_csv(path = file.path(paste0(tablesDir, "/hallmark_summary.csv")),
            col_names = TRUE)

#########################################################
# 12 Plot: panel of Hallmark modules across time points #
#########################################################

# Create data frame for pie object
# (the universe of 17,789 ortholog genes is the same for all time points,
# so it does not matter which one is used as the input character vector
# of gene symbols)
pie_hall <- tmodDecideTests(Wm1_HGNC$Ortholog.symbol,
                            lfc = cbind(Wm1_HGNC$logFC,
                                        W1_HGNC$logFC,
                                        W2_HGNC$logFC,
                                        W10_HGNC$logFC),
                            pval = cbind(Wm1_HGNC$FDR,
                                         W1_HGNC$FDR,
                                         W2_HGNC$FDR,
                                         W10_HGNC$FDR),
                            lfc.thr = 0.5,
                            pval.thr = 0.001,
                            mset = msig[hall])

names(pie_hall) <- c("Wm1", "W1", "W2", "W10")

# Run this chunk together
cairo_pdf(filename = file.path(paste0(imgDir, "/Hall_panel_plot.pdf")),
          width    = 10,
          height   = 12,
          family   = "Calibri",
          fallback_resolution = 300)
tmodPanelPlot(msiglists,
              pval.thr = 10^-3,
              pie = pie_hall,
              pie.style = "pie",
              pie.colors = c("#008837", "#e5e5e5", "#7b3294"),
              clust = "effect")
dev.off()

########################################################
# 13 Plots: evidence plot of selected Hallmark modules #
########################################################

# -1 wk
par(family = "Calibri")
evidencePlot(Wm1_HGNC$Ortholog.symbol,
             legend = "right",
             main = "-1 wk",
             lwd = 1.8,
             c("LI.M27.0", "LI.M27.1", "LI.M24",
               "LI.M29", "LI.M86.0", "LI.M115"),
             mset = msig[hall],
             col = c("#67001f", "#b2182b", "#d6604d",
                     "#5aae61", "#2166ac", "#053061"))

# +10 wk
par(family = "Calibri")
evidencePlot(W10_HGNC$Ortholog.symbol,
             legend = "right",
             main = "+10 wk",
             lwd = 1.8,
             c("LI.M27.0", "LI.M27.1", "LI.M24",
               "LI.M29", "LI.M86.0", "LI.M115"),
             mset = msig[hall],
             col = c("#67001f", "#b2182b", "#d6604d",
                     "#5aae61", "#2166ac", "#053061"))


#################################
# 14 Genes at each LI/DC module #
#################################

# -1 wk
View(rownames(tmod::showModule(Wm1_HGNC,
                               Wm1_HGNC$Ortholog.symbol,
                               module = c("LI.M27.0"))))
View(rownames(tmod::showModule(Wm1_HGNC,
                               Wm1_HGNC$Ortholog.symbol,
                               module = c("LI.M27.1"))))

View(rownames(tmod::showModule(Wm1_HGNC,
                               Wm1_HGNC$Ortholog.symbol,
                               module = c("LI.M29"))))
View(rownames(tmod::showModule(Wm1_HGNC,
                               Wm1_HGNC$Ortholog.symbol,
                               module = c("LI.M86.0"))))


LI.M11.0
LI.S5
LI.M86.0
LI.M29
LI.M37.0
LI.S10
LI.M4.3
LI.M43.1
LI.M118.0
DC.M5.14


######################################################
# 14 Apply tmod CERNO test with MSigDB KEGG pathways #
######################################################

# Apply the CERNO test at each time point with FDR < 0.001

# -1 wk
Wm1_kegg <- tmodCERNOtest(Wm1_HGNC$Ortholog.symbol,
                          qval = 0.001,
                          order.by = "pval",
                          mset = msig[kegg])

head(Wm1_kegg)
dim(Wm1_kegg)

# +1 wk
W1_kegg <- tmodCERNOtest(W1_HGNC$Ortholog.symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = msig[kegg])

head(W1_kegg)
dim(W1_kegg)

# +2 wk
W2_kegg <- tmodCERNOtest(W2_HGNC$Ortholog.symbol,
                         qval = 0.001,
                         order.by = "pval",
                         mset = msig[kegg])

head(W2_kegg)
dim(W2_kegg)

# +10 wk
W10_kegg <- tmodCERNOtest(W10_HGNC$Ortholog.symbol,
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

###############################################
# 15 Compare KEGG pathways across time points #
###############################################

# Named KEGG lists
kegglists <- list(Wm1_kegg, W1_kegg, W2_kegg, W10_kegg)
names(kegglists) <- c("Wm1", "W1", "W2", "W10")

# Summarise module information and order it by q-values
kegg_summary <- tmodSummary(kegglists, clust = "qval")

# Check data frame
head(kegg_summary)
dim(kegg_summary)

# Output data
kegg_summary %>%
  dplyr::rename(Module_id = ID) %>%
  write_csv(path = file.path(paste0(tablesDir, "/kegg_summary.csv")),
            col_names = TRUE)

######################################################
# 16 Plot: panel of KEGG pathways across time points #
######################################################

# Create data frame for pie object
# (the universe of 17,789 ortholog genes is the same for all time points,
# so it does not matter which one is used as the input character vector
# of gene symbols)
pie_kegg <- tmodDecideTests(Wm1_HGNC$Ortholog.symbol,
                            lfc = cbind(Wm1_HGNC$logFC,
                                        W1_HGNC$logFC,
                                        W2_HGNC$logFC,
                                        W10_HGNC$logFC),
                            pval = cbind(Wm1_HGNC$FDR,
                                         W1_HGNC$FDR,
                                         W2_HGNC$FDR,
                                         W10_HGNC$FDR),
                            lfc.thr = 0.5,
                            pval.thr = 0.001,
                            mset = msig[kegg])

names(pie_kegg) <- c("Wm1", "W1", "W2", "W10")

# Run this chunk together
cairo_pdf(filename = file.path(paste0(imgDir, "/kegg_panel_plot.pdf")),
          width    = 10,
          height   = 12,
          family   = "Calibri",
          fallback_resolution = 300)
tmodPanelPlot(kegglists,
              pval.thr = 10^-3,
              pie = pie_kegg,
              pie.style = "pie",
              pie.colors = c("#008837", "#e5e5e5", "#7b3294"),
              clust = "effect")
dev.off()

####################################################
# 17 Plots: evidence plot of selected KEGG modules #
####################################################

# +10 wk
par(family = "Calibri")
evidencePlot(W10_HGNC$Ortholog.symbol,
             legend = "right",
             main = "+10 wk",
             lwd = 1.8,
             c("LI.M27.0", "LI.M27.1", "LI.M24",
               "LI.M29", "LI.M86.0", "LI.M115"),
             mset = msig[kegg],
             col = c("#67001f", "#b2182b", "#d6604d",
                     "#5aae61", "#2166ac", "#053061"))


############################
# 18 Save Save .RData file #
############################

save.image("PPDb-RNA-seq-tmod.rda")

##########################
# 19 Save R session info #
##########################

devtools::session_info()

########################
# Proceed to heat maps #
########################

# File: 06-PPDb-RNA-seq_heatmaps.R
