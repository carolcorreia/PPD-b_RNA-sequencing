######################################################################
# PPD-b stimulated vs unstimulated peripheral blood RNA-seq Analysis #
#                           --- tmod ---                             #
######################################################################

# Author: Carolina N. Correia
# GitHub Repository DOI:
# Date: November 26th 2017

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/PPD-B-sti_vs_uns/edgeR_sense")

# Define variables for specific directories
workDir <- getwd()
tablesDir <- paste0(workDir, "/Tables")
modulesDir <- paste0(workDir, "/Tables/modules")
imgDir <- paste0(workDir, "/Figures")
imgModDir <- paste0(workDir, "/Figures/modules")

# Load previously saved data
load("PPDb-RNA-seq-tmod.rda")

############################################
# 02 Load and/or install required packages #
############################################

# Load packages
library(tidyverse)
library(magrittr)
library(stringr)
library(forcats)
library(tmod)
library(Cairo)
library(reshape2)

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
# 09 Get genes at each LI/DC module #
#####################################

# Define function
GetGenesModules <- function(x, y) {
  # x is a named list of data frames
  # y is character (module code)
  map(x,
      ~tmod::showModule(.x,
                        .x$Ortholog.symbol,
                        module = y,
                        mset = "all")) %>%
    melt() %>%
    mutate(mod_id = y) %>%
    dplyr::rename(time_point = L1) %>%
    mutate(time_point = str_replace(time_point, "_HGNC", "")) %>%
    as.tibble() -> result
}

# Define named list of time points
time_list <- list(Wm1_HGNC, W1_HGNC, W2_HGNC, W10_HGNC)
names(time_list) <- c("Wm1_HGNC", "W1_HGNC", "W2_HGNC", "W10_HGNC")

# All weeks
LI.M27.0 <- GetGenesModules(time_list, "LI.M27.0")
LI.M11.0 <- GetGenesModules(time_list, "LI.M11.0")
LI.S5 <- GetGenesModules(time_list, "LI.S5")
LI.M86.0 <- GetGenesModules(time_list, "LI.M86.0")
LI.M29 <- GetGenesModules(time_list, "LI.M29")
LI.M37.0 <- GetGenesModules(time_list, "LI.M37.0")
LI.M27.1 <- GetGenesModules(time_list, "LI.M27.1")

# Weeks -1, +1, +2
LI.S10 <- GetGenesModules(time_list[1:3], "LI.S10")
LI.M118.0 <- GetGenesModules(time_list[1:3], "LI.M118.0")

# Weeks -1, +2, +10
LI.M43.1 <- GetGenesModules(time_list[-(2)], "LI.M43.1")

# Weeks -1, +1
LI.M4.3 <- GetGenesModules(time_list[1:2], "LI.M4.3")

# Weeks +1, +2, +10
LI.M4.0 <- GetGenesModules(time_list[2:4], "LI.M4.0")
DC.M4.14 <- GetGenesModules(time_list[2:4], "DC.M4.14")

# Weeks +1, +2
LI.M20 <- GetGenesModules(time_list[2:3], "LI.M20")
LI.M112.0 <- GetGenesModules(time_list[2:3], "LI.M112.0")

# Weeks +2, +10
DC.M3.4 <- GetGenesModules(time_list[3:4], "DC.M3.4")
DC.M5.12 <- GetGenesModules(time_list[3:4], "DC.M5.12")
LI.M35.1 <- GetGenesModules(time_list[3:4], "LI.M35.1")
LI.M115 <- GetGenesModules(time_list[3:4], "LI.M115")
LI.M127 <- GetGenesModules(time_list[3:4], "LI.M127")

# Week +1
LI.M0 <- GetGenesModules(time_list[2], "LI.M0")
LI.M11.1 <- GetGenesModules(time_list[2], "LI.M11.1")

# Week +2
LI.M13 <- GetGenesModules(time_list[3], "LI.M13")
LI.M165 <- GetGenesModules(time_list[3], "LI.M165")

# Week +10
LI.M35.0 <- GetGenesModules(time_list[4], "LI.M35.0")
LI.M4.1 <- GetGenesModules(time_list[4], "LI.M4.1")
LI.M86.1 <- GetGenesModules(time_list[4], "LI.M86.1")
LI.M43.0 <- GetGenesModules(time_list[4], "LI.M43.0")
LI.M111.1 <- GetGenesModules(time_list[4], "LI.M111.1")
LI.M24 <- GetGenesModules(time_list[4], "LI.M24")

# Output data
mod_lists <- list(LI.M27.0, LI.M11.0, LI.S5, LI.M86.0, LI.M29,
                  LI.M37.0, LI.S10, LI.M27.1, LI.M4.3, LI.M118.0,
                  LI.M43.1, LI.M20, LI.M4.0, LI.M112.0, DC.M4.14,
                  LI.M0, LI.M11.1, DC.M3.4, DC.M5.12, LI.M35.1,
                  LI.M115, LI.M165, LI.M127, LI.M13, LI.M35.0,
                  LI.M4.1, LI.M86.1, LI.M43.0, LI.M111.1, LI.M24)

mod_files <- c(paste0(c("LIM270", "LIM110", "LIS5", "LIM860", "LIM29",
                        "LIM370", "LIS10", "LIM271", "LIM43", "LIM1180",
                        "LIM431", "LIM20", "LIM40", "LIM1120", "DCM414",
                        "LIM0", "LIM111", "DCM34", "DCM512", "LIM351",
                        "LIM115", "LIM165", "LIM127", "LIM13", "LIM350",
                        "LIM41", "LIM861", "LIM430", "LIM1111", "LIM24"),
                      ".csv"))

pwalk(list(mod_lists, file.path(modulesDir, mod_files)),
      write_csv,
      col_names = TRUE)

######################################
# 10 Plots: heat maps of all modules #
######################################

names(mod_lists) <- c("LIM270", "LIM110", "LIS5", "LIM860", "LIM29",
                      "LIM370", "LIS10", "LIM271", "LIM43", "LIM1180",
                      "LIM431", "LIM20", "LIM40", "LIM1120", "DCM414",
                      "LIM0", "LIM111", "DCM34", "DCM512", "LIM351",
                      "LIM115", "LIM165", "LIM127", "LIM13", "LIM350",
                      "LIM41", "LIM861", "LIM430", "LIM1111", "LIM24")

purpleGreen <- c("#00441B", "#1B7837", "#5AAE61",
                 "#C7E9C0", "white", "#F1B6DA",
                 "#9970AB", "#762A83", "#40004B")

map(mod_lists,
    ~ dplyr::filter(.x, variable == "logFC")) %>%
  map(~ mutate(.x, Ortholog.symbol = fct_rev(Ortholog.symbol))) %>%
  map(~ mutate(.x, time_point = fct_inorder(time_point))) %>%
  map(~ mutate(.x, is.DE = as.logical(is.DE))) %>%
  map(~ ggplot(.x) +
        geom_tile(aes(x = time_point, y = Ortholog.symbol,
                      fill = value),
                  colour = "black",
                  size = 0.2) +
        geom_point(data = .x[.x$is.DE, ],
                   aes(x = time_point, y = Ortholog.symbol,
                       size = as.numeric(.x$is.DE)),
                   colour = "black",
                   shape = 42,
                   size = 4) +
        coord_fixed() +
        ggtitle(unique(.x$mod_id)) +
        xlab(NULL) +
        ylab(NULL) +
        scale_x_discrete(labels = c("Wm1" = "-1 wk",
                                    "W1"  = "+1 wk",
                                    "W2"  = "+2 wk",
                                    "W10" = "+10 wk")) +
        scale_fill_gradientn(name = expression(paste(log[2], "FC")),
                             colours = purpleGreen) +
        theme(panel.border = element_rect(fill = NA,
                                          colour = "black", size = 1),
              axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5),
              axis.text.y = element_text(face = "italic",
                                         hjust = 1),
              text = element_text(size = 14, family = "Calibri"),
              plot.title = element_text(hjust = 0.5))) -> heatmaps

# Check heatmaps
heatmaps[1]

### Export high quality image for all heatmaps
heat_files <- paste0(c(names(heatmaps)), ".pdf")

purrr::pwalk(list(heat_files, heatmaps),
             ggsave,
             device    = cairo_pdf,
             path      = imgModDir ,
             limitsize = FALSE,
             dpi       = 300,
             height    = 15,
             width     = 8,
             units     = "in")

#####################################
# 11 Plots: ROC of LI/DC interferon #
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
# 12 Plots: ROC of LI/DC others #
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
# 13 Save Save .RData file #
############################

save.image("PPDb-RNA-seq-tmod.rda")

##########################
# 14 Save R session info #
##########################

devtools::session_info()

##############################
# Proceed to other heat maps #
##############################

# File: 06-PPDb-RNA-seq_heatmaps.R
