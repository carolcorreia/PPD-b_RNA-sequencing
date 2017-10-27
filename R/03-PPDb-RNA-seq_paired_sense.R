##########################################################################
# RNA-seq analysis: PPD-b stimulated vs unstimulated peripheral blood    #
#                           paired-end reads.                            #
#                                                                        #
#           --- R workflow for analyses of known sense genes ---         #
#                                 Part 3                                 #
##########################################################################

# Based on the workflow created by Nalpas, N.C. (2014)
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474

# Author of current version (4.0.0): Correia, C.N.
# DOI badge of current version:
# Last updated on 27/10/2017

##################################
# 30 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/PPD-B-sti_vs_uns/edgeR_sense")
getwd()

# Define variables for working and file directories
workDir <- getwd()
workDir
fileDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/PPD-B-sti_vs_uns/edgeR_sense/PPDbCounts"
imgDir <- paste0(workDir, "/Figures")
tablesDir <- paste0(workDir, "/Tables")

# Load previously saved data
load("PPDb-RNA-seq_paired_sense.RData")

############################################
# 31 Load and/or install required packages #
############################################

library(statmod)
library(edgeR)
library(devtools)
library(plyr)
library(tidyverse)
library(stringr)
library(magrittr)
library(purrr)
library(forcats)
library(biobroom)
library(ggjoy)
library(ggrepel)
library(Cairo)
library(cowplot)
library(extrafont)
library(VennDiagram)
library(treemap)

# Uncomment functions below to install packages in case you don't have them

#install.packages("cowplot")
#install.packages("extrafont")
#install.packages("statmod")
#install.packages("VennDiagram")
#install.packages("treemap")

###################
# 32 Set up fonts #
###################

# Registered fonts with R for the PDF output device
loadfonts()

#####################################################
# 33 Determine differential expression by fitting a #
# negative binomial GLM with Quasi-likelihood Tests #
#####################################################

# Test for differential expression between the different time points/treatments,
# using the coefficients from PPDbfit$design

# -1 wk
Wm1.QL <- glmQLFTest(PPDb_fitQL, coef = "PPDbStimulated.Wm1")
testDE.Wm1 <- topTags(object        = Wm1.QL,
                      n             = "inf",
                      adjust.method = "BH")

head(testDE.Wm1$table)

# +1 wk
W1.QL <- glmQLFTest(PPDb_fitQL, coef = "PPDbStimulated.W1")
testDE.W1 <- topTags(object        = W1.QL,
                     n             = "inf",
                     adjust.method = "BH")

head(testDE.W1$table)

# +2 wk
W2.QL <- glmQLFTest(PPDb_fitQL, coef = "PPDbStimulated.W2")
testDE.W2 <- topTags(object        = W2.QL,
                     n             = "inf",
                     adjust.method = "BH")

head(testDE.W2$table)

# +10 wk
W10.QL <- glmQLFTest(PPDb_fitQL, coef = "PPDbStimulated.W10")
testDE.W10 <- topTags(object        = W10.QL,
                      n             = "inf",
                      adjust.method = "BH")

head(testDE.W10$table)


### Output all genes tested for DE

# -1 wk
testDE.Wm1$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.001, "true", "false")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(workDir, "/Wm1_AllGenes.csv")),
            col_names = TRUE)

# +1 wk
testDE.W1$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.001, "true", "false")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(workDir, "/W1_AllGenes.csv")),
            col_names = TRUE)

# +2 wk
testDE.W2$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.001, "true", "false")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(workDir, "/W2_AllGenes.csv")),
            col_names = TRUE)

# +10 wk
testDE.W10$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.001, "true", "false")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(workDir, "/W10_AllGenes.csv")),
            col_names = TRUE)


### Filter genes considered DE (FDR < 0.001)

# -1 wk
testDE.Wm1$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::filter(FDR < 0.001) %>%
  dplyr::arrange(desc(logFC)) %>%
  as.tibble() -> Wm1_FDR_001

# +1 wk
testDE.W1$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::filter(FDR < 0.001) %>%
  dplyr::arrange(desc(logFC)) %>%
  as.tibble() -> W1_FDR_001

# +2 wk
testDE.W2$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::filter(FDR < 0.001) %>%
  dplyr::arrange(desc(logFC)) %>%
  as.tibble() -> W2_FDR_001

# +10 wk
testDE.W10$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::filter(FDR < 0.001) %>%
  dplyr::arrange(desc(logFC)) %>%
  as.tibble() -> W10_FDR_001

### Output genes considered DE (FDR < 0.001)
DElists <- list(Wm1_FDR_001, W1_FDR_001, W2_FDR_001, W10_FDR_001)
DEfiles <- c(paste0(c("Wm1_FDR_001", "W1_FDR_001", "W2_FDR_001", "W10_FDR_001"),
                    "_genes.csv"))
DEpaths <- file.path(tablesDir, DEfiles)

pwalk(list(DElists, DEpaths),
      write_csv,
      col_names = TRUE)

###############################################
# 34 Plot: Treemaps of DE genes (FDR < 0.001) #
###############################################

# Get numbers of up and down regulated genes
# at each time point
list_DE <- list(Wm1_FDR_001, W1_FDR_001, W2_FDR_001, W10_FDR_001)
names(list_DE) <- c("-1 wk", "+1 wk", "+2 wk", "+10 wk")

Up_Down <- map_df(list_DE,
                  ~ dplyr::count(.x,
                                 up = sum(logFC > 0),
                                 down = sum(logFC < 0),
                                 zero = sum(logFC == 0)),
                  .id = "time_point")

# Time point as factor
Up_Down$time_point %<>%
  factor() %>%
  fct_inorder()

# Plotting labels
Up_Down %<>% dplyr::mutate(labelsUp = paste(time_point, up, sep = ' '))
Up_Down %<>% dplyr::mutate(labelsDown = paste(time_point, down, sep = ' '))

# Check data frame
Up_Down

# Plot chart increased expression
# Run this chunk together
cairo_pdf(filename = file.path(paste0(imgDir, "/tree_up.pdf")),
          width    = 8,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
treemap(Up_Down,
        index             = "labelsUp",
        vSize             = "up",
        type              = "index",
        palette           = "PRGn",
        title             = "Increased expression",
        fontsize.title    = 14,
        fontfamily.title  = "Calibri",
        fontfamily.labels = "Calibri",
        fontsize.labels   = 16)
dev.off()

# Plot chart decreased expression
# Run this chunk together
cairo_pdf(filename = file.path(paste0(imgDir, "/tree_down.pdf")),
          width    = 8,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
treemap(Up_Down,
        index             = "labelsDown",
        vSize             = "down",
        type              = "index",
        palette           = "-PRGn",
        title             = "Decreased expression",
        fontsize.title    = 14,
        fontfamily.title  = "Calibri",
        fontfamily.labels = "Calibri",
        fontsize.labels   = 16)
dev.off()

###################################################
# 35 Plot: Venn diagram of DE genes (FDR < 0.001) #
###################################################

# Turn gene IDs into vectors
Wm1.vector <- Wm1_FDR_001$EntrezID
W1.vector <- W1_FDR_001$EntrezID
W2.vector <- W2_FDR_001$EntrezID
W10.vector <- W10_FDR_001$EntrezID

# Turn log files off from VennDiagram package
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Plot Venn diagram for all time points
venn.plot <- venn.diagram(list(A = Wm1.vector,
                               B = W10.vector,
                               C = W1.vector,
                               D = W2.vector),
                          filename        = file.path(paste0(imgDir, "/Venn_DE_FDR_001.png")),
                          imagetype       = "png",
                          col             = "transparent",
                          fill            = c("#ffffcc",
                                              "#225ea8",
                                              "#a1dab4",
                                              "#41b6c4"),
                          alpha           = 0.50,
                          label.col       = "#003333",
                          cex             = 3,
                          fontfamily      = "Calibri",
                          category.names  = c("-1 wk",
                                              "+10 wk",
                                              "+1 wk",
                                              "+2 wk"),
                          cat.col         = "black",
                          cat.cex         = 3,
                          cat.pos         = c(-11, 11, 0, 0),
                          cat.dist        = c(0.21, 0.21, 0.1, 0.1),
                          cat.fontfamily  = "Calibri",
                          rotation.degree = 360,
                          margin          = 0,
                          height          = 11,
                          width           = 11,
                          units           = 'in',
                          compression     = 'lzw',
                          resolution      = 300)

#######################################
# 36 DE genes common to all contrasts #
#######################################

# Join common genes into single data frame
Wm1_FDR_001 %>%
  dplyr::inner_join(W1_FDR_001,
                    by = "EntrezID",
                    suffix = c("_Wm1", "_W1")) %>%
  dplyr::inner_join(W2_FDR_001,
                    by = "EntrezID") %>%
  dplyr::inner_join(W10_FDR_001,
                    by = "EntrezID",
                    suffix = c("_W2", "_W10")) -> common_DE

# Check data frame
common_DE

# Output data
write_csv(common_DE,
          path = file.path(paste0(tablesDir, "/common_DE_genes.csv")),
          col_names = TRUE)

####################################################
# 37 Plot: Volcano of DE genes at each time point  #
####################################################

# -1 wk
testDE.Wm1$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.001, "true", "false")) %>%
  ggplot(aes(x = logFC,
             y = -log10(FDR),
             colour = is.DE)) +
  geom_point(alpha = 0.4, size  = 1.75) +
  geom_label_repel(aes(label = if_else(abs(logFC) > 3 & FDR < 0.001,
                                       gene_symbol, NULL),
                       fontface = "italic",
                       family = "Calibri",
                       size = 4,
                       fill = "gray",
                       alpha = 0.5),
                   show.legend = FALSE) +
  xlim(c(-5.5, 5.7)) +
  ylim(c(0, 10)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("-1 wk") +
  xlab(expression(paste(log[2], " fold-change"))) +
  ylab(expression(paste(-log[10], " FDR"))) +
  scale_colour_manual("FDR \n< 0.001",
                      labels = c("False", "True"),
                      values = c('#41b6c4', '#225ea8')) -> Wm1Volcano

Wm1Volcano

# +1 wk
testDE.W1$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.001, "true", "false")) %>%
  ggplot(aes(x = logFC,
             y = -log10(FDR),
             colour = is.DE)) +
  geom_point(alpha = 0.4, size  = 1.75) +
  geom_label_repel(aes(label = if_else(abs(logFC) > 3 & FDR < 0.001,
                                       gene_symbol, NULL),
                       fontface = "italic",
                       family = "Calibri",
                       size = 4,
                       fill = "gray",
                       alpha = 0.5),
                   show.legend = FALSE) +
  xlim(c(-7.5, 7.5)) +
  ylim(c(0, 10)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+1 wk") +
  xlab(expression(paste(log[2], " fold-change"))) +
  ylab(expression(paste(-log[10], " FDR"))) +
  scale_colour_manual("FDR \n< 0.001",
                      labels = c("False", "True"),
                      values = c('#41b6c4', '#225ea8')) -> W1Volcano

W1Volcano

# +2 wk
testDE.W2$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.001, "true", "false")) %>%
  ggplot(aes(x = logFC,
             y = -log10(FDR),
             colour = is.DE)) +
  geom_point(alpha = 0.4, size  = 1.75) +
  geom_label_repel(aes(label = if_else(abs(logFC) > 4 & FDR < 0.001,
                                       gene_symbol, NULL),
                       fontface = "italic",
                       family = "Calibri",
                       size = 4,
                       fill = "gray",
                       alpha = 0.5),
                   show.legend = FALSE) +
  xlim(c(-9.5, 9.5)) +
  ylim(c(0, 14.5)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+2 wk") +
  xlab(expression(paste(log[2], " fold-change"))) +
  ylab(expression(paste(-log[10], " FDR"))) +
  scale_colour_manual("FDR \n< 0.001",
                      labels = c("False", "True"),
                      values = c('#41b6c4', '#225ea8')) -> W2Volcano

W2Volcano

# +10 wk
testDE.W10$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.001, "true", "false")) %>%
  ggplot(aes(x = logFC,
             y = -log10(FDR),
             colour = is.DE)) +
  geom_point(alpha = 0.4, size  = 1.75) +
  geom_label_repel(aes(label = if_else(abs(logFC) > 5 & FDR < 0.001,
                                       gene_symbol, NULL),
                       fontface = "italic",
                       family = "Calibri",
                       size = 4,
                       fill = "gray",
                       alpha = 0.5),
                   show.legend = FALSE) +
  xlim(c(-9.5, 10)) +
  ylim(c(0, 24.5)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+10 wk") +
  xlab(expression(paste(log[2], " fold-change"))) +
  ylab(expression(paste(-log[10], " FDR"))) +
  scale_colour_manual("FDR \n< 0.001",
                      labels = c("False", "True"),
                      values = c('#41b6c4', '#225ea8')) -> W10Volcano

W10Volcano


### Export high quality image for all volcano plots
files_V <- paste0(c("Volcano_W_1", "Volcano_W1",
                    "Volcano_W2", "Volcano_W10"), ".pdf")
plots_V <- list(Wm1Volcano, W1Volcano, W2Volcano, W10Volcano)

purrr::pwalk(list(files_V, plots_V),
             ggsave,
             device    = cairo_pdf,
             path      = imgDir,
             limitsize = FALSE,
             dpi       = 300,
             height    = 12,
             width     = 14,
             units     = "in")

####################################################
# 38 Plot: Combine all volcanos into single figure #
####################################################

# Set grids
V_grid1 <- plot_grid(Wm1Volcano,
                     W1Volcano,
                     labels = c("A", "B"),
                     nrow = 2)

V_grid2 <- plot_grid(W2Volcano,
                     W10Volcano,
                     labels = c("C", "D"),
                     nrow = 2)
# Check plots
V_grid1
V_grid2

# Export high quality image for both grids
files_Vgrid <- paste0(c("V_grid1", "V_grid2"), ".pdf")
plots_Vgrid <- list(V_grid1, V_grid2)

purrr::pwalk(list(files_Vgrid, plots_Vgrid),
             ggsave,
             device    = cairo_pdf,
             path      = imgDir,
             limitsize = FALSE,
             dpi       = 300,
             height    = 15,
             width     = 10,
             units     = "in")

#######################
# 39 Save .RData file #
#######################

save.image(file = "PPDb-RNA-seq_paired_sense.RData")

##########################
# 40 Save R session info #
##########################

devtools::session_info()

############################
# Proceed to RNA-seq stats #
############################

# File: 04-PPDb-RNA-seq_stats.R
