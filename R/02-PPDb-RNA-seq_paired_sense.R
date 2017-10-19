##########################################################################
# RNA-seq analysis: PPD-b stimulated vs unstimulated peripheral blood    #
#                           paired-end reads.                            #
#                                                                        #
#           --- R workflow for analyses of known sense genes ---         #
#                                 Part 2                                 #
##########################################################################

# Based on the workflow created by Nalpas, N.C. (2014)
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474

# Author of current version (4.0.0): Correia, C.N.
# DOI badge of current version:
# Last updated on 19/10/2017

##################################
# 14 Working directory and RData #
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
# 15 Load and/or install required packages #
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
# 16 Set up fonts #
###################

# Import fonts
font_import()

# Registered fonts with R for the PDF output device
loadfonts()

##########################################################
# 17 Tidy DGElist for exploratory data analysis plotting #
##########################################################

tidy_PPDbNorm <-tidy(PPDb_norm, addSamples = TRUE)

# Correct PPDb stimulation info
tidy_PPDbNorm$group %<>%
  stringr::str_replace("PPDbStimulated", "PPD-b-stimulated") %>%
  stringr::str_replace("NonStimulated", "Non-stimulated") %>%
  factor(levels = c("Non-stimulated", "PPD-b-stimulated"))

# Check factors
levels(tidy_PPDbNorm$group)

# Clean animal IDs
tidy_PPDbNorm$animal %<>%
  stringr::str_replace("A", "") %>%
  fct_inorder()

# Check factors
levels(tidy_PPDbNorm$animal)

# Correct time point info
tidy_PPDbNorm$time.point %<>%
  stringr::str_replace("Wm1", "W-1") %>%
  factor(levels = c("W-1", "W1", "W2", "W10"))

# Check factors
levels(tidy_PPDbNorm$time.point)

# Combine animal and time point info for
# plotting labels
tidy_PPDbNorm %<>%
  dplyr::mutate(labels = paste0(time.point, "_", animal))

tidy_PPDbNorm$labels %<>%
  factor(levels = c("W-1_6511", "W-1_6514", "W-1_6520", "W-1_6522", "W-1_6526",
                    "W-1_6635", "W-1_6636", "W-1_6637", "W-1_6644", "W-1_6698",
                    "W1_6511", "W1_6514", "W1_6520", "W1_6522", "W1_6526",
                    "W1_6635", "W1_6636", "W1_6637", "W1_6644", "W1_6698",
                    "W2_6511", "W2_6514", "W2_6520", "W2_6522", "W2_6526",
                    "W2_6635", "W2_6636", "W2_6637", "W2_6644", "W2_6698",
                    "W10_6511", "W10_6514", "W10_6520", "W10_6522", "W10_6526",
                    "W10_6635", "W10_6636", "W10_6637", "W10_6644", "W10_6698"))

# Check factors
levels(tidy_PPDbNorm$labels)

# Check data frame
tidy_PPDbNorm

########################################################
# 18 Plot: density of filtered gene counts per library #
########################################################

ggplot(tidy_PPDbNorm, aes(x = log10(count + 1),
                          y = labels)) +
  scale_y_discrete(limits = rev(levels(tidy_PPDbNorm$labels))) +
  geom_joy(aes(fill = group), alpha = 0.5) +
  scale_fill_manual("Treatment", values = c("#b2b2b2", "#e06377")) +
  theme_bw(base_size = 12, base_family = "Calibri") +
  ggtitle("Density of filtered gene counts per sample") +
  facet_grid(. ~ group) +
  ylab("Time point_Animal number") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_norm


density_norm

# Export high quality image
ggsave("PPDb-density-filt.pdf",
       plot      = density_norm,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 9,
       units     = "in")

####################################
# 19 Plots: MDS at each time point #
####################################

### Plot MDS of Week -1 time point
mds_W_1 <- plotMDS.DGEList(x = PPDb_norm[ , grep(pattern = "_W.1_",
                                                 x = colnames(PPDb_norm))],
                           plot = FALSE,
                           method = "bcv")

names(mds_W_1)
W_1_coord <- mds_W_1$cmdscale.out # Get coords to plot with ggplot2

W_1_coord %<>% # Tidy coords
  tidy() %>%
  dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
  dplyr::mutate(animal = sample, group = sample)

W_1_coord$animal %<>% # Clean animal IDs for plotting
  str_replace("_W.+_\\w", "") %>%
  str_replace("A", "") %>%
  factor()

W_1_coord$group %<>% # Clean group info for plotting
  str_replace("A\\d\\d\\d\\d_.+_", "") %>%
  str_replace("P", "PPD-b-stimulated") %>%
  str_replace("U", "Non-stimulated") %>%
  factor(levels = c("Non-stimulated", "PPD-b-stimulated"))

# Check tidy coords
W_1_coord

# Plot MDS with ggplot2
MDS_W_1 <- ggplot(W_1_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 12, base_family = "Calibri") +
  ggtitle("-1 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

# Check W-1 MDS plot
MDS_W_1


### Plot MDS of Week +1 time point
mds_W1 <- plotMDS.DGEList(x = PPDb_norm[, grep(pattern = "_W1_",
                                               x = colnames(PPDb_norm))],
                          plot = FALSE,
                          method = "bcv")

names(mds_W1)
W1_coord <- mds_W1$cmdscale.out # Get coords to plot with ggplot2

W1_coord %<>% # Tidy coords
  tidy() %>%
  dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
  dplyr::mutate(animal = sample, group = sample)

W1_coord$animal %<>% # Clean animal IDs for plotting
  str_replace("_W.+_\\w", "") %>%
  str_replace("A", "") %>%
  factor()

W1_coord$group %<>% # Clean group info for plotting
  str_replace("A\\d\\d\\d\\d_.+_", "") %>%
  str_replace("P", "PPD-b-stimulated") %>%
  str_replace("U", "Non-stimulated") %>%
  factor(levels = c("Non-stimulated", "PPD-b-stimulated"))

# Check tidy coords
W1_coord

# Plot MDS with ggplot2
MDS_W1 <- ggplot(W1_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 12, base_family = "Calibri") +
  ggtitle("+1 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

# Check W1 MDS plot
MDS_W1


### Plot MDS of Week +2 time point
mds_W2 <- plotMDS.DGEList(x = PPDb_norm[, grep(pattern = "_W2_",
                                       x = colnames(PPDb_norm))],
                          plot = FALSE,
                          method = "bcv")

names(mds_W2)
W2_coord <- mds_W2$cmdscale.out # Get coords to plot with ggplot2

W2_coord %<>% # Tidy coords
  tidy() %>%
  dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
  dplyr::mutate(animal = sample, group = sample)

W2_coord$animal %<>% # Clean animal IDs for plotting
  str_replace("_W.+_\\w", "") %>%
  str_replace("A", "") %>%
  factor()

W2_coord$group %<>% # Clean group info for plotting
  str_replace("A\\d\\d\\d\\d_.+_", "") %>%
  str_replace("P", "PPD-b-stimulated") %>%
  str_replace("U", "Non-stimulated") %>%
  factor(levels = c("Non-stimulated", "PPD-b-stimulated"))

# Check tidy coords
W2_coord

# Plot MDS with ggplot2
MDS_W2 <- ggplot(W2_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 12, base_family = "Calibri") +
  ggtitle("+2 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

# Check W2 MDS plot
MDS_W2


### Plot MDS of Week +10 time point
mds_W10 <- plotMDS.DGEList(x = PPDb_norm[, grep(pattern = "_W10_",
                                        x = colnames(PPDb_norm))],
                           plot = FALSE,
                           method = "bcv")

names(mds_W10)
W10_coord <- mds_W10$cmdscale.out # Get coords to plot with ggplot2

W10_coord %<>% # Tidy coords
  tidy() %>%
  dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
  dplyr::mutate(animal = sample, group = sample)

W10_coord$animal %<>% # Clean animal IDs for plotting
  str_replace("_W.+_\\w", "") %>%
  str_replace("A", "") %>%
  factor()

W10_coord$group %<>% # Clean group info for plotting
  str_replace("A\\d\\d\\d\\d_.+_", "") %>%
  str_replace("P", "PPD-b-stimulated") %>%
  str_replace("U", "Non-stimulated") %>%
  factor(levels = c("Non-stimulated", "PPD-b-stimulated"))

# Check tidy coords
W10_coord

# Plot MDS with ggplot2
MDS_W10 <- ggplot(W10_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 12, base_family = "Calibri") +
  ggtitle("+10 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

# Check W10 MDS plot
MDS_W10


### Export high quality image for all MDS plots
files_MDS <- paste0(c("MDS_W_1", "MDS_W1", "MDS_W2", "MDS_W10"), ".pdf")
plots_MDS <- list(MDS_W_1, MDS_W1, MDS_W2, MDS_W10)

purrr::pwalk(list(files_MDS, plots_MDS),
      ggsave,
      device    = cairo_pdf,
      path      = imgDir,
      limitsize = FALSE,
      dpi       = 300,
      height    = 9,
      width     = 10,
      units     = "in")

#####################################################
# 20 Plot: Combine all MDS plots into single figure #
#####################################################

# Set grid
MDS_grid <- plot_grid((MDS_W_1 + theme(legend.position = "none")),
                      (MDS_W1  + theme(legend.position = "none")),
                      (MDS_W2 + theme(legend.position = "none")),
                      (MDS_W10 + theme(legend.position = "none")),
                      labels = c("A", "B", "C", "D"),
                      ncol = 2,
                      scale = .96)

# Check plot
MDS_grid

# Get legend from one MDS plot
legend <- get_legend(MDS_W_1)

# Add legend to grid
MDS_grid_leg <- plot_grid(MDS_grid, legend, ncol = 1,
                          rel_widths = c(1, .3),
                          rel_heights = c(1, .3))

# Check plot
MDS_grid_leg

# Export high quality image for both plots
files_MDSgrid <- paste0(c("MDS_grid", "MDS_grid_leg"), ".pdf")
plots_MDSgrid <- list(MDS_grid, MDS_grid_leg)

purrr::pwalk(list(files_MDSgrid, plots_MDSgrid),
             ggsave,
             device    = cairo_pdf,
             path      = imgDir,
             limitsize = FALSE,
             dpi       = 300,
             height    = 10,
             width     = 15,
             units     = "in")

##################################
# 21 Define experimental factors #
##################################

head(PPDb_norm$samples)

# PPD-b stimulation factor
condition <- PPDb_norm$samples$group

# Time point factor
time.point <- PPDb_norm$samples$time.point

# Animal factor
animal <- PPDb_norm$samples$animal

# Combine PPD-b stimulation and time point into one factor
# to simplify contrats
cond.time <- factor(paste(PPDb_norm$samples$group,
                          PPDb_norm$samples$time.point,
                          sep="."),
                    levels = c("NonStimulated.Wm1", "NonStimulated.W1",
                               "NonStimulated.W2", "NonStimulated.W10",
                               "PPDbStimulated.Wm1", "PPDbStimulated.W1",
                               "PPDbStimulated.W2", "PPDbStimulated.W10"))

#################################################
# 22 Create a design matrix for paired analysis #
#################################################

# Create a design matrix with animal as a blocking factor
block_animal <- model.matrix(~animal + cond.time,
                             data = PPDb_norm$samples)

dim(block_animal)
dim(PPDb_norm$samples)
head(block_animal)

# Rename design matrix columns for simplicity
colnames(block_animal) %<>%
  str_replace("animal", "") %>%
  str_replace("cond.time", "")

head(block_animal)

# Output the design matrix info
write_csv(as.data.frame(block_animal),
          path = file.path(paste0(workDir, "/PPDb_design-matrix.csv")),
          col_names = TRUE)

#########################################
# 23 Estimate the dispersion parameters #
#########################################

# Common and tended dispersions are estimated with the
# Cox-Reid method and tagwise dispersions with the
# empirical Bayes method
PPDb_disp <- estimateDisp.DGEList(y       = PPDb_norm,
                                  design  = block_animal,
                                  robust  = TRUE,
                                  verbose = TRUE)

names(PPDb_disp)

# Check the calculated dispersion
PPDb_disp$common.dispersion

# Check the calculated dispersion's square root,
# which corresponds to the biological coefficient of variation (BCV)
sqrt(PPDb_disp$common.dispersion)
sqrt(PPDb_disp$tagwise.dispersion)

# Create a matrix of the tagwise dispersion associated with gene information
Tagwisedisp <- cbind(PPDb_disp$genes, PPDb_disp$tagwise.dispersion)
head(Tagwisedisp)
dim(Tagwisedisp)

# Output tagwise dispersion values with gene info
Tagwisedisp <- as.data.frame(cbind(PPDb_disp$genes,
                                   PPDb_disp$tagwise.dispersion))
write_csv(Tagwisedisp,
          path = file.path(paste0(workDir, "/PPDb_Tagwise_dispersion.csv")),
          col_names = TRUE)

################################
# 24 Plot: BCV and dispersions #
################################

# Create a dataframe with the dispersion values
names(PPDb_disp)

Disp <- as.data.frame(cbind(PPDb_disp$genes,
                            PPDb_disp$tagwise.dispersion,
                            PPDb_disp$common.dispersion,
                            PPDb_disp$trended.dispersion,
                            PPDb_disp$AveLogCPM))

colnames(Disp) %<>%
  str_replace("PPDb_disp\\$", "")

Disp %<>%
  dplyr::mutate(type_point = "Tagwise dispersion") %>%
  dplyr::mutate(type_hline = "Common dispersion") %>%
  dplyr::mutate(type_smooth = "Trended dispersion")

# Plot all dispersions
PPDb_BCV <- ggplot(Disp) +
              geom_point(aes(x = AveLogCPM,
                             y = sqrt(tagwise.dispersion),
                             fill = type_point),
                         alpha = 0.5) +
              geom_hline(aes(yintercept = sqrt(common.dispersion),
                             colour = type_hline)) +
              geom_smooth(aes(x = AveLogCPM,
                              y = sqrt(trended.dispersion),
                              colour = type_smooth),
                              linetype = 2) +
              scale_fill_manual("", values = c("black")) +
              scale_colour_manual("", values = c("red", "blue")) +
              theme_bw(base_size = 12, base_family = "Calibri") +
              ggtitle("Estimated dispersions (Cox-Reid method)") +
              xlab(expression(paste(log[2], "Average CPM"))) +
              ylab("Biological Coefficient of Variation")

PPDb_BCV

# Output high resolution plot
ggsave("PPDb_BCV.pdf",
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 14,
       units     = "in")


#####################################################
# 25 Determine differential expression by fitting a #
# negative binomial GLM with Quasi-likelihood Tests #
#####################################################

# Fit a quasi-likelihood negative binomial generalized log-linear model
# for each tag using the design matrix and calculated dispersion
PPDb_fitQL <- glmQLFit(y = PPDb_disp,
                       design = block_animal,
                       robust = TRUE)

names(PPDb_fitQL)
colnames(PPDb_fitQL$design)

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
# 26 Plot: Treemaps of DE genes (FDR < 0.001) #
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
# 27 Plot: Venn diagram of DE genes (FDR < 0.001) #
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
# 28 DE genes common to all contrasts #
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
# 29 Plot: Volcano of DE genes at each time point  #
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
# 30 Plot: Combine all volcanos into single figure #
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
# 31 Save .RData file #
#######################

save.image(file = "PPDb-RNA-seq_paired_sense.RData")

##########################
# 32 Save R session info #
##########################

devtools::session_info()

####################################
# Proceed to Part 3: RNA-seq stats #
####################################

# File: 03-PPDb-RNA-seq_stats.R








