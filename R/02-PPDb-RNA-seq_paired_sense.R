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
# Last updated on 27/10/2017

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

# Uncomment functions below to install packages in case you don't have them

#install.packages("cowplot")
#install.packages("extrafont")
#install.packages("statmod")

###################
# 16 Set up fonts #
###################

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

# Common and trended dispersions are estimated with the
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
              theme_bw(base_size = 14, base_family = "Calibri") +
              ggtitle("Estimated dispersions (NB model)") +
              xlab(expression(paste("Average ", log[2],"CPM"))) +
              ylab("Biological Coefficient of Variation")

PPDb_BCV

# Output high resolution plot
ggsave("PPDb_BCV.pdf",
       plot = PPDb_BCV,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 14,
       units     = "in")

####################
# 25 Fit GLM model #
####################

# Fit a quasi-likelihood negative binomial generalized log-linear model
# for each tag using the design matrix and calculated dispersion
PPDb_fitQL <- glmQLFit(y = PPDb_disp,
                       design = block_animal,
                       robust = TRUE)

names(PPDb_fitQL)
colnames(PPDb_fitQL$design)

###########################
# 26 Plot: QL dispersions #
###########################

# Create a dataframe with the dispersion values
names(PPDb_fitQL)

DispQL <- as.data.frame(cbind(AveLogCPM = PPDb_fitQL$AveLogCPM,
                              deviance = PPDb_fitQL$deviance,
                              df.residual.zeros = PPDb_fitQL$df.residual.zeros,
                              var.prior = PPDb_fitQL$var.prior,
                              var.post = PPDb_fitQL$var.post))

DispQL %<>%
  dplyr::mutate(type_point = "Raw dispersion (NB)") %>%
  dplyr::mutate(type_point2 = "Squeezed EB dispersion") %>%
  dplyr::mutate(type_smooth = "Trended EB dispersion")

head(DispQL)

# Plot all dispersions
PPDb_BCVQL <- ggplot(DispQL) +
                geom_point(aes(x = AveLogCPM,
                               y = sqrt(sqrt(deviance/df.residual.zeros)),
                               fill = type_point),
                           alpha = 0.5) +
                geom_point(aes(x = AveLogCPM,
                               y = sqrt(sqrt(var.post)),
                               colour = type_point2),
                           alpha = 0.5) +
                geom_smooth(aes(x = AveLogCPM,
                                y = sqrt(sqrt(var.prior)),
                                colour = type_smooth),
                            linetype = 2) +
                scale_fill_manual("", values = c("black")) +
                scale_colour_manual("", values = c("indianred4", "blue")) +
                theme_bw(base_size = 14, base_family = "Calibri") +
                ggtitle("Estimated QL dispersions") +
                xlab(expression(paste("Average ", log[2],"CPM"))) +
                ylab("Quarter-Root Mean Deviance")

PPDb_BCVQL

# Output high resolution plot
ggsave("PPDb_BCVQL.pdf",
       plot = PPDb_BCVQL,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 14,
       units     = "in")

################################
# 27 Combine dispersions plots #
################################

# Set grids
BCV_grid <- plot_grid(PPDb_BCV,
                      PPDb_BCVQL,
                      labels = c("A)", "B)"),
                      nrow = 2)

# Check plot
BCV_grid

# Export high quality image
ggsave("BCV_grid.pdf",
       plot      = BCV_grid,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 8,
       units     = "in")

#######################
# 28 Save .RData file #
#######################

save.image(file = "PPDb-RNA-seq_paired_sense.RData")

##########################
# 29 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 3 of this analysis #
######################################

# File: 03-PPDb-RNA-seq_paired_sense.R








