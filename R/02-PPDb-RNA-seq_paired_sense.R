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
# Last updated on 04/10/2017

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

# Load previously saved data
load("PPDb-RNA-seq_paired_sense.RData")

############################################
# 15 Load and/or install required packages #
############################################

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

#library(VennDiagram)
#library(grDevices)
#library(grid)
#library(gridExtra)
#library(svglite)
#library(tools)

# Uncomment functions below to install packages in case you don't have them

#install.packages("cowplot")
#install.packages("extrafont")


#install.packages("gdata")
#install.packages("gridExtra")
#install.packages("svglite")
#install.packages("VennDiagram")

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
mds_W_1 <- plotMDS(x = PPDb_norm[ , grep(pattern = "_W.1_",
                                         x = colnames(PPDb_norm))],
                   plot = FALSE)

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
  xlab("Dimension 1") +
  ylab("Dimension 2")

# Check W-1 MDS plot
MDS_W_1


### Plot MDS of Week +1 time point
mds_W1 <- plotMDS(x = PPDb_norm[, grep(pattern = "_W1_",
                             x = colnames(PPDb_norm))],
                  plot = FALSE)

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
  xlab("Dimension 1") +
  ylab("Dimension 2")

# Check W1 MDS plot
MDS_W1


### Plot MDS of Week +2 time point
mds_W2 <- plotMDS(x = PPDb_norm[, grep(pattern = "_W2_",
                                       x = colnames(PPDb_norm))],
                  plot = FALSE)

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
  xlab("Dimension 1") +
  ylab("Dimension 2")

# Check W2 MDS plot
MDS_W2


### Plot MDS of Week +10 time point
mds_W10 <- plotMDS(x = PPDb_norm[, grep(pattern = "_W10_",
                                        x = colnames(PPDb_norm))],
                   plot = FALSE)

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
  xlab("Dimension 1") +
  ylab("Dimension 2")

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
             height    = 14,
             width     = 12,
             units     = "in")

###############################
# 21 Define experimental factors # ----
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
# 17 Create a design matrix for paired analysis # ----
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
# 17 Estimate the dispersion parameter for each tag using the Cox-Reid method # ----
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
# 18 Determine differential expression using negative binomial GLMs # ----
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
# 19 Output lists of DE genes: FDR < 0.05  # ----
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
# 20 Venn diagram of DE genes: FDR < 0.05 # ----
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
                          filename        = "Venn_DE_FDR_0-05.pdf",
                          imagetype       = cairo_pdf,
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
# 21 Volcano plots of DE genes: FDR < 0.05 and absolute logFC > 1 # ----
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
# 22 Combine all volcano plots into one figure # ----
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

devtools::session_info()

####################################
# Proceed to Part 3: RNA-seq stats #
####################################

# File: 03-PPDb-RNA-seq_stats.R








