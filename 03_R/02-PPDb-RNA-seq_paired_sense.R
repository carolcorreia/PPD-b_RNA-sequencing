##########################################################################
# RNA-seq analysis: PPD-b stimulated vs unstimulated peripheral blood    #
#                           paired-end reads.                            #
#                                                                        #
#           --- R workflow for analyses of known sense genes ---         #
#                                 Part 2                                 #
##########################################################################

# Author: Carolina N. Correia
# Last updated on 02/06/2021

############################################
# 14 Load and/or install required packages #
############################################

library(statmod)
library(edgeR)
library(devtools)
library(tidyverse)
library(magrittr)
library(biobroom)
library(ggridges)
library(ggrepel)
library(Cairo)
library(cowplot)
library(extrafont)
library(here)
library(ggfortify)

# Uncomment functions below to install packages in case you don't have them

#install.packages("statmod")
#install.packages("cowplot")
#install.packages("ggfortify")
#install.packages("ggrepel")
#install.packages("ggridges")

##############################################
# 15 Working directory, fonts, and time zone #
##############################################

# Check working directory
here()

# Load previously saved R data and environment
load("PPDb-RNA-seq_paired_sense.RData")

# Check variables for subdirectories
imgDir
tablesDir

# Set time zone
Sys.setenv(TZ = "Europe/Dublin")

# Load registered fonts with R for the PDF output device
loadfonts()

##########################################################
# 16 Tidy DGElist for exploratory data analysis plotting #
##########################################################

tidy_PPDbNorm <- broom::tidy(PPDb_norm, addSamples = TRUE)

# Correct PPD-b stimulation info
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
  stringr::str_replace("Pre1", "W-1") %>%
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
# 17 Plot: density of filtered gene counts per library #
########################################################

ggplot(tidy_PPDbNorm, aes(x = log10(count + 1),
                          y = labels)) +
  scale_y_discrete(limits = rev(levels(tidy_PPDbNorm$labels))) +
  geom_density_ridges(aes(fill = group), alpha = 0.5) +
  scale_fill_manual("Treatment", values = c("#808080", "#e06377")) +
  theme_bw(base_size = 16, base_family = "Calibri") +
  ggtitle("Density of filtered gene counts per sample") +
  facet_grid(. ~ group) +
  ylab("Time point_Animal number") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_norm


density_norm

# Export high quality image
ggsave("PPDb-density-filt.pdf",
       path      = imgDir,
       plot      = density_norm,
       device    = cairo_pdf,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 9,
       units     = "in")

####################################
# 18 Plots: PCA at each time point #
####################################

############ One week pre-infection
# Calculate PCs using log2CPM filtered counts
Pre1_pca <- prcomp(t(cpm(PPDb_norm[ , grep(pattern = "_W_1_", x = colnames(PPDb_norm))],
                    log = TRUE)),
              scale. = TRUE)

# Get sample names for plotting
Pre1_samples <- PPDb_norm$samples
Pre1_samples %<>%
  rownames_to_column(var = "sample_name") %>%
  dplyr::filter(stringr::str_detect(sample_name, "_W_1_"))

# Plot PCA
Pre1_plot <- autoplot(Pre1_pca,
                      data = Pre1_samples,
                      colour = "group",
                      shape = "group",
                      size = 3) +
  theme_bw(base_size = 16, base_family = "Calibri") +
  scale_colour_manual("Treatment", values = c("#808080", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17)) +
  geom_text_repel(aes(label = animal), size = 4) +
  ggtitle("-1 wk")

Pre1_plot


############ One week post-infection
# Calculate PCs using log2CPM filtered counts
W1_pca <- prcomp(t(cpm(PPDb_norm[ , grep(pattern = "_W1_", x = colnames(PPDb_norm))],
                       log = TRUE)),
                 scale. = TRUE)

# Get sample names for plotting
W1_samples <- PPDb_norm$samples
W1_samples %<>%
  rownames_to_column(var = "sample_name") %>%
  dplyr::filter(stringr::str_detect(sample_name, "_W1_"))

# Plot PCA
W1_plot <- autoplot(W1_pca,
                    data = W1_samples,
                    colour = "group",
                    shape = "group",
                    size = 3) +
  theme_bw(base_size = 16, base_family = "Calibri") +
  scale_colour_manual("Treatment", values = c("#808080", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17)) +
  geom_text_repel(aes(label = animal), size = 4) +
  ggtitle("+1 wk")

W1_plot


############ Two weeks post-infection
# Calculate PCs using log2CPM filtered counts
W2_pca <- prcomp(t(cpm(PPDb_norm[ , grep(pattern = "_W2_", x = colnames(PPDb_norm))],
                       log = TRUE)),
                 scale. = TRUE)

# Get sample names for plotting
W2_samples <- PPDb_norm$samples
W2_samples %<>%
  rownames_to_column(var = "sample_name") %>%
  dplyr::filter(stringr::str_detect(sample_name, "_W2_"))

# Plot PCA
W2_plot <- autoplot(W2_pca,
                    data = W2_samples,
                    colour = "group",
                    shape = "group",
                    size = 3) +
  theme_bw(base_size = 16, base_family = "Calibri") +
  scale_colour_manual("Treatment", values = c("#808080", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17)) +
  geom_text_repel(aes(label = animal), size = 4) +
  ggtitle("+2 wk")

W2_plot


############ Ten weeks post-infection
# Calculate PCs using log2CPM filtered counts
W10_pca <- prcomp(t(cpm(PPDb_norm[ , grep(pattern = "_W10_", x = colnames(PPDb_norm))],
                        log = TRUE)),
                  scale. = TRUE)

# Get sample names for plotting
W10_samples <- PPDb_norm$samples
W10_samples %<>%
  rownames_to_column(var = "sample_name") %>%
  dplyr::filter(stringr::str_detect(sample_name, "_W10_"))

# Plot PCA
W10_plot <- autoplot(W10_pca,
                     data = W10_samples,
                     colour = "group",
                     shape = "group",
                     size = 3) +
  theme_bw(base_size = 16, base_family = "Calibri") +
  scale_colour_manual("Treatment", values = c("#808080", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17)) +
  geom_text_repel(aes(label = animal), size = 4) +
  ggtitle("+10 wk")

W10_plot


### Export high quality image for all MDS plots
files_PCA <- paste0(c("Pre1_plot", "W1_plot", "W2_plot", "W10_plot"), ".pdf")
plots_PCA <- list(Pre1_plot, W1_plot, W2_plot, W10_plot)

purrr::pwalk(list(files_PCA, plots_PCA),
      ggsave,
      device    = cairo_pdf,
      path      = imgDir,
      limitsize = FALSE,
      dpi       = 300,
      height    = 9,
      width     = 10,
      units     = "in")

#####################################################
# 19 Plot: Combine all PCA plots into single figure #
#####################################################

# Set grid
PCA_grid <- plot_grid((Pre1_plot + theme(legend.position = "none")),
                      (W1_plot  + theme(legend.position = "none")),
                      (W2_plot + theme(legend.position = "none")),
                      (W10_plot + theme(legend.position = "none")),
                      labels = c("A)", "B)", "C)", "D)"),
                      ncol = 2,
                      scale = .96)

# Check plot
PCA_grid

# Get legend from one MDS plot
legend <- get_legend(Pre1_plot)

# Add legend to grid
PCA_grid_leg <- plot_grid(PCA_grid, legend, ncol = 1,
                          rel_widths = c(1, .3),
                          rel_heights = c(1, .3))

# Check plot
PCA_grid_leg

# Export high quality image for both plots
files_PCAgrid <- paste0(c("PCA_grid", "PCA_grid_leg"), ".pdf")
plots_PCAgrid <- list(PCA_grid, PCA_grid_leg)

purrr::pwalk(list(files_PCAgrid, plots_PCAgrid),
             ggsave,
             device    = cairo_pdf,
             path      = imgDir,
             limitsize = FALSE,
             dpi       = 300,
             height    = 10,
             width     = 15,
             units     = "in")

##################################
# 20 Review experimental factors #
##################################

head(PPDb_norm$samples)

# Treatment (Non-stimulated or PPD-b stimulated)
levels(PPDb_norm$samples$group)

# Time points (infection status is nested within time points)
levels(PPDb_norm$samples$time.point)

# Animal ID
levels(PPDb_norm$samples$animal)

##############################
# 21 Create a design matrix  #
##############################

# This is based on Section 3.5 (Comparisons both between and within subjects)
# from the edgeR User's Guide (Last revised 12 May 2021)

# Initialize the design matrix with animal effects to adjust for baseline
# differences between the animals (effectively adding animal ID as
# a blocking factor due to paired nature of samples)
design <- model.matrix(~animal, data = PPDb_norm$samples)

# Now define the treatment-specific effects at each time point
# Pre1.PPDbStimulated: genes responding to PPDb at -1 wk pre-infection
Pre1.PPDbStimulated <- PPDb_norm$samples$time.point == "Pre1" &
                    PPDb_norm$samples$group == "PPDbStimulated"
# W1.PPDbStimulated: genes responding to PPDb at +1 wk post-infection
W1.PPDbStimulated <- PPDb_norm$samples$time.point == "W1" &
                    PPDb_norm$samples$group == "PPDbStimulated"
# W2.PPDbStimulated: genes responding to PPDb at +2 wk post-infection
W2.PPDbStimulated <- PPDb_norm$samples$time.point == "W2" &
                    PPDb_norm$samples$group == "PPDbStimulated"
# W10.PPDbStimulated: genes responding to PPDb at +10 wk post-infection
W10.PPDbStimulated <- PPDb_norm$samples$time.point == "W10" &
  PPDb_norm$samples$group == "PPDbStimulated"

# Append them to the design matrix
design <- cbind(design, Pre1.PPDbStimulated, W1.PPDbStimulated,
                W2.PPDbStimulated, W10.PPDbStimulated)

# Check design matrix
dim(design)
dim(PPDb_norm$samples)
head(design)

# Output the design matrix info
as.data.frame(design) %>%
  rownames_to_column(var = "rownames") %>%
  write_csv(file = here(tablesDir, "PPDb_design-matrix.csv"),
          col_names = TRUE)

#########################################
# 22 Estimate the dispersion parameters #
#########################################

# Common and trended dispersions are estimated with the
# Cox-Reid method and tagwise dispersions with the
# weighted likelihood empirical Bayes method to obtain posterior
# dispersion estimates
PPDb_disp <- estimateDisp(y       = PPDb_norm,
                          design  = design,
                          robust  = TRUE,
                          verbose = TRUE)

names(PPDb_disp)

# Check the calculated dispersion
PPDb_disp$common.dispersion

# Check the calculated dispersion's square root,
# which corresponds to the biological coefficient of variation (BCV)
sqrt(PPDb_disp$common.dispersion)

####################
# 23 Fit GLM model #
####################

# Fit a quasi-likelihood negative binomial generalized log-linear model
# for each tag using the design matrix and calculated dispersion
PPDb_fitQL <- glmQLFit(y = PPDb_disp,
                       design = design,
                       robust = TRUE)

names(PPDb_fitQL)
colnames(PPDb_fitQL$design)

##########################
# 24 Save R session info #
##########################

devtools::session_info()

#######################
# 25 Save .RData file #
#######################

# Detach all loaded packages (except base R packages):
require(nothing, quietly = TRUE)

# Save all environment objects:
save.image(file = "PPDb-RNA-seq_paired_sense.RData")

######################################
# Proceed to Part 3 of this analysis #
######################################

# File: 03-PPDb-RNA-seq_paired_sense.R








