######################################################################
# PPD-b stimulated vs unstimulated peripheral blood RNA-seq Analysis #
#                        --- heat maps ---                           #
######################################################################

# Author: Carolina N. Correia
# GitHub Repository DOI:
# Date: November 20th 2017

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
load("PPDb-RNA-seq-heatmaps.rda")

############################################
# 02 Load and/or install required packages #
############################################

# Load packages
library(tidyverse)
library(magrittr)
library(stringr)
library(reshape2)
library(Cairo)

##################
# 03 Import data #
##################

# Import DE genes
commonDE <- read_csv(file.path(paste0(tablesDir,
                                      "/common_DE_genes.csv")),
                     col_types = cols(EntrezID = "c"),
                     col_names = TRUE)

commonDEinfc <- read_csv(file.path(paste0(tablesDir,
                                          "/common_DE_infected_genes.csv")),
                         col_types = cols(EntrezID = "c"),
                         col_names = TRUE)

commonDEheal <- read_csv(file.path(paste0(tablesDir,
                                          "/healthy_DE_genes.csv")),
                         col_types = cols(EntrezID = "c"),
                         col_names = TRUE)

#############################
# 04 Tidy data for plotting #
#############################

# Subset logFC and FDR
commonDE %<>%
  dplyr::select(gene_symbol_Wm1, starts_with("logFC"),
                starts_with("FDR")) %>%
  dplyr::rename(gene_symbol = gene_symbol_Wm1)


commonDEinfc %<>%
  dplyr::select(gene_symbol_W1, starts_with("logFC"),
                starts_with("FDR")) %>%
  dplyr::rename(gene_symbol = gene_symbol_W1,
                logFC_W10   = logFC,
                FDR_W10     = FDR)

commonDEheal %<>%
  dplyr::select(gene_symbol, starts_with("logFC"),
                starts_with("FDR")) %>%
  dplyr::rename(logFC_Wm1 = logFC,
                FDR_Wm1   = FDR)

# Check dataframes
commonDE
commonDEinfc
commonDEheal

# Tidy dataframes with logFC data
commonDE %>%
  melt(id.vars = "gene_symbol") %>%
  dplyr::filter(stringr::str_detect(variable, "^log")) %>%
  dplyr::rename(logFC_timepoint = variable,
                log2FC_value    = value) -> logFC_DE
logFC_DE$gene_symbol %<>%
  as.factor()

commonDEinfc %>%
  melt(id.vars = "gene_symbol") %>%
  dplyr::filter(stringr::str_detect(variable, "^log")) %>%
  dplyr::rename(logFC_timepoint = variable,
                log2FC_value    = value) -> logFC_DEinfec
logFC_DEinfec$gene_symbol %<>%
  as.factor()

commonDEheal %>%
  melt(id.vars = "gene_symbol") %>%
  dplyr::filter(stringr::str_detect(variable, "^log")) %>%
  dplyr::rename(logFC_timepoint = variable,
                log2FC_value    = value) -> logFC_DEheal
logFC_DEheal$gene_symbol %<>%
  as.factor()


# Check tidy logFC dataframes
head(logFC_DE)
head(logFC_DEinfec)
head(logFC_DEheal)

##############################################################
# 05 Plot: heatmap of common DE genes across all time points #
##############################################################

# Palette colours
purpleGreen <- c("#00441B", "#1B7837", "#5AAE61",
                 "#C7E9C0", "white", "#F1B6DA",
                 "#9970AB", "#762A83", "#40004B")

# Order gene symbols alphabetically and reverse factor order
logFC_DE %<>% dplyr::arrange(gene_symbol)
logFC_DE$gene_symbol %<>% factor(levels = rev(levels(logFC_DE$gene_symbol)))

# Split data
split1 <- logFC_DE[1:208, ]
split2 <- logFC_DE[209:416, ]
split3 <- logFC_DE[417:nrow(logFC_DE), ]

# Plot first chunk
split1 %>%
  ggplot() +
  geom_tile(aes(x = logFC_timepoint, y = gene_symbol,
                fill = log2FC_value),
            colour = "black",
            size = 0.2) +
  coord_fixed() +
  scale_x_discrete(labels = c("logFC_Wm1" = "-1 wk",
                              "logFC_W1"  = "+1 wk",
                              "logFC_W2"  = "+2 wk",
                              "logFC_W10" = "+10 wk")) +
  scale_fill_gradientn(name = expression(paste(log[2], "FC")),
                       colours = purpleGreen,
                       limits = c(-10, 10),
                       breaks = c(-10, -8, -6, -4, -2, 0,
                                  2, 4, 6, 8, 10),
                       labels = c("-9", "-8", "-6", "-4", "-2", " 0",
                                  " 2", " 4", " 6", " 8", " 9")) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("DE genes across all time points") +
  theme(panel.border = element_rect(fill = NA,
                                    colour = "black", size = 1),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(face = "italic",
                                   hjust = 1),
        text = element_text(size = 14, family = "Calibri"),
        plot.title = element_text(hjust = 0.5)) -> heatmap_common1

# Plot second chunk
split2 %>%
  ggplot() +
  geom_tile(aes(x = logFC_timepoint, y = gene_symbol,
                fill = log2FC_value),
            colour = "black",
            size = 0.2) +
  coord_fixed() +
  scale_x_discrete(labels = c("logFC_Wm1" = "-1 wk",
                              "logFC_W1"  = "+1 wk",
                              "logFC_W2"  = "+2 wk",
                              "logFC_W10" = "+10 wk")) +
  scale_fill_gradientn(name = expression(paste(log[2], "FC")),
                       colours = purpleGreen,
                       limits = c(-10, 10),
                       breaks = c(-10, -8, -6, -4, -2, 0,
                                  2, 4, 6, 8, 10),
                       labels = c("-9", "-8", "-6", "-4", "-2", " 0",
                                  " 2", " 4", " 6", " 8", " 9")) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("DE genes across all time points") +
  theme(panel.border = element_rect(fill = NA,
                                    colour = "black", size = 1),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(face = "italic",
                                   hjust = 1),
        text = element_text(size = 14, family = "Calibri"),
        plot.title = element_text(hjust = 0.5)) -> heatmap_common2

# Plot third chunk
split3 %>%
  ggplot() +
  geom_tile(aes(x = logFC_timepoint, y = gene_symbol,
                fill = log2FC_value),
            colour = "black",
            size = 0.2) +
  coord_fixed() +
  scale_x_discrete(labels = c("logFC_Wm1" = "-1 wk",
                              "logFC_W1"  = "+1 wk",
                              "logFC_W2"  = "+2 wk",
                              "logFC_W10" = "+10 wk")) +
  scale_fill_gradientn(name = expression(paste(log[2], "FC")),
                       colours = purpleGreen,
                       limits = c(-10, 10),
                       breaks = c(-10, -8, -6, -4, -2, 0,
                                  2, 4, 6, 8, 10),
                       labels = c("-9", "-8", "-6", "-4", "-2", " 0",
                                  " 2", " 4", " 6", " 8", " 9")) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("DE genes across all time points") +
  theme(panel.border = element_rect(fill = NA,
                                    colour = "black", size = 1),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(face = "italic",
                                   hjust = 1),
        text = element_text(size = 14, family = "Calibri"),
        plot.title = element_text(hjust = 0.5)) -> heatmap_common3

### Export high quality image for all heatmap chunks
heat_files <- paste0(c("heatmap_common1", "heatmap_common2", "heatmap_common3"),
                     ".pdf")
heat_plots <- list(heatmap_common1, heatmap_common2, heatmap_common3)

purrr::pwalk(list(heat_files, heat_plots),
             ggsave,
             device    = cairo_pdf,
             path      = imgDir,
             limitsize = FALSE,
             dpi       = 300,
             height    = 10,
             width     = 8,
             units     = "in")

###############################################################
# 06 Plot: heatmap of common DE genes in infected time points #
###############################################################

# Order gene symbols alphabetically and reverse factor order
logFC_DEinfec %<>% dplyr::arrange(gene_symbol)
logFC_DEinfec$gene_symbol %<>%
  factor(levels = rev(levels(logFC_DEinfec$gene_symbol)))

# Split data
spl1 <- logFC_DEinfec[1:324, ]
spl2 <- logFC_DEinfec[325:645, ]
spl3 <- logFC_DEinfec[646:966, ]
spl4 <- logFC_DEinfec[967:nrow(logFC_DEinfec), ]

# Plot first chunk
spl1 %>%
  ggplot() +
  geom_tile(aes(x = logFC_timepoint, y = gene_symbol,
                fill = log2FC_value),
            colour = "black",
            size = 0.2) +
  coord_fixed() +
  scale_x_discrete(labels = c("logFC_W1"  = "+1 wk",
                              "logFC_W2"  = "+2 wk",
                              "logFC_W10" = "+10 wk")) +
  scale_fill_gradientn(name = expression(paste(log[2], "FC")),
                       colours = purpleGreen,
                       limits = c(-8, 9),
                       breaks = c(-8, -6, -4, -2, 0,
                                  2, 4, 6, 8.5),
                       labels = c("-8", "-6", "-4", "-2", " 0",
                                  " 2", " 4", " 6", " 8.5")) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("DE genes across infected time points") +
  theme(panel.border = element_rect(fill = NA,
                                    colour = "black", size = 1),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(face = "italic",
                                   hjust = 1),
        text = element_text(size = 14, family = "Calibri"),
        plot.title = element_text(hjust = 0.5)) -> heatmap_infec1

# Plot second chunk
spl2 %>%
  ggplot() +
  geom_tile(aes(x = logFC_timepoint, y = gene_symbol,
                fill = log2FC_value),
            colour = "black",
            size = 0.2) +
  coord_fixed() +
  scale_x_discrete(labels = c("logFC_W1"  = "+1 wk",
                              "logFC_W2"  = "+2 wk",
                              "logFC_W10" = "+10 wk")) +
  scale_fill_gradientn(name = expression(paste(log[2], "FC")),
                       colours = purpleGreen,
                       limits = c(-8, 9),
                       breaks = c(-8, -6, -4, -2, 0,
                                  2, 4, 6, 8.5),
                       labels = c("-8", "-6", "-4", "-2", " 0",
                                  " 2", " 4", " 6", " 8.5")) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("DE genes across infected time points") +
  theme(panel.border = element_rect(fill = NA,
                                    colour = "black", size = 1),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(face = "italic",
                                   hjust = 1),
        text = element_text(size = 14, family = "Calibri"),
        plot.title = element_text(hjust = 0.5)) -> heatmap_infec2

# Plot third chunk
spl3 %>%
  ggplot() +
  geom_tile(aes(x = logFC_timepoint, y = gene_symbol,
                fill = log2FC_value),
            colour = "black",
            size = 0.2) +
  coord_fixed() +
  scale_x_discrete(labels = c("logFC_W1"  = "+1 wk",
                              "logFC_W2"  = "+2 wk",
                              "logFC_W10" = "+10 wk")) +
  scale_fill_gradientn(name = expression(paste(log[2], "FC")),
                       colours = purpleGreen,
                       limits = c(-8, 9),
                       breaks = c(-8, -6, -4, -2, 0,
                                  2, 4, 6, 8.5),
                       labels = c("-8", "-6", "-4", "-2", " 0",
                                  " 2", " 4", " 6", " 8.5")) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("DE genes across infected time points") +
  theme(panel.border = element_rect(fill = NA,
                                    colour = "black", size = 1),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(face = "italic",
                                   hjust = 1),
        text = element_text(size = 14, family = "Calibri"),
        plot.title = element_text(hjust = 0.5)) -> heatmap_infec3

# Plot fourth chunk
spl4 %>%
  ggplot() +
  geom_tile(aes(x = logFC_timepoint, y = gene_symbol,
                fill = log2FC_value),
            colour = "black",
            size = 0.2) +
  coord_fixed() +
  scale_x_discrete(labels = c("logFC_W1"  = "+1 wk",
                              "logFC_W2"  = "+2 wk",
                              "logFC_W10" = "+10 wk")) +
  scale_fill_gradientn(name = expression(paste(log[2], "FC")),
                       colours = purpleGreen,
                       limits = c(-8, 9),
                       breaks = c(-8, -6, -4, -2, 0,
                                  2, 4, 6, 8.5),
                       labels = c("-8", "-6", "-4", "-2", " 0",
                                  " 2", " 4", " 6", " 8.5")) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("DE genes across infected time points") +
  theme(panel.border = element_rect(fill = NA,
                                    colour = "black", size = 1),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(face = "italic",
                                   hjust = 1),
        text = element_text(size = 14, family = "Calibri"),
        plot.title = element_text(hjust = 0.5)) -> heatmap_infec4

### Export high quality image for all heatmap chunks
heat_files2 <- paste0(c("heatmap_infec1", "heatmap_infec2",
                        "heatmap_infec3", "heatmap_infec4"),
                     ".pdf")
heat_plots2 <- list(heatmap_infec1, heatmap_infec2,
                    heatmap_infec3, heatmap_infec4)

purrr::pwalk(list(heat_files2, heat_plots2),
             ggsave,
             device    = cairo_pdf,
             path      = imgDir,
             limitsize = FALSE,
             dpi       = 300,
             height    = 13,
             width     = 8,
             units     = "in")

########################################################
# 07 Plot: heatmap of DE genes in healthy samples only #
########################################################

# Palette colours
purpleGreen2 <- c("#00441B", "#1B7837", "white", "#762A83", "#40004B")

logFC_DEheal %>%
    ggplot() +
    geom_tile(aes(x = logFC_timepoint, y = gene_symbol,
                  fill = log2FC_value),
              colour = "black",
              size = 0.2) +
    coord_fixed() +
    scale_y_discrete(limits = rev(levels(logFC_DEheal$gene_symbol))) +
    scale_x_discrete(labels = c("logFC_Wm1" = "-1 wk")) +
    scale_fill_gradientn(name = expression(paste(log[2], "FC")),
                         colours = purpleGreen2,
                         limits = c(-1, 2),
                         breaks = c(-1, 0, 1, 2),
                         labels = c("-1", " 0", " 1", " 2")) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle("DE genes unique to non-infected samples") +
    theme(panel.border = element_rect(fill = NA,
                                      colour = "black", size = 1),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.text.y = element_text(face = "italic",
                                     hjust = 1),
          text = element_text(size = 14, family = "Calibri"),
          plot.title = element_text(hjust = 0.5)) -> heatmap_heal

heatmap_heal

ggsave("heatmap_healthy_DE.pdf",
       plot      = heatmap_heal,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 8,
       units     = "in")

############################
# 08 Save Save .RData file #
############################

save.image("PPDb-RNA-seq-heatmaps.rda")

##########################
# 09 Save R session info #
##########################

devtools::session_info()

#######
# END #
#######