######################################################################
# PPD-b stimulated vs unstimulated peripheral blood RNA-seq Analysis #
#                    --- RNA-seq statistics ---                      #
######################################################################

# Author: Carolina N. Correia
# GitHub Repository DOI:
# Date: October 7th 2017

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/PPD-B-sti_vs_uns/edgeR_sense")

# Define variables for specific directories
workDir <- getwd()
tablesDir <- paste0(workDir, "/Tables")
imgDir <- paste0(workDir, "/Figures")

############################################
# 02 Load and/or install required packages #
############################################

# Load packages
library(plyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(devtools)
library(waffle)
library(extrafont)
library(Cairo)
library(skimr)
library(scales)
library(cowplot)

###################################
# 03 Import ngsShoRT summary file #
###################################

# Import file into a tibble
ngsShortAll <- read_table2("ngsshort_cattle.txt")

View(ngsShortAll)

#######################
# 04 Clean data frame #
#######################

# Remove unnecessary columns
ngsShortAll %<>%
  dplyr::select(-contains("X"))

# Remove fresh samples
ngsShortPU <- dplyr::filter(ngsShortAll, !grepl("_F_", Sample_name))

# Check data frame
View(ngsShortPU)

########################
# 05 Add labels column #
########################

# Add new column for plotting labels
ngsShortPU$labels <- ngsShortPU$Sample_name

# Correct labels
ngsShortPU$labels %<>%
  str_replace("A", "") %>%
  str_replace("_00\\d", "") %>%
  factor()

# Check data frame
ngsShortPU

###########################
# 06 Summarise lane stats #
###########################

ngsShortPU %<>%
  dplyr::group_by(labels) %>%
  dplyr::summarise_at(vars(Read_pair_count:Removed_PE_pair_count), sum)

# Check data frame
ngsShortPU

################################################################
# 07 Add animal ID, time point, and experimental group columns #
################################################################

# Animal ID
ngsShortPU$Animal_identifier <- ngsShortPU$labels
ngsShortPU$Animal_identifier %<>%
  str_replace("_W.+_\\w", "") %>%
  factor()

# Experimental group
ngsShortPU$Experimental_group <- ngsShortPU$labels
ngsShortPU$Experimental_group %<>%
  str_replace("\\d\\d\\d\\d_.+_", "") %>%
  str_replace("P", "PPD-b") %>%
  str_replace("U", "Control") %>%
  factor(levels = c("PPD-b", "Control"))

# Time point
ngsShortPU$Infection_time_point <- ngsShortPU$labels
ngsShortPU$Infection_time_point %<>%
  str_replace("\\d\\d\\d\\d", "") %>%
  str_replace("(P|U)", "") %>%
  str_replace("_W-1_", "-1 wk") %>%
  str_replace("_W1_", "+1 wk") %>%
  str_replace("_W2_", "+2 wk") %>%
  str_replace("_W10_", "+10 wk") %>%
  factor(levels = c("-1 wk", "+1 wk", "+2 wk", "+10 wk"))

# Check data frame
ngsShortPU

#######################
# 08 Order data frame #
#######################

# Order rows
ngsShortPU %<>% dplyr::arrange(Animal_identifier, Experimental_group)

# Rename and order columns
ngsShortPU %<>%
  dplyr::rename(`No. reads pre-filtering` = Read_pair_count) %>%
  dplyr::select(Animal_identifier,
                Experimental_group,
                Infection_time_point,
                everything())

# Check data frame
ngsShortPU

##############################
# 09 Add more ngsShort stats #
##############################

# Number and proportion of input reads
ngsShortPU %<>%
  dplyr::mutate(`No. reads post-filtering (input)` = `No. reads pre-filtering`
                                                      - Removed_PE_pair_count,
                `% reads post-filtering` = `No. reads post-filtering (input)`
                                            / `No. reads pre-filtering` * 100) %>%
  dplyr::select(-Removed_PE_pair_count)

# Check data frame
View(ngsShortPU)

###############################
# 10 Import STAR summary file #
###############################

# Import file and convert to tibble
STARlogAll <- as.tibble(read.delim("All_star_log_final_out.txt"))

View(STARlogAll)

#######################
# 11 Clean data frame #
#######################

# Remove unnecessary columns
STARlogAll %<>%
  dplyr::select(-Number.of.input.reads,
                -Number.of.reads.mapped.to.too.many.loci,
                -X..of.reads.mapped.to.too.many.loci)

# Remove fresh samples
STARlogAll <- dplyr::filter(STARlogAll, !grepl("_F", Sample_id))

# Check data frame
View(STARlogAll)

################################
# 12 Rename data frame columns #
################################

STARlogAll %<>%
  dplyr::rename(labels = Sample_id,
                `No. uniquely mapped reads` = Uniquely.mapped.reads.number,
                `% uniquely mapped reads` = Uniquely.mapped.reads..,
                `Mean mapped length` = Average.mapped.length,
                `No. reads mapped to multiple loci` = Number.of.reads.mapped.to.multiple.loci,
                `% reads mapped to multiple loci` = X..of.reads.mapped.to.multiple.loci,
                `No. unmapped reads` = Unmapped.reads,
                `% unmapped reads` = X..Unmapped.reads)

# Check data frame
View(STARlogAll)

#########################################
# 13 Join ngsShort and STAR data frames #
#########################################

# Clean labels
STARlogAll$labels %<>%
  str_replace("A", "") %>%
  factor(levels = levels(ngsShortPU$labels))

# Join data frames
Supp_tableIV2 <- inner_join(ngsShortPU, STARlogAll)

# Check data frame
View(Supp_tableIV2)

###########################
# 14 Export RNA-seq stats #
###########################

write_csv(Supp_tableIV2,
          path = file.path(paste0(tablesDir, "/Supp_tableIV2.csv")),
          col_names = TRUE)

######################################
# 15 Calculate RNA-seq summary stats #
######################################

# Remove percent symbol and convert data to numeric type
Supp_tableIV2$`% uniquely mapped reads` %<>%
  str_replace("%", "") %>%
  as.numeric()

Supp_tableIV2$`% reads mapped to multiple loci` %<>%
  str_replace("%", "") %>%
  as.numeric()

# Calculate summary stats
RNAseqStats <- skim(Supp_tableIV2)

# Tidy summary stats
stats_to_keep <- c("mean", "sd", "max", "min", "median")

RNAseqStats %>%
  dplyr::filter(stat %in% stats_to_keep) %>%
  tidyr::spread(stat, value) %>%
  dplyr::select(-type, -level) -> RNAseqStats_tidy

###################
# 16 Set up fonts #
###################

# Import fonts
font_import()

# Registered fonts with R for the PDF output device
loadfonts()

###########################################
# 17 Plot: Waffle charts of RNA-seq stats #
###########################################

# Filtering
RNAseqStats_tidy %>%
  dplyr::filter(var == "% reads post-filtering") %>%
  dplyr::select(var, mean) %>%
  tidyr::spread(var, mean) %>%
  dplyr::mutate(`% removed reads` = 100 - `% reads post-filtering`) %>%
  round()

filt_reads <- c(`Reads post-filtering (82%)` = 82,
                `Removed reads (18%)` = 18)

filt_waffle <- waffle(filt_reads,
       rows = 5,
       size = 0.5,
       title = "Filtering of reads",
       colors = c("#5ab4ac", "#d8b365"))

filt_waffle <- filt_waffle + theme(text=element_text(size = 16, family = "Calibri"))


# Alignment
vars1 <- c("% uniquely mapped reads",
           "% reads mapped to multiple loci",
           "% unmapped reads")

RNAseqStats_tidy %>%
  dplyr::filter(var %in% vars1) %>%
  dplyr::select(var, mean) %>%
  tidyr::spread(var, mean) %>%
  round() %>%
  View()


align_reads <- c(`Uniquely mapped reads (85%)` = 85,
                 `Unmapped reads (11%)` = 11,
                 `Reads mapped to multiple loci (4%)` = 4)
align_waffle <- waffle(align_reads,
                       rows = 5,
                       size = 0.5,
                       title = "Alignment of reads",
                       colors = c("#5e3c99", "#b2abd2", "#e66101"))

align_waffle <- align_waffle + theme(text=element_text(size = 16, family = "Calibri"))

# Combine plots
grid_waffle <- plot_grid(filt_waffle,
                         align_waffle,
                         labels = c("A", "B"),
                         nrow = 2)

# Save combined plot
ggsave("grid_waffle.pdf",
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 9,
       width     = 16,
       units     = "in")

##########################
# 17 Save R session info # ----
##########################

devtools::session_info()

#######
# END #
#######