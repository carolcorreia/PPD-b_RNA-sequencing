######################################################################
# PPD-b stimulated vs unstimulated peripheral blood RNA-seq Analysis #
#                    --- RNA-seq statistics ---                      #
######################################################################

# Author: Carolina N. Correia
# Last updated on 28/07/2021

############################################
# 01 Load and/or install required packages #
############################################

# Load packages
library(tidyverse)
library(magrittr)
library(devtools)
library(waffle)
library(extrafont)
library(Cairo)
library(skimr)
library(scales)
library(cowplot)
library(here)

# Uncomment functions below to install packages in case you don't have them

# GitHub package
#devtools::install_github("hrbrmstr/waffle")

# CRAN package
#install.packages("skimr")
#install.packages("scales")
#install.packages("cowplot")

##############################################
# 02 Working directory, fonts, and time zone #
##############################################

# Check working directory
here()

# Define variables for subdirectories
imgDir <- here("Figures")
tablesDir <- here("Tables")

# Set time zone
Sys.setenv(TZ = "Europe/Dublin")

# Register fonts with R for the PDF output device
font_import()

# Load fonts
loadfonts()

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
STARlogAll <- read_delim("All_star_log_final_out.txt", delim = "\t")

View(STARlogAll)

#######################
# 11 Clean data frame #
#######################

# Remove unnecessary column
STARlogAll %<>%
  dplyr::select(-`Number of input reads`)

# Check data frame
View(STARlogAll)

#########################################
# 12 Join ngsShort and STAR data frames #
#########################################

# Create column with keys for joining data frames
STARlogAll$labels <- STARlogAll$Sample_id

# Clean labels
STARlogAll$labels %<>%
  str_replace("A", "") %>%
  factor(levels = levels(ngsShortPU$labels))

# Join data frames
Supp_table2 <- inner_join(ngsShortPU, STARlogAll)

# Check data frame
View(Supp_table2)

###########################
# 13 Export RNA-seq stats #
###########################

write_csv(Supp_table2,
          file = here(tablesDir, "Supp_table2.csv"),
          col_names = TRUE)

######################################
# 14 Calculate RNA-seq summary stats #
######################################

# Remove percent symbol and convert data to numeric type
Supp_table2$`Uniquely mapped reads %` %<>%
  str_replace("%", "") %>%
  as.numeric()

Supp_table2$`% of reads mapped to multiple loci` %<>%
  str_replace("%", "") %>%
  as.numeric()

Supp_table2$`% of reads mapped to too many loci` %<>%
  str_replace("%", "") %>%
  as.numeric()

# Calculate summary stats
RNAseqStats <- skim(Supp_table2)

###########################################
# 15 Plot: Waffle charts of RNA-seq stats #
###########################################

# Get the mean values for read filtering
RNAseqStats %>%
  dplyr::filter(skim_variable == "% reads post-filtering")

filt_reads <- c(`Reads post-filtering (82%)` = 82,
                `Removed reads (18%)` = 18)

filt_waffle <- waffle(filt_reads,
       rows = 5,
       size = 0.5,
       title = "Filtering of reads",
       colors = c("#5ab4ac", "#d8b365"))

filt_waffle <- filt_waffle +
  theme(text=element_text(size = 16, family = "Calibri"))


# Get the mean values for read alignment
RNAseqStats %>%
  dplyr::filter(skim_variable == "Uniquely mapped reads %")

RNAseqStats %>%
  dplyr::filter(skim_variable == "% Unmapped reads")

RNAseqStats %>%
  dplyr::filter(skim_variable == "% of reads mapped to multiple loci")

RNAseqStats %>%
  dplyr::filter(skim_variable == "% of reads mapped to too many loci")

align_reads <- c(`Uniquely mapped reads (87.3%)` = 87.3,
                 `Unmapped reads (9.5%)` = 9.5,
                 `Reads mapped to multiple or too many loci (3.2%)` = 3.2)
align_waffle <- waffle(align_reads,
                       rows = 5,
                       size = 0.5,
                       title = "Alignment of reads",
                       colors = c("#5e3c99", "#b2abd2", "#e66101"))

align_waffle <- align_waffle +
  theme(text=element_text(size = 16, family = "Calibri"))

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
# 17 Save R session info #
##########################

devtools::session_info()

