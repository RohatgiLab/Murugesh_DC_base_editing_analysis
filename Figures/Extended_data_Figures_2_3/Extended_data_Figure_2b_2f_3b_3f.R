# zscore vs protein length plots: Extended Figure 2b, 2f, 3b, 3f:
# Please use supplementary data table 01 for making these figures:

library(readxl)
library(openxlsx)
library(janitor)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Protein_length vs z-score plots :

dc_cbe <- read_excel("Z_score_vs_protein_length_plots/CBE_bot_zscore_plots/dc_cbe_for_zscore_plots.xlsx")

# Function I used to categorize mutation_category: 

categorize_mutation_cbe <- function(data_column) {
  data_column <- str_to_lower(data_column)  
  
  
  case_when(
    str_detect(data_column, "nonsense") ~ "Nonsense",  
    str_detect(data_column, "intron|splice-acceptor|splice-donor|utr") ~ "Intron", 
    str_detect(data_column, "missense") ~ "Missense",  
    str_detect(data_column, "silent") ~ "Silent",  
    TRUE ~ NA  
  )
}

custom_colors <- c(
  "Missense" = "#5EB0B3",  
  "Nonsense"      = "#FA6C44",  
  "Splice-site"    = "#FDBA47",
  "Intron" = "#FDBA47",
  "UTR" = "#FDBA47",
  "Silent"   = "#F8F8F6"
  
)


fade_color <- rgb(245/255, 245/255, 245/255, alpha = 0.005)


# beta_catenin z_score vs protein_length: Extended Figure 2b:

ids_to_annotate <- dc_cbe[dc_cbe$z_score_bot >= 3.0 & dc_cbe$gene_names == "CTNNB1", ]

ctnnb1_cbe_bot_with_labels <- ggplot(
  data = dc_cbe[dc_cbe$gene_names == "CTNNB1",],
  aes(
    x = residue_number, y = z_score_bot, 
    fill = ifelse(z_score_bot > -2.0 & z_score_bot < 2.0, fade_color, categorized_mutation),
    alpha = ifelse(z_score_bot > -2.0 & z_score_bot < 2.0, 0.01, 1)
  )
) + 
  geom_point(size = 2.0, shape = 21, stroke = 0.175) + 
  scale_fill_manual(values = c(custom_colors, "white"="white")) +
  scale_alpha_identity() + 
  theme_pubr() + 
  theme(legend.position = "none") +
  geom_hline(yintercept = c(-2.0, 2.0), linetype = "dotted", color = "black") + 
  geom_text_repel(data = ids_to_annotate, aes(label = ids),
                  color = "black",
                  size = 1.5, 
                  box.padding = 0.0, 
                  point.padding =0.2,
                  segment.length = 5, 
                  segment.size = 0.2,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  nudge_x = 0.5, 
                  nudge_y = 2.0) + 
  theme(axis.ticks.length = unit(1.5, "pt")) + 
  scale_y_continuous(limits = c(-10, 12)) +
  scale_x_continuous(limits = c(0, 790),expand = expansion(mult = c(0.05, 0)), breaks = c(seq(0, 781, 200), 781))

ggsave(filename = "Z_score_vs_protein_length_plots/CBE_bot_zscore_plots/beta_catenin_cbe_bot_with_labels.pdf",
       plot = ctnnb1_cbe_bot_with_labels,
       dpi = 3000,
       width = 4.21,
       height = 1.46,
       units = "in")

# Axin1 zscore_plot vs protein_length: Extended Figure 2f:

ids_to_annotate <- dc_cbe[dc_cbe$z_score_bot >= 3.0 & dc_cbe$gene_names == "Axin1", ]

Axin1_cbe_bot_with_names <- ggplot(
  data = dc_cbe[dc_cbe$gene_names == "Axin1",],
  aes(
    x = residue_number, y = z_score_bot, 
    fill = ifelse(z_score_bot > -2.0 & z_score_bot < 2.0, fade_color, categorized_mutation),
    alpha = ifelse(z_score_bot > -2.0 & z_score_bot < 2.0, 0.01, 1)
  )
) + 
  geom_point(size = 2.0, shape = 21, stroke = 0.175) + 
  scale_fill_manual(values = c(custom_colors, "white"="white")) +
  scale_alpha_identity() + 
  theme_pubr() + 
  theme(legend.position = "none") +
  geom_hline(yintercept = c(-2.0, 2.0), linetype = "dotted", color = "black") + 
  geom_text_repel(data = ids_to_annotate, aes(label = ids),
                  color = "black",
                  size = 1.5, 
                  box.padding = 0.0, 
                  point.padding = 0.2,
                  segment.length = 5, 
                  segment.size = 0.2,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  nudge_x = 0.5, 
                  nudge_y = 2.0) + 
  scale_y_continuous(limits = c(-10, 15)) +
  theme(axis.ticks.length = unit(1.5, "pt")) +
  scale_x_continuous(limits = c(0, 870), expand =expansion(mult = c(0.05, 0)), breaks = c(seq(0, 870, 200), 862))


ggsave(filename = "Z_score_vs_protein_length_plots/CBE_bot_zscore_plots/Axin1_cbe_bot_with_labels.pdf",
       plot = Axin1_cbe_bot_with_names,
       dpi = 3000,
       width = 4.56,
       height = 1.46,
       units = "in")

# APC zscore vs protein_length: Extended Figure 3b:

ids_to_annotate <- dc_cbe[dc_cbe$z_score_bot >= 3.0 & dc_cbe$gene_names == "APC", ]

APC_cbe_bot_with_names <- ggplot(
  data = dc_cbe[dc_cbe$gene_names == "APC",],
  aes(
    x = residue_number, y = z_score_bot, 
    fill = ifelse(z_score_bot > -2.0 & z_score_bot < 2.0, fade_color, categorized_mutation),
    alpha = ifelse(z_score_bot > -2.0 & z_score_bot < 2.0, 0.008, 1)
  )
) + 
  geom_point(size = 2, shape = 21, stroke = 0.175) + 
  scale_fill_manual(values = c(custom_colors, "white"="white")) +
  scale_alpha_identity() + 
  theme_pubr() + 
  theme(legend.position = "none") +
  geom_hline(yintercept = c(-2.0, 2.0), linetype = "dotted", color = "black") + 
  geom_text_repel(data = ids_to_annotate, aes(label = ids),
                  color = "black",
                  size =1.0, 
                  box.padding = 0.0, 
                  point.padding = 0.2,
                  segment.length = 5, 
                  segment.size = 0.5,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  nudge_x = 0.5, 
                  nudge_y = 2.0) +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_x_continuous(limits = c(0, 2850), expand =expansion(mult = c(0.05, 0.01)), breaks = c(seq(0, 2850, 500), 2843)) +
  theme(axis.ticks.length = unit(1.5, "pt"))


ggsave(filename = "Z_score_vs_protein_length_plots/CBE_bot_zscore_plots/APC_cbe_bot_with_names.pdf",
       plot = APC_cbe_bot_with_names,
       dpi = 3000,
       width = 6.2,
       height = 1.46, 
       units = "in")


# GSK3B z_score vs protein_length: Extended Figure 3f:


ids_to_annotate <- dc_cbe[dc_cbe$z_score_bot >= 3.0 & dc_cbe$gene_names == "Gsk3b", ]

Gsk3b_cbe_bot_with_names <- ggplot(
  data = dc_cbe[dc_cbe$gene_names == "Gsk3b",],
  aes(
    x = residue_number, y = z_score_bot, 
    fill = ifelse(z_score_bot > -2.0 & z_score_bot < 2.0, fade_color, categorized_mutation),
    alpha = ifelse(z_score_bot > -2.0 & z_score_bot < 2.0, 0.008, 1)
  )
) + 
  geom_point(size = 2, shape = 21, stroke = 0.175) + 
  scale_fill_manual(values = c(custom_colors, "white"="white")) +
  scale_alpha_identity() + 
  theme_pubr() + 
  theme(legend.position = "none") +
  geom_hline(yintercept = c(-2.0, 2.0), linetype = "dotted", color = "black") + 
  geom_text_repel(data = ids_to_annotate, aes(label = ids),
                  color = "black",
                  size = 3, 
                  box.padding = 0.0, 
                  point.padding = 0.2,
                  segment.length = 5,
                  segment.size = 0.5,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  nudge_x = 0.5, 
                  nudge_y = 2.0) + 
  scale_y_continuous(limits = c(-8, 8)) +
  scale_x_continuous(limits = c(0, 430), expand =expansion(mult = c(0.05, 0.01)), breaks = c(seq(0, 430, 200), 420)) +
  theme(axis.ticks.length = unit(1.5, "pt"))



ggsave(filename = "Z_score_vs_protein_length_plots/CBE_bot_zscore_plots/Gsk3b_cbe_bot_with_names.pdf",
       plot = Gsk3b_cbe_bot_with_names,
       dpi = 3000,
       width = 2.05,
       height = 1.46,
       units = "in")

