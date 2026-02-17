# False negative rate: Extended data Figure 1e, 
library(readxl)
library(ggridges)
library(ggplot2)
library(ggridges)
library(openxlsx)
library(janitor)
library(ggpubr)
library(data.table)
library(stringr)

# Extended Data Figure 1e: 

dc_abe_data <- read_excel("False_Negative_rate/dc_abe_for_z_score_plots.xlsx") |> data.table()

ctnnb1_abe <- dc_abe_data |> filter(
  gene_names %in% c("CTNNB1", "control_abe")
)

ctnnb1_abe$categorized_mutation <- factor(
  ctnnb1_abe$categorized_mutation, 
  levels = c("non_targeting","Intron","Missense")
)

abe_outcomes <- ggplot(
  ctnnb1_abe,
  aes(x = z_score_bot, y = categorized_mutation, fill = categorized_mutation)
) +geom_density_ridges(scale = 0.95) +
  theme_ridges() +
  scale_x_continuous(breaks = seq(-6, 12, by = 2), limits = c(-6, 12)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.2))) + theme_pubr() +
  theme(legend.position = "none")

abe_outcomes
ggsave(
  filename = "Categories of mutations/control_abe_bot.pdf",
  plot = abe_outcomes,
  height = 10, 
  width = 10, 
  units = "in",
  dpi = 600
) 

# To calculate True positives: 
table(ctnnb1_abe$categorized_mutation)
nrow(ctnnb1_abe |> filter(gene_names == "CTNNB1" & categorized_mutation == "Missense" & z_score_bot >= 2.0))
nrow(ctnnb1_abe |> filter(gene_names == "CTNNB1" & categorized_mutation == "Intron" & z_score_bot >= 2.0))

#  Beta-Catenin from the CBE bottom gate: Extended Data Figure 1e

dc_cbe_data <-read_excel("False_Negative_rate/dc_cbe_for_zscore_plots.xlsx") |> data.table()
ctnnb1_cbe <- dc_cbe_data |> filter(
  gene_names %in% c("CTNNB1", "control")
)
table(ctnnb1_cbe$categorized_mutation)

ctnnb1_cbe$categorized_mutation <- factor(
  ctnnb1_cbe$categorized_mutation, 
  levels = c("non_targeting", "Missense","Intron","Nonsense")
)

nrow(ctnnb1_cbe |> filter(gene_names == "CTNNB1" & categorized_mutation == "Intron" & z_score_bot >= 2.0 ))

cbe_outcomes <- ggplot(
  ctnnb1_cbe,
  aes(x = z_score_bot, y = categorized_mutation, fill = categorized_mutation)
) +geom_density_ridges(scale = 0.95) +
  theme_ridges() +
  scale_x_continuous(breaks = seq(-6, 12, by = 2), limits = c(-6, 12)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.2))) + theme_pubr() +
  theme(legend.position = "none")

cbe_outcomes
ggsave(
  filename = "Categories of mutations/control_cbe_bot.pdf",
  plot = cbe_outcomes,
  height = 10, 
  width = 10, 
  units = "in",
  dpi = 600
)

# To calculate true positives: 
table(ctnnb1_cbe$categorized_mutation)
nrow(ctnnb1_cbe |> filter(gene_names == "CTNNB1" & categorized_mutation == "Nonsense" & z_score_bot >= 2.0))
nrow(ctnnb1_cbe |> filter(gene_names == "CTNNB1" & categorized_mutation == "Missense" & z_score_bot >= 2.0))
nrow(ctnnb1_cbe |> filter(gene_names == "CTNNB1" & categorized_mutation == "Intron" & z_score_bot >= 2.0))

# Extended data figure 1F: 

# APC-ABE top-gate: 

apc_abe <- dc_abe_data |> filter(
  gene_names %in% c("APC", "control_abe")
)

apc_abe <- apc_abe |> filter(
  !str_detect(predicted_mutations, "NC")
)

apc_abe$categorized_mutation <- factor(
  apc_abe$categorized_mutation, 
  levels = c("non_targeting","Intron","Missense")
)


table(apc_abe$categorized_mutation)
table(apc_abe$categorized_mutation != "non_targeting")
nrow(apc_abe |> filter(gene_names == "control_abe" & categorized_mutation == "non_targeting" & z_score_top >= 2.0))
nrow(apc_abe |> filter(gene_names == "APC" & categorized_mutation == "Intron" & residue_number <= 1570))


nrow(apc_abe |> filter(gene_names == "APC" & categorized_mutation == "Intron" & z_score_top >= 2.0 & residue_number <= 1570))

apc_abe_plot <- ggplot(
  apc_abe,
  aes(x = z_score_top, y = categorized_mutation, fill = categorized_mutation)
) +geom_density_ridges(scale = 0.95) +
  theme_ridges() +
  scale_x_continuous(breaks = seq(-6, 12, by = 2), limits = c(-6, 12)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.2))) + theme_pubr() +
  theme(legend.position = "none")

apc_abe_plot
ggsave(
  filename = "Categories of mutations/control_abe_top.pdf",
  plot = apc_abe_plot,
  height = 10, 
  width = 10, 
  units = "in",
  dpi = 600
)

# To calculate true positives: 
table(apc_abe$categorized_mutation)
nrow(apc_abe |> filter(gene_names == "APC" & categorized_mutation == "Intron" & z_score_top >= 2.0))

# Extended data figure 1F: APC_CBE top gate: 

apc_cbe <- dc_cbe_data |> filter(
  gene_names %in% c("APC", "control")
)

table(apc_cbe$categorized_mutation)

apc_cbe$categorized_mutation <- factor(
apc_cbe$categorized_mutation, 
  levels = c("Nonsense", "Missense","Intron","non_targeting")
)

table(apc_cbe$categorized_mutation == "Intron")

nrow(apc_cbe |> filter(gene_names == "APC" & residue_number <= 1570 & categorized_mutation == "Intron"))

nrow(apc_cbe |> filter(gene_names == "APC" & categorized_mutation == "Intron" & z_score_top >= 2.0 & residue_number <= 1570 ))

apc_cbe_outcomes <- ggplot(
  apc_cbe,
  aes(x = z_score_top, y = categorized_mutation, fill = categorized_mutation)
) +geom_density_ridges(scale = 0.95) +
  theme_ridges() +
  scale_x_continuous(breaks = seq(-6, 12, by = 2), limits = c(-6, 12)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.2))) + theme_pubr() +
  theme(legend.position = "none")

apc_cbe_outcomes
ggsave(
  filename = "Categories of mutations/apc_top_cbe.pdf",
  plot = apc_cbe_outcomes,
  height = 10, 
  width = 10, 
  units = "in",
  dpi = 600
)

# To calculate outcomes: To calculate Nonsense true positives, we only considered 
# APC residues (1-1570) because nonsense mutation after SAMP1 repeat has no phenotype (i.e doesn't elevate the pathway);

table(apc_cbe$categorized_mutation)

nrow(apc_cbe |> filter(gene_names == "APC" &  categorized_mutation == "Intron"))
nrow(apc_cbe |> filter(gene_names == "APC" & categorized_mutation == "Intron" & z_score_top >= 2.0 ))
nrow(apc_cbe |> filter(gene_names == "APC" & categorized_mutation == "Nonsense" & residue_number <= 1570))
nrow(apc_cbe |> filter(gene_names == "APC" & categorized_mutation == "Nonsense" & z_score_top >= 2.0))



# Extended data Figure 1g: 

apc_cbe_after_1570 <- apc_cbe |> 
  filter(gene_names == "APC" & residue_number >= 1570 & categorized_mutation == "Nonsense")

apc_cbe_after_1570_plot <- ggplot(
  apc_cbe_after_1570,
  aes(x = z_score_top, y = categorized_mutation, fill = categorized_mutation)
) +geom_density_ridges(scale = 0.95) +
  theme_ridges() +
  scale_x_continuous(breaks = seq(-6, 12, by = 2), limits = c(-6, 12)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.2))) + theme_pubr() +
  theme(legend.position = "none")

apc_cbe_after_1570_plot
ggsave(
  filename = "Categories of mutations/apc_cbe_after_1570_false_positives.pdf",
  plot = apc_cbe_after_1570_plot, 
  height = 10, 
  width = 10, 
  units = "in",
  dpi = 600
)

nrow(apc_cbe_after_1570 |> filter(categorized_mutation == "Nonsense" & residue_number >= 1570))
nrow(apc_cbe_after_1570 |> filter(categorized_mutation == "Nonsense" & residue_number >= 1570 & z_score_top >= 2.0))
