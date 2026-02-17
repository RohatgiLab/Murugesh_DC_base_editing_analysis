# Hap1 validation data plotting with subsampled values: Figure 2B, Figure 5A:
library(viridis)
library(viridisLite)
library(ggbreak)
library(plotly)
library(readxl)
library(openxlsx)
library(janitor)
library(dplyr)
library(tidyr)

# Selected CTNNB1 guides were validated in two sets: set1 and set2: 

ctnnb1_set1 <- read_excel("ctnnb1_set1_hap1_appended.xlsx") |> data.table()

# hard coded 20 comes from sampling 20 mean values from the fluorescence data from hap1 validation:
control_mean <- ctnnb1_set1$means_PE_A[ctnnb1_set1$sample == "CTNNB1_Set1_NTC1_ABE_WNT_treated"]
control_mean
control_sd <- ctnnb1_set1$stdev[ctnnb1_set1$sample == "CTNNB1_Set1_NTC1_ABE_WNT_treated" ]
control_sd
ctnnb1_set1 <- ctnnb1_set1 |> 
  mutate(
    ratio = means_PE_A/control_mean
  )
View(ctnnb1_set1)


ctnnb1_set1 <- ctnnb1_set1 |> 
  mutate(
    percent_fluorescence = (means_PE_A / control_mean - 1)* 100, 
    percent_stdev = ratio * sqrt(
      (stdev/means_PE_A)^2 +
        (control_sd/control_mean)^2
    )
  ) |> mutate(
    percent_stdev = (percent_stdev *100), .after = ratio
  )

View(ctnnb1_set1)



ctnnb1_set2 <- read_excel("ctnnb1_se_wnt_set2_hap1_appended.xlsx") |> data.table()
ctnnb1_set2 <- ctnnb1_set2 |> mutate(
  stdev = se_PE_A * sqrt(20)
)

set2_control_mean <- ctnnb1_set2$means_PE_A[ctnnb1_set2$sample == "CTNNB1_Set2_NTC2_CBE+"]
set2_control_sd <- ctnnb1_set2$stdev[ctnnb1_set2$sample == "CTNNB1_Set2_NTC2_CBE+" ]


ctnnb1_set2 <- ctnnb1_set2 |> 
  mutate(
    ratio = (means_PE_A / set2_control_mean), 
    percent_fluorescence = ((means_PE_A/set2_control_mean) -1)*100,
    percent_stdev = ratio * sqrt(
      (stdev/means_PE_A)^2 +
        (set2_control_sd/set2_control_mean)^2)
  ) |> rename(ids_hap1 = ids) |> mutate(
    percent_stdev = (percent_stdev *100), .after = ratio
  )

View(ctnnb1_set2)

ctnnb1_full_set <- rbind(ctnnb1_set1, ctnnb1_set2)
View(ctnnb1_full_set)

write.xlsx(ctnnb1_full_set, "Subsampling_FACS_files/ctnnb1_full_set.xlsx")

ctnnb1_for_plotting <- ctnnb1_full_set |> 
  filter(!str_detect(sample, "ABE|CBE"))


ctnnb1_for_plotting <- ctnnb1_for_plotting |>
  mutate(residue_number = str_extract(amino_acid_edits, "\\d+")) |>
  mutate(residue_number = as.integer(residue_number))

View(ctnnb1_for_plotting)

ctnnb1_for_plotting <- ctnnb1_for_plotting |>
  filter(!str_detect(amino_acid_edits, "[Tt][Ee][Rr]"))


# merging dataset with validation data: 

hap1_data_Ganesh <- read_excel("Supplementary_Table 2_Validation_Guides.xlsx", 
                               sheet = "CTNNB1") |> select(1, 3, 6, 7) |> 
  rename(Screen = "Base_Editing_Screen", ids_hap1 = "guideRNA_identification_number", e_hap1_7s_validation_phenotype = "eHAP1-7S validation phenotype"  )



ctnnb1_for_plotting <- left_join(ctnnb1_for_plotting, hap1_data_Ganesh, by = c("Screen", "ids_hap1"))

ctnnb1_for_plotting <- ctnnb1_for_plotting |> filter(!is.na(enriched_gate))

ctnnb1_for_plotting <- ctnnb1_for_plotting |> mutate(categories = case_when(
  str_detect(e_hap1_7s_validation_phenotype, "reduction") & str_detect(enriched_gate, "Top") ~ "false_positive", 
  str_detect(e_hap1_7s_validation_phenotype, "elevation") & str_detect(enriched_gate, "Bot") ~ "false_positive",
  .default = "validated"
  
))

ctnnb1_for_plotting <- ctnnb1_for_plotting |> 
  filter(!categories == "false_positive")

View(ctnnb1_for_plotting)
domain_boundaries <- c(148,180, 223, 264, 306, 349, 390, 429, 473, 519, 582, 623, 664, 691, 781)

ctnnb1_validation_plot <- ggplot(ctnnb1_for_plotting, aes(x = residue_number, y = percent_fluorescence, fill = enriched_gate)) +
  geom_col(width = 5, color = "black", size = 0.1) + geom_errorbar(
    aes(
      ymin = percent_fluorescence - percent_stdev,
      ymax = percent_fluorescence + percent_stdev
    ),
    width = 1,
    linewidth = 0.2
  ) +
  theme_pubr() + scale_x_continuous(expand = expansion(0), breaks = seq(0, 800, 100), limits = c(0, 800)) +
  scale_fill_brewer(palette = "Set1", direction = -1) + theme(legend.position = "none") +
  geom_hline(yintercept = 0, size = 0.1) + scale_y_continuous(breaks = seq(-100, 800, 50), limits =c(-100, 800)) +
  scale_y_break(c(30, 500),scales = 0.1) +
  geom_vline(xintercept = domain_boundaries, linetype = "dashed", color = "green")

ctnnb1_validation_plot

ggsave(
  filename = "Subsampling_FACS_files/CTNNB1_fluorescence_data_with_error_bars.pdf",
  plot = ctnnb1_validation_plot,
  units = "in",
  height = 4,
  width = 10,
  dpi = 600
)

write.xlsx(ctnnb1_for_plotting, "Subsampling_FACS_files/ctnnb1_for_plotting.xlsx")


# Axin1 guides were also validated in two sets: set1 and set2:  

axin1_set1 <- read_excel("axin1_set01_hap1_appended.xlsx") |> data.table()
axin1_set1 <- axin1_set1 |> 
  mutate(
    stdev = se_PE_A * sqrt(20)
  )

# hard coded 20 comes from sampling 20 mean values from the fluorescence data from hap1 validation:
set1_axin1_mean <- axin1_set1$means_PE_A[axin1_set1$sample == "AXIN1_Set1_78+_1stSet"]
set1_axin1_mean

set1_axin1_sd <- axin1_set1$stdev[axin1_set1$sample == "AXIN1_Set1_78+_1stSet"]
set1_axin1_sd

axin1_set1 <- axin1_set1 |> 
  mutate(
    ratio = means_PE_A / set1_axin1_mean
  )

View(axin1_set1)

axin1_set1 <- axin1_set1 |> 
  mutate(
    percent_fluorescence = (means_PE_A / set1_axin1_mean - 1) * 100, 
    percent_stdev = ratio * sqrt(
      (stdev/means_PE_A)^2 +
        (set1_axin1_sd/set1_axin1_mean)^2
    )
  ) |> mutate(
    percent_stdev = (percent_stdev * 100), .after = ratio
  )
axin1_set1 <- axin1_set1 |> select(-ids)

View(axin1_set1)

axin1_set2 <- read_excel("axin1_set02_hap1_appended.xlsx") |> data.table()

axin1_set2 <- axin1_set2 |> mutate(
  stdev = se_PE_A * sqrt(20)
)

View(axin1_set2)

set2_axin1_mean <- axin1_set2$means_PE_A[axin1_set2$sample == "AXIN1_Set2_76+_2ndSet"]
set2_axin1_mean

set2_axin1_sd <- axin1_set2$stdev[axin1_set2$sample == "AXIN1_Set2_76+_2ndSet"]
set2_axin1_sd

axin1_set2 <- axin1_set2 |>
  mutate(
    ratio = means_PE_A / set2_axin1_mean,
    percent_fluorescence = (ratio - 1) * 100,
    percent_stdev = ratio * sqrt(
      (stdev / means_PE_A)^2 +
        (set2_axin1_sd / set2_axin1_mean)^2
    ),
    percent_stdev = percent_stdev * 100,
    .after = ratio
  )

axin1_set2 <- axin1_set2 |> select(-ids) 
names(axin1_set2)
axin1_set2 <- axin1_set2 |> dplyr::relocate(label_on_tubes, .after = se_PE_A)

names(axin1_set1) == names(axin1_set2)

axin1_full_set <- rbind(axin1_set1, axin1_set2)
View(axin1_full_set)
write.xlsx("Subsampling_FACS_files/axin1_full_set.xlsx")

axin1_for_plotting <- axin1_full_set |> 
  filter(!str_detect(Screen, "NTC"))

View(axin1_for_plotting)

axin1_for_plotting <- axin1_for_plotting |> 
  mutate(residue_number = as.numeric(str_extract(amino_acid_edits, "\\d+")))


axin1_for_plotting <- axin1_for_plotting |> filter(!str_detect(amino_acid_edits,"[Tt][Ee][Rr]" ))




hap1_data_Ganesh_axin1 <- read_excel("Supplementary_Table 2_Validation_Guides.xlsx", 
                               sheet = "AXIN1") |> select(1, 3, 6, 7) |> 
  rename(Screen = "Base_Editing_Screen", ids_hap1 = "guideRNA_identification_number", ehap1_7S_validation_phenotype = "eHAP1_7S_validation_phenotype")


axin1_for_plotting <- axin1_for_plotting |> rename(ids_hap1 = "ids.y")
axin1_for_plotting <- left_join(axin1_for_plotting, hap1_data_Ganesh, by = c("Screen", "ids_hap1"))
write.xlsx(axin1_for_plotting, "Subsampling_FACS_files/axin1_for_plotting.xlsx")
axin1_for_plotting <- axin1_for_plotting |> filter(!is.na(enriched_gate))

axin1_for_plotting <- read_excel("Subsampling_FACS_files/axin1_for_plotting.xlsx")
View(axin1_for_plotting)

axin1_for_plotting <- axin1_for_plotting |> filter(gene != "CTNNB1")
View(axin1_for_plotting)

axin1_for_plotting <- axin1_for_plotting |> mutate(categories = case_when(
  str_detect(ehap1_7S_validation_phenotype, "reduction") & str_detect(enriched_gate, "Top") ~ "false_positive", 
  str_detect(ehap1_7S_validation_phenotype, "elevation") & str_detect(enriched_gate, "Bot") ~ "false_positive",
  .default = "validated"
  
))

View(axin1_for_plotting)

axin1_for_plotting <- axin1_for_plotting |> 
  filter(!categories == "false_positive")


axin1_domain_boundaries <- c(1,81, 213, 383, 402, 462, 480, 780, 862)
axin1_validation_plot <- ggplot(axin1_for_plotting, aes(x = residue_number, y = percent_fluorescence, fill = enriched_gate)) +
  geom_col(width = 5, color = "black", size = 0.1) + geom_errorbar(
    aes(
      ymin = percent_fluorescence - percent_stdev,
      ymax = percent_fluorescence + percent_stdev
    ),
    width = 1,
    linewidth = 0.2
  ) + theme_pubr() + scale_x_continuous(expand = expansion(0), breaks = seq(0, 862, 100), limits = c(0, 862)) +
  scale_y_continuous(limits = c(-100, 50)) +
  scale_fill_brewer(palette = "Set1", direction = -1) + theme(legend.position = "none") +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_vline(xintercept = axin1_domain_boundaries, linetype = "dashed", color = "green")

axin1_validation_plot
ggsave(
  filename = "Subsampling_FACS_files/Axin1_validation_plot.pdf",
  plot = axin1_validation_plot,
  units = "in",
  height = 4,
  width = 10,
  dpi = 600
)

