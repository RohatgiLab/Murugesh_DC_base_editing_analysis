# Z_score_correlations: Extended Data Figure 1c:
# Please use raw read counts file and the supplementary table 1: 


library(dplyr)
library(janitor)
library(readxl)
library(openxlsx)
library(ggplot2)
library(data.table)
library(ggpubr)
library(plotly)
library(ggstatsplot)

# DC_cbe library raw read counts

dc_cbe_wnt <- read_excel("dc_cbe_wnt_raw_countsavg_LFC.xlsx") |> data.table()

master_file_cbe <- read_excel("CBE_ABE_libraries_analysis_revist/Supplementary Table 1_Master File_All Guides.xlsx", 
                          sheet = "WNT_DC_CBE Screen") |> data.table()

master_file_cbe <- master_file_cbe[, c(1:3, 14:15)]


# log2 transformation of normalized reads:

log2_transform <- function(df){
  df_log <- df |> mutate(across(ends_with("rpm"), function(x) log(x+2, 2), .names = "{.col}_log"))
  return(df_log)
}

dc_cbe_wnt_log <- log2_transform(dc_cbe_wnt)

View(dc_cbe_wnt_log)

# Calculating log fold change accounting for replicates:

dc_cbe_wnt_lfc <- dc_cbe_wnt_log |> 
  mutate(B_1_LFC = (B_1_rpm_log - U_1_rpm_log), 
         T_1_LFC = (T_1_rpm_log - U_1_rpm_log), 
         B_2_LFC = (B_2_rpm_log - U_2_rpm_log), 
         T_2_LFC = (T_2_rpm_log - U_2_rpm_log))

View(dc_cbe_wnt_lfc)

# splitting the cbe dataset with replicates 01 and replicates 02: 

dc_cbe_wnt_rep_01 <- dc_cbe_wnt_lfc[, c(1, 2, 4,6, 8, 10, 12, 14, 16, 18, 20, 21)]

dc_cbe_wnt_rep_02 <- dc_cbe_wnt_lfc[, c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 22, 23)]


# appending gene column into both data frames: dc_cbe_wnt_rep_01 and dc_cbe_wnt_rep_02:  

cbe_genes <- c("CTNNB1", "APC", "Axin1", "Gsk3b", "control")

dc_cbe_wnt_rep_01 <- dc_cbe_wnt_rep_01 |> 
  mutate(gene_names = rep(cbe_genes, c(868,2449,1350,446,282)), .after = ids)

dc_cbe_wnt_rep_02 <- dc_cbe_wnt_rep_02 |> 
  mutate(gene_names = rep(cbe_genes, c(868, 2449, 1350, 446, 282)), .after = ids)

plot_01 <- ggplot(dc_cbe_wnt_rep_01, aes(y= U_1_rpm_log)) + geom_boxplot()

ggplotly(plot_01)

plot_02 <- ggplot(dc_cbe_wnt_rep_02, aes(y = U_2_rpm_log)) + geom_boxplot()

ggplotly(plot_02)


# Filtering outliers based on IQR values:

dc_cbe_wnt_rep_01 <- dc_cbe_wnt_rep_01 |> filter(U_1_rpm_log >= 5.56 & U_1_rpm_log <= 9.22)
nrow(dc_cbe_wnt_rep_01)

dc_cbe_wnt_rep_02 <- dc_cbe_wnt_rep_02 |> filter(U_2_rpm_log >= 5.53 & U_2_rpm_log <= 9.28)
nrow(dc_cbe_wnt_rep_02)

# Okay, now I can get control statistics for both dc_cbe_wnt_rep_01 and dc_cbe_wnt_rep_02 as well:-

control_stats <- function(data, gene_colname, columns_to_mean_sd){
  
  control_data <- data |> filter(.data[[gene_colname]] == "control")
  
  control_stats <- control_data |> summarize(across(
    all_of(columns_to_mean_sd), 
    list(mean = function(x) mean(x, na.rm = TRUE),
         sd = function(x) sd(x, na.rm = TRUE))
    
  ))
  return(control_stats)
}

cbe_rep_01_control <- control_stats(
  data = dc_cbe_wnt_rep_01, 
  gene_colname = "gene_names",
  columns_to_mean_sd = c("B_1_LFC", "T_1_LFC")
)

# Let me Calculate z-scores for both B_1_LFC and T_1_LFC columns now:

dc_cbe_wnt_rep_01 <- dc_cbe_wnt_rep_01 |> mutate(
  B_1_zscores = (B_1_LFC - cbe_rep_01_control$B_1_LFC_mean) /cbe_rep_01_control$B_1_LFC_sd, 
  T_1_zscores = (T_1_LFC - cbe_rep_01_control$T_1_LFC_mean)/cbe_rep_01_control$T_1_LFC_sd
)


# Calculating z-scores in case of dc_cbe_wnt_rep_02: 

cbe_rep_02_control <- control_stats(
  data = dc_cbe_wnt_rep_02, 
  gene_colname = "gene_names",
  columns_to_mean_sd = c("B_2_LFC", "T_2_LFC")
)

cbe_rep_02_control

dc_cbe_wnt_rep_02 <- dc_cbe_wnt_rep_02 |> mutate(
  B_2_zscores = (B_2_LFC - cbe_rep_02_control$B_2_LFC_mean) / cbe_rep_02_control$B_2_LFC_sd, 
  T_2_zscores = (T_2_LFC - cbe_rep_02_control$T_2_LFC_mean) / cbe_rep_02_control$T_2_LFC_sd
)

dc_cbe_rep_01_02 <- inner_join(dc_cbe_wnt_rep_01, dc_cbe_wnt_rep_02, by = "ids")

# Plotting z_score correlation plots for CBE library:

cbe_bottom_zscores_plot <- ggplot(dc_cbe_rep_01_02, aes(x = B_1_zscores, y = B_2_zscores)) + xlim (-10, 15) + ylim (-10, 10) +geom_point(color = "palegreen3", size = 1.0) + theme_pubr()
cbe_bottom_zscores_plot
cbe_top_zscores_plot <- ggplot(dc_cbe_rep_01_02, aes(x = T_1_zscores, y = T_2_zscores)) +xlim(-10, 10) + ylim(-10, 10) + geom_point(color = "cornflowerblue", size = 1.0) + theme_pubr() 
cbe_top_zscores_plot

write.xlsx(dc_cbe_rep_01_02, "CBE_ABE_libraries_analysis_revist/dc_cbe_rep_01_02_zscores.xlsx")


# Analyzing DC_ABE_library_zscores of each replicates: 

dc_abe_data <- read_excel("dc_abe_wnt_raw_counts_avg_LFC.xlsx") |> data.table()

dc_abe_rep_01 <- dc_abe_data |> select(ids, B_1, T_1, U_1)
dc_abe_rep_02 <- dc_abe_data |> select(ids, B_2, T_2, U_2)

View(dc_abe_rep_01)
View(dc_abe_rep_02)

# Normalizing to reads per million reads and log2 transformation:

rpm_normalization <- function(df){
  df_rpm <- df |> mutate(across(-1, function(x) x/sum(x) * 10^6, .names = "{.col}_rpm"))
  return(df_rpm)
}

dc_abe_rep_01 <- rpm_normalization(dc_abe_rep_01)
View(dc_abe_rep_01)

dc_abe_rep_02 <- rpm_normalization(dc_abe_rep_02)
View(dc_abe_rep_02)

dc_abe_rep_01_log <- log2_transform(dc_abe_rep_01)
dc_abe_rep_02_log <- log2_transform(dc_abe_rep_02)


# Appending gene_names : 
abe_genes <- c("CTNNB1", "APC", "Axin1", "Gsk3b", "control")

dc_abe_rep_01_log <- dc_abe_rep_01_log |> 
  mutate(gene_names = rep(abe_genes, c(795,2556,896,416,282)), .after = "ids")

dc_abe_rep_02_log <- dc_abe_rep_02_log |> 
  mutate(gene_names = rep(abe_genes, c(795,2556,896,416,282)), .after = "ids")

abe_plot_01 <- ggplot(dc_abe_rep_01_log, aes(y = U_1_rpm_log)) + geom_boxplot()
ggplotly(abe_plot_01)

abe_plot_02 <- ggplot(dc_abe_rep_02_log, aes(y = U_2_rpm_log)) + geom_boxplot()

ggplotly(abe_plot_02)

# Removing outliers: 

dc_abe_rep_01_log <- dc_abe_rep_01_log |> 
  filter(U_1_rpm_log >= 3.87 & U_1_rpm_log <= 10.62)

dc_abe_rep_02_log <- dc_abe_rep_02_log |> 
  filter(U_2_rpm_log >= 3.89 & U_2_rpm_log <= 10.62)

# Calculating log fold change:

dc_abe_rep_01_lfc <- dc_abe_rep_01_log |> 
  mutate(
    B_1_LFC = (B_1_rpm_log - U_1_rpm_log), 
    T_1_LFC = (T_1_rpm_log - U_1_rpm_log)
  )

dc_abe_rep_02_lfc <- dc_abe_rep_02_log |> 
  mutate(
    B_2_LFC = (B_2_rpm_log - U_2_rpm_log),
    T_2_LFC = (T_2_rpm_log - U_2_rpm_log)
  )



# Control statistics: 

abe_rep_01_control <- control_stats(
  data = dc_abe_rep_01_lfc, 
  gene_colname = "gene_names",
  columns_to_mean_sd = c("B_1_LFC", "T_1_LFC")
)

dc_abe_rep_01_zscores <- dc_abe_rep_01_lfc |> 
  mutate(
    B_1_zscores = (B_1_LFC - abe_rep_01_control$B_1_LFC_mean) / abe_rep_01_control$B_1_LFC_sd, 
    T_1_zscores = (T_1_LFC - abe_rep_01_control$T_1_LFC_mean) / abe_rep_01_control$T_1_LFC_sd
  )

abe_rep_02_control <- control_stats(
  data = dc_abe_rep_02_lfc, 
  gene_colname = "gene_names", 
  columns_to_mean_sd = c("B_2_LFC", "T_2_LFC")
)


dc_abe_rep_02_zscores <- dc_abe_rep_02_lfc |> 
  mutate(
    B_2_zscores = (B_2_LFC - abe_rep_02_control$B_2_LFC_mean) / abe_rep_02_control$B_2_LFC_sd, 
    T_2_zscores = (T_2_LFC - abe_rep_02_control$T_2_LFC_mean) / abe_rep_02_control$T_2_LFC_sd
  )

View(dc_abe_rep_01_zscores)
View(dc_abe_rep_02_zscores)

dc_abe_rep_01_02_zscores <- inner_join(dc_abe_rep_01_zscores, dc_abe_rep_02_zscores, by = "ids")

# Plotting z_score correlation plots: 

abe_bottom_plot <- ggplot(dc_abe_rep_01_02_zscores, aes(x = B_1_zscores, y = B_2_zscores)) + geom_jitter(width = 2, height = 1.5, color = "palegreen3", size = 1) + xlim(-10, 10) + ylim(-10, 10) + theme_pubr()
abe_bottom_plot
abe_top_plot <- ggplot(dc_abe_rep_01_02_zscores, aes(x = T_1_zscores, y = T_2_zscores)) + geom_point(color = "cornflowerblue", size = 1) + xlim(-10, 10) + ylim(-10, 10) + theme_pubr()
abe_top_plot

ggsave(
  filename = "CBE_ABE_libraries_analysis_revist/Extended_data_figure_01/cbe_bottom_replicates_zscores.pdf", 
  plot = cbe_bottom_zscores_plot, 
  dpi = 600,
  height = 4, 
  width = 5, 
  units = "in"
  
)

ggsave(
  filename = "CBE_ABE_libraries_analysis_revist/Extended_data_figure_01/cbe_top_replicates_zscores.pdf", 
  plot = cbe_top_zscores_plot, 
  dpi = 600,
  height = 4, 
  width = 5, 
  units = "in"
  
)

ggsave(
  filename = "CBE_ABE_libraries_analysis_revist/Extended_data_figure_01/abe_bottom_replicates_zscores_03.pdf", 
  plot = abe_bottom_plot, 
  dpi = 600,
  height = 4, 
  width = 5, 
  units = "in"
  
)

ggsave(
  filename = "CBE_ABE_libraries_analysis_revist/Extended_data_figure_01/abe_top_replicates_zscores_03.pdf", 
  plot = abe_top_plot, 
  dpi = 600,
  height = 4, 
  width = 5, 
  units = "in"
  
)



# correlation coefficients (Spearman rank) were computed between z-score profiles of replicate 1 and replicate 2.

cbe_correlations <- ggcorrmat(
  data = dc_cbe_rep_01_02, 
  cor.vars = c("B_1_zscores", "B_2_zscores", "T_1_zscores", "T_2_zscores")
)
cbe_correlations


abe_correlations <-  ggcorrmat(
  data = dc_abe_rep_01_02_zscores, 
  cor.vars = c("B_1_zscores", "B_2_zscores", "T_1_zscores", "T_2_zscores")
)

abe_correlations

ggsave(
  filename = "CBE_ABE_libraries_analysis_revist/cbe_correlations.pdf", 
  plot = cbe_correlations, 
  height = 3,
  width = 3, 
  units = "in"
)

ggsave(
  filename = "CBE_ABE_libraries_analysis_revist/abe_correlations.pdf", 
  plot = abe_correlations, 
  height = 3,
  width = 3, 
  units = "in"
)



