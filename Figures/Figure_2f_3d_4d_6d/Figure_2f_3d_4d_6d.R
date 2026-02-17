# Figures 2f, 3d, 4d, 6d
# datasets attached: 

library(readxl)
library(openxlsx)
library(janitor)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gt)
library(cowplot)
library(emmeans)
library(kableExtra)
library(viridis)
library(tibble)


# Figure 2F: 

fig2f <- read_excel("BLI_data_consolidated/Figure_2F_data.xlsx") |> data.table()
  

## tidy data format for plotting:

fig2f_longer <- fig2f |> 
  pivot_longer(
    cols = -1, 
    names_to = "variants",
    values_to = "nm_shift"
  )

fig2f_longer$variants <- factor(fig2f_longer$variants, levels = c("WT", "E462K", "C466Y/A467T", "K435G", "R469H"))
names(fig2f_longer)
## Plotting: 

fig2f_plot <- ggplot(fig2f_longer, aes(x = Time, y = nm_shift, color = variants)) + 
  geom_line(linewidth = 0.9) + 
  scale_color_viridis(discrete = TRUE, option = "D") + 
  theme_pubr()  +  scale_y_continuous(limits =c(-0.02, 0.4)) +geom_vline(xintercept = 240, linetype = "dotted", color = "magenta", linewidth = 0.3) +
  scale_x_continuous(limits = c(220,360), breaks = seq(220,360, 40))

fig2f_plot # colors and x-axis label (transformed) in illustrator: 



# Figure 3d_ ICAT plot:

figure_3d_icat <- read_excel("BLI_data_consolidated/Figure_3d_ICAT_data.xlsx") |> data.table()
names(figure_3d)
# tidy format:

fig3d_longer_icat <- figure_3d_icat |> 
  pivot_longer(
    cols = -1, 
    names_to = "variants",
    values_to = "nm_shift"
  ) 

fig3d_longer_icat$variants <- factor(fig3d_longer_icat$variants, levels = c("WT", "E571K", "K435G/M436V"))

fig3d_icat_plot <- ggplot(fig3d_longer_icat, aes(x = time, y = nm_shift, color = variants)) + 
  geom_line(linewidth = 0.9) + 
  scale_color_viridis(discrete = TRUE, option = "D") + 
  theme_pubr()  +  scale_y_continuous(limits =c(-0.02, 0.5)) +geom_vline(xintercept = 240, linetype = "dotted", color = "magenta", linewidth = 0.3) +
  scale_x_continuous(limits = c(220,360), breaks = seq(220,360, 40))

fig3d_icat_plot # colors and x-axis label (transformed) in illustrator: 

# Figure 3d_TCF4 plot: 

figure_3d_tcf4 <- read_excel("BLI_data_consolidated/Figure_3d_TCF4_data.xlsx") |> data.table()

names(figure_3d_tcf4)
# tidy format:

fig3d_longer_tcf4 <- figure_3d_tcf4 |> 
  pivot_longer(
    cols = -1, 
    names_to = "variants",
    values_to = "nm_shift"
  ) 

fig3d_longer_tcf4$variants <- factor(fig3d_longer_tcf4$variants, levels = c("WT", "E571K", "K435G/M436V"))

fig3d_tcf4_plot <- ggplot(fig3d_longer_tcf4, aes(x = time, y = nm_shift, color = variants)) + 
  geom_line(linewidth = 0.9) + 
  scale_color_viridis(discrete = TRUE, option = "D") + 
  theme_pubr()  +  scale_y_continuous(limits =c(-0.02, 0.4)) +geom_vline(xintercept = 245, linetype = "dotted", color = "magenta", linewidth = 0.3) +
  scale_x_continuous(limits = c(220,360), breaks = seq(220,360, 40))

fig3d_tcf4_plot # colors and x-axis label (transformed) in illustrator: 


# Figure 4d_ICAT:

fig4d_icat <- read_excel("BLI_data_consolidated/Figure_4d_ICAT_data.xlsx") |> data.table()
## tidy format
fig4d_longer_icat <- fig4d_icat |> 
  pivot_longer(
    cols = -1, 
    names_to = "variants",
    values_to = "nm_shift"
  ) 

fig4d_longer_icat$variants <- factor(fig4d_longer_icat$variants, levels = c("WT", "C619Y/E620K", "D624N/E626K/E629K"))

fig4d_icat_plot <- ggplot(fig4d_longer_icat, aes(x = time, y = nm_shift, color = variants)) + 
  geom_line(linewidth = 0.9) + 
  scale_color_viridis(discrete = TRUE, option = "D") + 
  theme_pubr()  +  scale_y_continuous(limits =c(-0.02, 0.5)) +geom_vline(xintercept = 240, linetype = "dotted", color = "magenta", linewidth = 0.3) +
  scale_x_continuous(limits = c(220,360), breaks = seq(220,360, 40))

fig4d_icat_plot # colors and x-axis label (transformed) in illustrator: 


# Fig4d_TCF4 plot: 

fig4d_tcf4 <- read_excel("BLI_data_consolidated/Figure_4d_TCF_data.xlsx") |> data.table()
names(fig4d_tcf4)
## tidy format
fig4d_longer_tcf4 <- fig4d_tcf4 |> 
  pivot_longer(
    cols = -1, 
    names_to = "variants",
    values_to = "nm_shift"
  ) 

fig4d_longer_tcf4$variants <- factor(fig4d_longer_tcf4$variants, levels = c("WT", "C619Y/E620K", "D624N/E626K/E629K"))

fig4d_tcf4_plot <- ggplot(fig4d_longer_tcf4, aes(x = time, y = nm_shift, color = variants)) + 
  geom_line(linewidth = 0.9) + 
  scale_color_viridis(discrete = TRUE, option = "D") + 
  theme_pubr()  +  scale_y_continuous(limits =c(-0.02, 0.4)) +geom_vline(xintercept = 245, linetype = "dotted", color = "magenta", linewidth = 0.3) +
  scale_x_continuous(limits = c(220,360), breaks = seq(220,360, 40))

fig4d_tcf4_plot # colors and x-axis label (transformed) in illustrator: 



# Figure 6d:

fig6d  <- read_excel("BLI_data_consolidated/Fig6D_data_consolidated.xls") |> data.table()

## tidy data: 

fig6d_longer <- fig6d |> 
  pivot_longer(
    cols = -1, 
    names_to = "variants",
    values_to = "nm_shift",
  )

fig6d_longer$variants <- factor(fig6d_longer$variants, levels = c(
  "WT","L479W","V475I","V475A", "V475I/V479W" , "L471A/D472A/H474A"  
))

fig6d_plot <- ggplot(fig6d_longer, aes(x = Time, y = nm_shift, color = variants)) + 
  geom_line(linewidth = 0.9) + 
  scale_color_viridis(discrete = TRUE, option = "D") + 
  theme_pubr()  +  scale_y_continuous(limits =c(-0.02, 0.5)) +geom_vline(xintercept = 240, linetype = "dotted", color = "magenta", linewidth = 0.3) +
  scale_x_continuous(limits = c(220,360), breaks = seq(220,360, 40))

fig6d_plot # colors and x-axis label (transformed) in illustrator: 
