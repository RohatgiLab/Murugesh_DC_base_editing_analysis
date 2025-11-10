#### β-catenin Destruction Complex (BDC) Base editing Screens
Base editing screens of Destruction Complex Components
The WNT/β-catenin pathway regulates tissue development and homeostasis, and its dysregulation drives many human diseases, including colorectal cancer. Central to this pathway is the β-catenin destruction complex (BDC), which controls β-catenin stability through a network of protein–protein interactions.
To better understand how specific sequence elements within BDC components regulate signaling in vivo, we conducted CRISPR base editing screens targeting multiple BDC genes. These screens enabled systematic mutational analysis at endogenous loci, revealing activating and regulatory mutations that refine our understanding of WNT/β-catenin pathway control.

The R script file contains the computational workflow used to process and analyze base editing screen data. The analysis includes:
 1. Processing of raw next-generation sequencing (NGS) data
 2. Quantification of raw read counts of guideRNAs
 3. Statistical modeling to identify variants that alter WNT signaling activity
 4. Annotation of predicted mutations based on genomic position and codon context
 5. Categorization of mutations by predicted functional class (e.g., synonymous, missense, nonsense, or splice site)

Together, these analyses provide a comprehensive view of how specific mutations affect β-catenin destruction complex function and WNT/β-catenin pathway regulation.

#### System requirements
1. Operating system tested:
   macOS Sequoia 15.6.1
2. R version 4.3.2.
3. Packages (CRAN / Bioconductor)
   dplyr 1.1.4, tibble 3.2.1, janitor 2.2.1, readxl 1.4.5
   openxlsx 4.2.8, ggplot2 3.5.2, ggpubr, data.table 1.17.0
   Biostrings 2.70.3, GenomicAlignments 1.38.2 ,Rsamtools 2.18.0, 
   QuasR 1.42.1, ShortRead 1.60.0, Rsubread 2.16.1
   stringr 1.5.1, plotly 4.10.4
4. Non-standard hardware
   None required. The pipeline runs on standard personal computers or servers.

#### Installation guide: Instructions
A. Prerequisites
 R version 4.3.2 installed
 Internet access to install R packages
 RStudio for a user-friendly interface
B. Install CRAN packages
install.packages(c(
  "dplyr", "tibble", "janitor", "readxl", "openxlsx",
  "ggplot2", "ggpubr", "data.table", "stringr", "plotly"
))
C. Install Bioconductor packages:
Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "Biostrings", "GenomicAlignments", "GenomicFeatures",
  "Rsamtools", "QuasR", "ShortRead", "Rsubread"
))

#### Expected run time: 
Full pipeline (all steps, full dataset): ~2–3 hours on a standard desktop computer (16 GB RAM, 4–8 CPU cores).


