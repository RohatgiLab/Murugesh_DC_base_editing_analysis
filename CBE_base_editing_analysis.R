# Base editing screen analysis of Destruction complex CBE libraries:

## libraries loaded:
library(dplyr)
library(tibble)
library(janitor)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(data.table)
library(Biostrings)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(QuasR)
library(ShortRead)
library(Rsubread)
library(stringr)
library(plotly)


## Processing reads from CBE-library; getting rid of backbone sequence, to keep only the guideRNA sequence:
### Function to process reads: process_fastq_files function trims reads based on defining flanking patterns
### (Lpattern and Rpattern) corresponding to the vector backbone sequences. Each input Fastq file is scanned
### for these sequence motifs and only the intervening region (the gRNA) is retained. Processed reads are 
### written to a new subdirectory (processed_reads/). These processed sequences can then be used for alignment
### quantification, and downstream base editing enrichment analysis. 

process_fastq_files <- function(input_files,
                                input_dir, 
                                output_suffix ="processed.fastq.gz",
                                nBases = 2,
                                Lpattern = "GCTTTATATATCTTGTGGAAAGGACGAAA" 
                                ){
  output_dir <- file.path(input_dir, "processed_reads")
  dir.create(output_dir, showWarnings = FALSE)
  
  for(file in input_files){
    input_path <- file.path(input_dir, file)
    output_path <- file.path(output_dir, sub("fastq.gz", output_suffix, file))
    preprocessReads(input_path, output_path, nBases = nBases,Lpattern = Lpattern)
  }
  message("Processing completed, processed files area save in:", output_dir)
}

#### Input : gzipped FASTQ files stored in CBE_fastq_files/
#### Output: processed FASTQ files (gRNA-only) stored in CBE_fastq_files/processed_reads/ folder. 

cbe_fastq_files <- list.files("CBE_fastq_files/", pattern= "fastq.gz", all.files = TRUE) # listed files

process_fastq_files(
  input_files = cbe_fastq_files, 
  input_dir = "CBE_fastq_files/", 
  output_suffix = "processed.fastq.gz", 
  nBases = 2, 
  Lpattern = "GCTTTATATATCTTGTGGAAAGGACGAAA"
)

## Checking the output of processed reads
### This step verifies that the FASTQ pre-processing step completed successfully and the
### the expected guide RNA (gRNA) sequences were retained. The script lists all processed FASTQ files 
### in the processed_reads/ sub directory that matches suffix processed.fastq.gz. The resulting
### sequences are stored as elements in a list (reads_list), with each list element corresponding to
### one processed FASTQ file. 

cbe_processed_files <- list.files("CBE_fastq_files/processed_reads/", pattern = "processed.fastq.gz", full.names = TRUE)

reads_list <- list()
for(file in cbe_processed_files){
  fastq_file <- readFastq(file)
  reads_list[[file]] <- sread(fastq_file)
}
reads_list

## Aligning reads to a reference file:

### This function aligns processed guide (gRNA) sequencing reads to a specified reference file
### to generate alignment files suitable for downstream quantification.The process involves building
### a reference index, aligning each FASTQ file, producing sorted, indexed BAM files. 

align_fastq_files <- function(reference_file, fastq_directory, index_name = "abe_library", memory = 8000){
  buildindex(index_name, reference_file, memory = memory)
  fastq_files <- list.files(fastq_directory, pattern = "*.processed.fastq.gz", full.names = TRUE)
  for(fastq_file in fastq_files){
    base_name <- tools::file_path_sans_ext(basename(fastq_file))
    output_folder <- file.path(fastq_directory, base_name)
    dir.create(output_folder, showWarnings = FALSE)
    output_file <- file.path(output_folder, paste0(base_name, "aligned.bam"))
    align(index_name, fastq_file, output_format = "BAM", output_file = output_file)
    sorted_bam_file <- file.path(output_folder, paste0("sorted_", base_name))
    sortBam(output_file, sorted_bam_file)
    indexBam(paste0(sorted_bam_file, ".bam"))
    
  }
}

align_fastq_files(
  reference_file = "CBE_fastq_files/Reference_FASTA_FILE/DC_cbe_sgRNAs.fa", 
  fastq_directory = "CBE_fastq_files/processed_reads/", 
  index_name = "dc_cbe_library"
)

## This function collects all sorted BAM files matching the pattern ^sorted.bam$ from
## multiple subdirectories into a single destination folder. 
### Sorted BAM files are required for downstream analysis, such as counting
### reads per guide RNA or calculating read frequencies. 

move_sorted_bam_files <- function(main_folder, destination_folder){
  subfolders <- list.dirs(path = main_folder, full.names = TRUE, recursive = FALSE)
  dir.create(destination_folder, showWarnings = FALSE)
  for(subfolder in subfolders){
    sorted_bam_files <- list.files(subfolder, pattern ="^sorted.*\\.bam$", full.names = TRUE)
    
    for(bam_file in sorted_bam_files){
      file.copy(bam_file, destination_folder, overwrite = FALSE)
    }
    
  }
}

move_sorted_bam_files(
  main_folder = "CBE_fastq_files/processed_reads/", 
  destination_folder = "CBE_fastq_files/cbe_sorted_bam_files"
)

## This function processes all BAM files in a specified directory to
## calculate read counts for each guide RNA. For each BAM file, the resulting
## read frequencies are saved as individual excel files in a subdirectory named
## "read_counts".

process_and_save_bams <- function(bam_files_directory){
  output_folder <- file.path(bam_files_directory, "read_counts")
  dir.create(output_folder, showWarnings = FALSE)
  bam_files <- list.files(bam_files_directory, pattern = "\\.bam$", full.names = TRUE)
  for(bam_file in bam_files){
    alignments <- readGAlignments(bam_file)
    alignments_df <- as.data.frame(alignments)
    seqnames_table <- table(alignments_df$seqnames)
    seqnames_df <- as.data.frame(seqnames_table)
    output_excel <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(bam_file)), ".xlsx"))
    write.xlsx(seqnames_df, output_excel)
    message("processing and saving completed successfully", output_excel)
  }
  
}

process_and_save_bams(
  bam_files_directory = "CBE_fastq_files/cbe_sorted_bam_files/"
)

## This function reads multiple excel files, extracts a specific identifier
## based on the regular expression I used from each filename, and combines
## selected columns from all files into a single data table. It then saves the
## merged data as combined_data.xlsx. I moved to a folder and renamed as 
## CBE_raw_read_counts. 

extract_cbind_cols <- function(filenames){
  
  extract_identifier <- function(file_name){
    str_extract(file_name, "(?<=sorted_)DC_CBE_WT_[BTU]_[0-9]")
  }
  
  data_list <- list()
  
  for(i in 1:length(filenames)){
    file <- filenames[i]
    data <- read_excel(file) |> data.table()
    
    identifier <- extract_identifier(file)
    
    if(i == 1){
      data <- data[, c(1, 2)]
      colnames(data) <- c("ids", identifier)
    }else{
      data <- data[, 2, drop =FALSE]
      colnames(data) <- identifier
    }
    
    
    data_list[[i]] <- data
  }
  final_data <- do.call(cbind, data_list)
  save_path <- file.path(getwd() , "combined_data.xlsx")
  write.xlsx(final_data, save_path)
  print(paste("Data has been successfully combined and saved as 'combined_data.xlsx' in", getwd()))
}

cbe_excel_files <- list.files(
  "CBE_fastq_files/cbe_sorted_bam_files/read_counts/", pattern = ".xlsx",  full.names = TRUE)

extract_cbind_cols(filenames = cbe_excel_files)


## Load raw read count data generated from sequencing of CBE screens.
## Column names are manually simplified to clearly represent the different
## sorted populations: Bottom (B1, B2), Top (T1, T2), and Unsorted (U1, U2).

cbe_read_counts <- read_excel("CBE_fastq_files/Raw_read_counts/CBE_raw_read_counts.xlsx") |> 
  clean_names()
colnames(cbe_read_counts) <- c("ids", "cbe_B1", "CBE_B2", "CBE_T1", "CBE_T2", "CBE_U1", "CBE_U2")


## Calculate_rpm() function calculate reads-per-million (RPM) normalization of sequencing read
## counts for each condition. 

calculate_rpm <- function(df, output_file = "df_rpm.xlsx"){
  df_rpm <- df |> mutate(across(-1, function(x) x/sum(x)*10^6, .names = "{.col}_rpm"))
  write.xlsx(df_rpm, output_file)
}

calculate_rpm(
  df = cbe_read_counts, 
  output_file ="CBE_fastq_files/Raw_read_counts/CBE_reads_counts_rpm.xlsx"
)

## Average read counts were calculated across biological replicates for each
## experimental condition. 

cbe_reads_rpm <- read_excel("CBE_fastq_files/Raw_read_counts/CBE_reads_counts_rpm.xlsx") |> clean_names()
names(cbe_reads_rpm)

cbe_reads_rpm_avg <- cbe_reads_rpm |> mutate(
  cbe_bot_rpm_avg = (cbe_b1_rpm + cbe_b2_rpm) / 2, 
  cbe_top_rpm_avg = (cbe_t1_rpm + cbe_t2_rpm) / 2, 
  cbe_unsort_rpm_avg = (cbe_u1_rpm + cbe_u2_rpm) / 2
)

write.xlsx(cbe_reads_rpm_avg, "CBE_fastq_files/Raw_read_counts/CBE_reads_rpm_avg.xlsx")

View(cbe_reads_rpm_avg)

## applying log2 transformation to averaged RPM values: 

log2_transform <- function(df, output_file = "log2_transformed.xlsx"){
  df_log <- df |> mutate(across(ends_with("avg"), function(x) log(x+1, 2), .names = "{.col}_log"))
  write.xlsx(df_log, output_file)
}

log2_transform(
  df = cbe_reads_rpm_avg, 
  output_file = "CBE_fastq_files/Raw_read_counts/CBE_reads_rpm_avg_log.xlsx"
)

cbe_reads_rpm_avg_log <- read_excel("CBE_fastq_files/Raw_read_counts/CBE_reads_rpm_avg_log.xlsx") |> 
  clean_names()
colnames(cbe_reads_rpm_avg_log)

# Calculating log2 fold-change values relative to unsorted population: 

cbe_reads_rpm_avg_log_lfc <- cbe_reads_rpm_avg_log |> 
  mutate(
    cbe_bot_lfc = cbe_bot_rpm_avg_log - cbe_unsort_rpm_avg_log,
    cbe_top_lfc = cbe_top_rpm_avg_log - cbe_unsort_rpm_avg_log
  )

write.xlsx(cbe_reads_rpm_avg_log_lfc, 
           "CBE_fastq_files/Raw_read_counts/CBE_reads_rpm_avg_log_lfc.xlsx")
colnames(cbe_reads_rpm_avg_log_lfc)

## The reference FASTA file containing sgRNA sequences was converted into a tabular format, 
## with each unique identifier linked to its corresponding 20-nt protospacer region. Gene annotations
## (CTNNB1, APC, Axin1, Gsk3b, and control) were assigned based on the experimental design, and
## the resulting dataframe is used for downstream integration with read count and z-score datasets.

cbe_fasta_file <- readDNAStringSet("CBE_fastq_files/Reference_FASTA_FILE/DC_cbe_sgRNAs.fa") |> 
  as.data.frame() |> rownames_to_column(var = "ids")
colnames(cbe_fasta_file) <- c("ids", "sequences")

cbe_fasta_file <- cbe_fasta_file |> 
  dplyr::mutate(gRNA = str_sub(sequences, 6, 25)) |> dplyr::select(1, 3)
cbe_genes <- c("CTNNB1", "APC", "Axin1", "Gsk3b", "control")
cbe_fasta_file <- cbe_fasta_file |> mutate(gene_names = rep(cbe_genes, c(868,2449,1350,446,282)))
colnames(cbe_fasta_file)

cbe_reads_rpm_avg_log_lfc <- left_join(cbe_reads_rpm_avg_log_lfc, cbe_fasta_file, by = "ids") |> 
  relocate(gRNA, gene_names, .after = "ids")

View(cbe_reads_rpm_avg_log_lfc)

## Outliers filtering: To remove extreme outliers from the unsorted population, 
## we inspected the log2_transformed RPM values using a boxplot. Guides with log2(RPM) values
## outside this range defined by 1.5 * interquartile range (IQR) were filtered out. Specifically
##, only guides with log2(RPM) values between 5.53 and 9.26 were retained for calculating 
## z-scores. 


outliers_plot <- ggplot(cbe_reads_rpm_avg_log_lfc, aes(y = cbe_unsort_rpm_avg_log)) + geom_boxplot()
ggplotly(outliers_plot)

cbe_filtered_outliers <- cbe_reads_rpm_avg_log_lfc |> filter(cbe_unsort_rpm_avg_log > 5.53 & cbe_unsort_rpm_avg_log < 9.26)

# Control statistics now:

control_stats <- function(data, gene_colname, columns_to_mean_sd){
  
  control_data <- data |> filter(.data[[gene_colname]] == "control")
  
  control_stats <- control_data |> summarize(across(
    all_of(columns_to_mean_sd), 
    list(mean = function(x) mean(x, na.rm = TRUE),
         sd = function(x) sd(x, na.rm = TRUE))
    
  ))
  return(control_stats)
}

control_stats_cbe <- control_stats(
  data = cbe_filtered_outliers, 
  gene_colname = "gene_names", 
  columns_to_mean_sd = c("cbe_bot_lfc", 
                         "cbe_top_lfc")
)

control_stats_cbe

# Calculating Z-scores:

dc_cbe_zscores <- cbe_filtered_outliers |> mutate(
  cbe_bot_zscores = (cbe_bot_lfc - control_stats_cbe$cbe_bot_lfc_mean)/ control_stats_cbe$cbe_bot_lfc_sd, 
  cbe_top_zscores = (cbe_top_lfc - control_stats_cbe$cbe_top_lfc_mean) / control_stats_cbe$cbe_top_lfc_sd
)

View(dc_cbe_zscores)

## Appending predicted mutations;

dc_cbe_pred <- read_excel("CBE_fastq_files/predicted_mutations_3_10_window/dc_cbe_aa_3_10_distinct.xlsx")
setnames(dc_cbe_pred, old = "sg_rna_sequence", new = "gRNA")
dc_cbe_pred <- dc_cbe_pred[, c(1, 3, 4)]

## Merging predicted mutations with dc_cbe_zscores file:

dc_cbe_zscores_annotations <- left_join(dc_cbe_zscores, dc_cbe_pred, by = "gRNA")
View(dc_cbe_zscores_annotations)

## Saving dc_cbe_zscores and dc_cbe_zscores_annotations:

write.xlsx(dc_cbe_zscores, "CBE_fastq_files/Raw_read_counts/dc_cbe_zscores.xlsx")
write.xlsx(dc_cbe_zscores_annotations, "CBE_fastq_files/Raw_read_counts/dc_cbe_zscores_annotations.xlsx")

## Extracting residue number for z-score vs protein length plots and categorizing
## mutations based on predicted mutation category outcomes: 
categorize_mutations_cbe <- function(data_column){
  data_column <- str_to_lower(data_column)
  case_when(
    str_detect(data_column, "nonsense") ~ "Nonsense",
    str_detect(data_column, "splice-acceptor|splice-donor|intron|utr") ~ "Splice-site", 
    str_detect(data_column, "missense") ~ "Missense",
    TRUE ~ "Silent"
  )
}


dc_cbe_zscores_annotations <- dc_cbe_zscores_annotations |> 
  mutate(residue_number = as.numeric(str_extract(cleaned_mutations, "(?<=\\D)(\\d+)(?=[A-Za-z])")),
         categorized_mutations = categorize_mutations_cbe(mutation_category))

View(dc_cbe_zscores_annotations)
write.xlsx(dc_cbe_zscores_annotations, "CBE_fastq_files/Raw_read_counts/dc_cbe_zscores_annotations.xlsx")

## This final excel sheet contains z-scores, residue numbers and categorized mutations for
## making downstream visualization plots. 














