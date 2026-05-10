# Functions to process Fastq files and align them against the reference,
# sort and index bam files, and get readcounts for each library: 
library(Rsubread)
library(Biostrings)
library(GenomicAlignments)
library(GenomicFeatures)
library(QuasR)
library(Rsamtools)
library(ShortRead)
library(readxl)
library(openxlsx)
library(stringr)
library(dplyr)
library(data.table)


# Individual functions: 

## Step 1: setup function with directories: 

setup_grna_pipeline <- function(input_dir){
  
  processed_dir <- file.path(input_dir, "processed_reads")
  bam_dir <- file.path(input_dir, "sorted_bam_files")
  counts_dir <- file.path(input_dir, "read_counts")
  
  dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(bam_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(counts_dir, showWarnings = FALSE, recursive = TRUE)
  
  return(list(
    processed_dir = processed_dir,
    bam_dir = bam_dir,
    counts_dir = counts_dir
  ))
}

# The key arguments for FASTQ processing include truncateEndBases, 
# which allows trimming of a specified number of bases from the 3′ end of reads, 
# and truncateStartBases, which performs trimming from the 5′ end. In addition, 
# pattern-based trimming can be applied using Lpattern and Rpattern, where defined 
# sequence motifs are used to identify positions within the read for more targeted 
# processing at the 5′ and 3′ ends, respectively. These parameters can be flexibly 
# adjusted depending on the design of the sequencing library and the specific 
# requirements of the experimental setup.

process_fastq_files <- function(
    input_dir,
    output_suffix = "processed.fastq.gz",
    nBases = 2,
    Lpattern = "",
    Rpattern = "GAGAGCTAGAAATAGCAAGTT",
    truncateEndBases = 127,
    truncateStartBases = NULL
){
  
  processed_dir <- file.path(input_dir, "processed_reads")
  
  fastq_files <- list.files(
    input_dir,
    pattern = "\\.fastq\\.gz$",
    full.names = TRUE
  )
  
  message("Found ", length(fastq_files), " FASTQ files.")
  
  for(file in fastq_files){
    
    output_file <- file.path(
      processed_dir,
      sub("\\.fastq\\.gz$", output_suffix, basename(file))
    )
    
    message("Processing: ", basename(file))
    
    preprocessReads(
      file,
      output_file,
      nBases = nBases,
      Rpattern = Rpattern,
      Lpattern = Lpattern,
      truncateEndBases = truncateEndBases,
      truncateStartBases = truncateStartBases
    )
  }
  
  message("FASTQ preprocessing completed.")
  
  return(processed_dir)
}


## Step 3: Index building; The reference file used for index construction 
## should be in FASTA format and contain two columns: guide RNA identifiers 
## and their corresponding gRNA sequences as used in the library.

build_grna_index <- function(
    reference_file,
    index_name = "grna_library",
    memory = 8000
    ){
      
      buildindex(
        basename = index_name,
        reference = reference_file,
        memory = memory
      )
      
      message("Index built successfully")
      
      return(index_name)
}

## Step 4: Alignment: Aligning each processed FASTQ file against the reference 
## index to generate read counts for downstream analysis.

align_processed_fastqs <- function(
    processed_dir,
    bam_dir,
    index_name
){
  
  processed_fastqs <- list.files(
    processed_dir,
    pattern = "processed.fastq.gz$",
    full.names = TRUE
  )
  
  message("Found ", length(processed_fastqs), " processed FASTQ files.")
  
  for(fastq_file in processed_fastqs){
    
    base_name <- tools::file_path_sans_ext(
      tools::file_path_sans_ext(basename(fastq_file))
    )
    
    bam_output <- file.path(
      bam_dir,
      paste0(base_name, ".bam")
    )
    
    message("Aligning: ", basename(fastq_file))
    
    align(
      index = index_name,
      readfile1 = fastq_file,
      output_file = bam_output,
      output_format = "BAM"
    )
  }
  
  message("Alignment completed.")
}


# BAM sorting and indexing: 

sort_and_index_bams <- function(bam_dir){
  
  sorted_dir <- file.path(bam_dir, "sorted_bams")
  
  if(!dir.exists(sorted_dir)){
    dir.create(sorted_dir, recursive = TRUE)
  }
  
  bam_files <- list.files(
    bam_dir,
    pattern = "\\.bam$",
    full.names = TRUE
  )
  
  for(bam_file in bam_files){
    
    base_name <- tools::file_path_sans_ext(basename(bam_file))
    
    sorted_prefix <- file.path(
      sorted_dir,
      paste0("sorted_", base_name)
    )
    
    message("Sorting: ", basename(bam_file))
    
    sortBam(
      bam_file,
      destination = sorted_prefix
    )
  }
  
  return(sorted_dir)
}


## Step 6: Counting: This step generates a read count table from the aligned data,
## summarizing the number of reads mapped to each guide RNA or reference sequence 
## across all samples. The resulting matrix serves as the primary input for downstream 
## normalization and statistical analysis.

count_grna_reads <- function(
    bam_dir,
    counts_dir,
    sample_regex = "(?<=sorted_).*(?=\\.bam)"
){
  
  sorted_bams <- list.files(
    bam_dir,
    pattern = "^sorted_.*\\.bam$",
    full.names = TRUE
  )
  
  count_tables <- list()
  
  for(i in seq_along(sorted_bams)){
    
    bam_file <- sorted_bams[i]
    
    message("Counting: ", basename(bam_file))
    
    alignments <- readGAlignments(bam_file)
    
    seq_table <- table(seqnames(alignments))
    
    count_df <- data.frame(
      ids = names(seq_table),
      counts = as.numeric(seq_table)
    )
    
    sample_name <- stringr::str_extract(
      basename(bam_file),
      sample_regex
    )
    
    colnames(count_df)[2] <- sample_name
    
    output_excel <- file.path(
      counts_dir,
      paste0(sample_name, ".xlsx")
    )
    
    write.xlsx(count_df, output_excel)
    
    count_tables[[i]] <- count_df
  }
  
  return(count_tables)
}

## Step 7: Merge counts: This step combines read count tables from all libraries 
## into a single consolidated Excel file, enabling unified downstream analysis 
## across all conditions and replicates.

merge_count_tables <- function(
    count_tables,
    counts_dir,
    output_file = "combined_read_counts.xlsx"
){
  
  final_data <- Reduce(
    function(x, y) merge(x, y, by = "ids", all = TRUE),
    count_tables
  )
  
  final_data[is.na(final_data)] <- 0
  
  final_data$order_num <- as.numeric(gsub(".*_(\\d+)$", "\\1", final_data$ids))
  final_data <- final_data[order(final_data$order_num), ]
  final_data$order_num <- NULL
  
  final_path <- file.path(counts_dir, output_file)
  
  write.xlsx(final_data, final_path)
  
  message("Combined count table saved at: ", final_path)
  
  return(final_data)
}


## The run_grna_pipeline function provides a fully automated workflow for 
## processing CRISPR gRNA screening data from raw FASTQ files to final read count 
## tables. It orchestrates the entire analysis in sequential steps, beginning with 
## directory setup, followed by FASTQ preprocessing that includes trimming and 
## pattern-based filtering. Next, it builds a gRNA reference index from the 
## provided library file and aligns processed reads against this index. 
## The resulting alignments are then sorted and indexed in BAM format to 
## ensure efficient downstream access. Read counts are subsequently generated 
## for each sample and finally merged into a single consolidated output table. 
## Overall, this wrapper function integrates all individual modules into a 
## streamlined, end-to-end pipeline for gRNA screen analysis.

run_grna_pipeline <- function(
    input_dir,
    reference_file,
    index_name = "grna_library",
    nBases = 2,
    Lpattern = "",
    Rpattern = "GAGAGCTAGAAATAGCAAGTT",
    truncateEndBases = 127,
    truncateStartBases = NULL
){
   # Step1 function: Setting up directories
  
  dirs <- setup_grna_pipeline(input_dir)
  
  processed_dir <- dirs$processed_dir
  bam_dir <- dirs$bam_dir
  counts_dir <- dirs$counts_dir
  
  # Step2 function: Fastq processing
  process_fastq_files(
    input_dir = input_dir,
    nBases = nBases,
    Lpattern = Lpattern,
    Rpattern = Rpattern,
    truncateEndBases = truncateEndBases
  )
  # Step3 function: Building index
  build_grna_index(
    reference_file = reference_file,
    index_name = index_name
  )
  # Step4 function: Alignment against reference
  align_processed_fastqs(
    processed_dir = processed_dir,
    bam_dir = bam_dir,
    index_name = index_name
  )
  #Step5 function: sorting and indexing bam files
  sorted_bam_dir <- sort_and_index_bams(bam_dir)
  
  #Step6 function: Generating readcounts
  count_tables <- count_grna_reads(
    bam_dir = sorted_bam_dir,
    counts_dir = counts_dir
  )
  #Step7: Merging count tables
  final_counts <- merge_count_tables(
    count_tables = count_tables,
    counts_dir = counts_dir
  )
  
  message("Pipeline completed successfully.")
  
  return(final_counts)
}



