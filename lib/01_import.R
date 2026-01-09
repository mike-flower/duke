# ==============================================================================
# Duke Pipeline - Import Functions
# ==============================================================================
# Functions for importing FASTQ, FASTA, and BAM files
# ==============================================================================

# Function to import fastq including quality data
fastq_import <- function(path) {
  readfastq_file <- readFastq(path)
  read_names <- as.character(ShortRead::id(readfastq_file))
  # Strip tags from read names by keeping only the first field before the first space
  clean_read_names <- sapply(strsplit(read_names, " "), `[`, 1)
  quality_scaled_dnastringset <- 
      QualityScaledDNAStringSet(sread(readfastq_file), as(quality(readfastq_file), "PhredQuality"))
  names(quality_scaled_dnastringset) <- clean_read_names
  return(quality_scaled_dnastringset)
}

# Function to import fasta
fasta_import <- function(path) {
  readDNAStringSet(path, format = "fasta")
}

# Function to import bam
bam_import <- function(path) {
  bam_file <- BamFile(path)
  bam_data <- scanBam(bam_file)
  bam_sequences <- bam_data[[1]]$seq
  bam_qualities <- bam_data[[1]]$qual
  bam_read_names <- bam_data[[1]]$qname
  quality_scaled_dnastringset <-
    QualityScaledDNAStringSet(bam_sequences, bam_qualities)
  names(quality_scaled_dnastringset) <- bam_read_names
  return(quality_scaled_dnastringset)
}

# Function to import sequencing files
import_sequencing <- function(path, format) {
  if (format == "fastq" || format == "fastq.gz") {
    return(fastq_import(path))
  } else if (format == "fasta" || format == "fa") {
    return(fasta_import(path))
  } else if (format == "bam") {
    return(bam_import(path))
  } else {
    stop("Unsupported file format")
  }
}
