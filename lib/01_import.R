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

# ------------------------------------------------------------------------------
# PacBio HiFi detection and per-read aux-tag extraction (v2.3.0)
# ------------------------------------------------------------------------------

# detect_pacbio: Cheap auto-detection of PacBio HiFi BAMs.
#   Primary check: the np:i: tag on the first record (one-record scan).
#   Secondary check (only on failure): @PG header lines with PN matching
#   ccs / skera / lima (case-insensitive).
#   Returns FALSE for non-BAM formats, missing files, and any error.
#
#   Usage:
#     is_pacbio <- detect_pacbio(path = "/path/to/sample.bam", format = "bam")
detect_pacbio <- function(path, format) {

  # Non-BAM inputs are never PacBio
  if (!is.character(format) || length(format) != 1 || tolower(format) != "bam") {
    return(FALSE)
  }
  if (!is.character(path) || length(path) != 1 || !file.exists(path)) {
    return(FALSE)
  }

  tryCatch({

    # Primary: read the np:i: tag from the first record
    bf <- BamFile(path, yieldSize = 1)
    open(bf)
    on.exit(close(bf), add = TRUE)

    first_record <- scanBam(bf, param = ScanBamParam(tag = "np"))
    np_tag <- first_record[[1]]$tag$np
    if (!is.null(np_tag) && length(np_tag) > 0 && !all(is.na(np_tag))) {
      return(TRUE)
    }

    # Secondary: scan @PG header lines for ccs / skera / lima
    hdr <- scanBamHeader(path)
    pg_block <- hdr[[1]]$text
    if (!is.null(pg_block) && length(pg_block) > 0) {
      pg_entries <- pg_block[names(pg_block) == "@PG"]
      if (length(pg_entries) > 0) {
        pg_flat <- unlist(pg_entries, use.names = FALSE)
        pn_values <- grep("^PN:", pg_flat, value = TRUE, ignore.case = TRUE)
        if (length(pn_values) > 0) {
          pn_clean <- sub("^PN:", "", pn_values, ignore.case = TRUE)
          if (any(grepl("^(ccs|skera|lima)$", pn_clean, ignore.case = TRUE))) {
            return(TRUE)
          }
        }
      }
    }

    FALSE

  }, error = function(e) FALSE)
}


# extract_pacbio_metadata: Single-pass extraction of PacBio HiFi per-read aux
#   tags from a BAM. Returns a data.frame with one row per read and columns
#   read_name, np, rq, ec, zm, bq. Absent tags are NA.
#
#   Usage:
#     md <- extract_pacbio_metadata(path = "/path/to/sample.bam")
extract_pacbio_metadata <- function(path) {

  bf <- BamFile(path)
  open(bf)
  on.exit(close(bf), add = TRUE)

  bam_data <- scanBam(
    bf,
    param = ScanBamParam(
      what = "qname",
      tag  = c("np", "rq", "ec", "zm", "bq")
    )
  )

  qnames <- bam_data[[1]]$qname
  tags   <- bam_data[[1]]$tag
  n_reads <- length(qnames)

  # Pull a tag with NA padding when absent or wrong length
  pull_tag <- function(tag_name, expected_type = NA_real_) {
    v <- tags[[tag_name]]
    if (is.null(v) || length(v) == 0) {
      return(rep(expected_type, n_reads))
    }
    if (length(v) != n_reads) {
      # Tag present on some records only — re-fill to length with NA padding
      warning(sprintf(
        "Tag '%s' length (%d) does not match read count (%d); padding with NA",
        tag_name, length(v), n_reads))
      return(rep(expected_type, n_reads))
    }
    v
  }

  data.frame(
    read_name = qnames,
    np = pull_tag("np", NA_integer_),
    rq = pull_tag("rq", NA_real_),
    ec = pull_tag("ec", NA_real_),
    zm = pull_tag("zm", NA_integer_),
    bq = pull_tag("bq", NA_integer_),
    stringsAsFactors = FALSE
  )
}
