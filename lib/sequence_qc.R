# ==============================================================================
# Duke Pipeline - Sequence QC Functions
# ==============================================================================
# Functions for sequence quality control, deduplication, and downsampling
# Used primarily in Module 1 (Import and QC)
# ==============================================================================

# Function to find duplicate read names and select the longest
select_longest_duplicate_readname <- function(dna_stringset) {
  
  # Table of read names and lengths
  sequence_table <- data.frame(read_index = seq_along(dna_stringset),
                               read_name = names(dna_stringset),
                               sequence = as.character(dna_stringset),
                               read_length = width(dna_stringset))
  
  # Count duplicate read names
  read_name_count <- sequence_table %>%
    group_by(read_name) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::arrange(-count)
  
  # Select the longest read for each read_name
  longest_read <- sequence_table %>%
    group_by(read_name) %>%
    dplyr::arrange(desc(read_length), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  # Subset the dna_stringset by the indices of the longest reads
  longest_dna_stringset <- dna_stringset[longest_read$read_index]
  
  # Output
  return(list(read_name_count = read_name_count,
              dna_stringset = longest_dna_stringset))
}


# Function to count reads and compare between datasets
compare_read_count <- function(dna_stringset1, dna_stringset2, name1 = NULL, name2 = NULL) {
  
  # Function to count reads
  read_count_function <- function(dna_stringset) {
    data.frame(read_count = unlist(lapply(dna_stringset, length))) %>%
      tibble::rownames_to_column("file_stem")
  }
  
  # Count reads
  read_count1 <- read_count_function(dna_stringset1)
  read_count2 <- read_count_function(dna_stringset2)
  
  # Join
  read_count <- read_count1 %>%
    left_join(read_count2, 
              by = join_by(file_stem),
              suffix = c("1", "2"))
  
  # Rename columns
  if (!is.null(name1)) {
    read_count <- read_count %>%
      dplyr::rename(!!sym(name1) := read_count1)
  }
  if (!is.null(name2)) {
    read_count <- read_count %>%
      dplyr::rename(!!sym(name2) := read_count2)
  }
  
  return(read_count)
}


# Downsample sequences to a specified size
downsample <- function(dna_stringset, sample_size) {
  read_count <- length(dna_stringset)
  if (read_count > sample_size) {
    set.seed(123) # Allow random processes to be reproducible
    dna_stringset <- dna_stringset[sample(seq_along(dna_stringset), sample_size)]
  }
  return(dna_stringset)
}
