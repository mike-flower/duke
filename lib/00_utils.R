# ==============================================================================
# Duke Pipeline - Utility Functions
# ==============================================================================
# General utilities: Excel export, threading wrappers
# ==============================================================================

# Excel export function
write_excel_function <- function(named_list, path, overwrite_sheet = TRUE) {
  if (file.exists(path)) { wb <- loadWorkbook(path) } else { wb <- createWorkbook() }
  
  truncate_name <- function(name, max_length = 31) {
    if (nchar(name) > max_length) {
      return(paste0(substr(name, 1, max_length - 2), ".."))
    } else {
      return(name)
    }
  }
  
  for (table_name in names(named_list)) {
    truncated_name <- truncate_name(table_name)
    if (truncated_name %in% names(wb)) {
      if (overwrite_sheet) {
        removeWorksheet(wb, truncated_name)
        export_name <- truncated_name
      } else {
        suffix <- 1
        new_name <- paste0(truncated_name, "_", suffix)
        while (new_name %in% names(wb)) {
          suffix <- suffix + 1
          new_name <- paste0(truncated_name, "_", suffix)
        }
        export_name <- new_name
      }
    } else {
      export_name <- truncated_name
    }
    addWorksheet(wb, export_name)
    writeData(wb, export_name, named_list[[table_name]])
  }
  saveWorkbook(wb, file = path, overwrite = TRUE)
}

# Threading wrapper
apply_fn <- function(...) {
  if (!is.na(params$threads) && params$threads > 1) {
    pbmclapply(..., mc.cores = params$threads)
  } else {
    pblapply(...)
  }
}

# ------------------------------------------------------------------------------
# Timing utilities
# ------------------------------------------------------------------------------

# Format a byte count as a human-readable string with appropriate units
format_size <- function(bytes) {
  bytes <- as.numeric(bytes)
  if (bytes < 1024) {
    sprintf("%.0f B", bytes)
  } else if (bytes < 1024^2) {
    sprintf("%.1f KB", bytes / 1024)
  } else if (bytes < 1024^3) {
    sprintf("%.1f MB", bytes / 1024^2)
  } else {
    sprintf("%.2f GB", bytes / 1024^3)
  }
}

# Format elapsed seconds as mm:ss or hh:mm:ss
format_elapsed <- function(secs) {
  secs <- round(as.numeric(secs))
  if (secs < 3600) {
    sprintf("%02d:%02d", secs %/% 60, secs %% 60)
  } else {
    sprintf("%02d:%02d:%02d", secs %/% 3600, (secs %% 3600) %/% 60, secs %% 60)
  }
}

# Build a timing summary table from a named numeric vector of elapsed seconds
build_timing_table <- function(timings) {
  elapsed_secs <- as.numeric(timings)
  total <- sum(elapsed_secs)
  df <- data.frame(
    section    = c(names(timings), "Total"),
    secs       = c(round(elapsed_secs, 1), round(total, 1)),
    elapsed    = c(sapply(elapsed_secs, format_elapsed), format_elapsed(total)),
    cumulative = c(sapply(cumsum(elapsed_secs), format_elapsed), format_elapsed(total)),
    stringsAsFactors = FALSE
  )
  df
}

# ------------------------------------------------------------------------------
# Safe cache loading utilities
# ------------------------------------------------------------------------------

# safe_load_cache: Load a temp/ cache file with automatic fallback on corruption.
#   Returns TRUE if the file loaded successfully, FALSE if it was missing or
#   corrupt (in which case the corrupt file is deleted so the caller recomputes).
#   Call with envir = parent.env(environment()) so loaded objects are visible
#   in the calling scope.
#
#   Usage pattern:
#     if (params$resume && safe_load_cache(cache_path, environment())) {
#       message("Loaded from cache")
#     } else {
#       # recompute ...
#       save(result, file = cache_path)
#     }
safe_load_cache <- function(path, envir = parent.frame()) {
  if (!file.exists(path)) return(FALSE)
  if (file.size(path) == 0) {
    warning("Cache file is empty (zero bytes), will recompute: ", basename(path))
    file.remove(path)
    return(FALSE)
  }
  tryCatch({
    load(path, envir = envir)
    TRUE
  }, error = function(e) {
    warning("Cache file is corrupt, will recompute: ", basename(path))
    file.remove(path)
    FALSE
  })
}

# safe_load_module: Load a module_data/ RData file with a clear error message
#   if the file is missing, empty, or corrupt. These files cannot be recovered
#   within the current module — the user must delete the file and rerun the
#   source module.
#
#   Usage:
#     safe_load_module(module3_path, module_num = 3, envir = environment())
safe_load_module <- function(path, module_num, envir = parent.frame()) {
  if (!file.exists(path)) {
    stop("Module ", module_num, " output not found. Has Module ", module_num,
         " been run?\n  Expected: ", path)
  }
  if (file.size(path) == 0) {
    stop("Module ", module_num, " output is empty (zero bytes) — the previous run",
         " was likely interrupted during save.\n",
         "  Delete and rerun Module ", module_num, ":\n  ", path)
  }
  tryCatch(
    load(path, envir = envir),
    error = function(e) {
      stop("Module ", module_num, " output is corrupt — the previous run was",
           " likely interrupted during save.\n",
           "  Delete and rerun Module ", module_num, ":\n  ", path, "\n",
           "  Original error: ", e$message)
    }
  )
}

# Extract result
extract_apply_fn_result <- function(result) {
  if ("value" %in% names(result)) {
    if (!is.null(result$warning)) {
      warning(result$warning$message)
    }
    result <- result$value
  }
  return(result)
}
