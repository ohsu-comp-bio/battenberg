########################################################################################
# Concatenate files
########################################################################################
#' Function to concatenate Impute output
#' @noRd
concatenateImputeFiles <- function(inputStart, boundaries) {
  # Generate the list of potential filenames
  # Using paste0 and vectorized division for a bit more speed
  infiles <- paste0(inputStart, "_", boundaries[, 1] / 1000, "K_", boundaries[, 2] / 1000, "K.txt_haps")

  # Filter for existing files with data
  existing_files <- infiles[file.exists(infiles) & file.info(infiles)$size > 0]

  if (length(existing_files) == 0) {
    return(NULL)
  }
  # Impute files (.haps) have no headers
  result <- vroom::vroom(
    existing_files,
    delim = " ",
    col_names = FALSE,
    show_col_types = FALSE
  )
  return(data.table::as.data.table(result))
}

#' Function to concatenate allele counter output
#' @noRd
concatenateAlleleCountFiles <- function(inputStart, inputEnd, chr_names) {
  # Vectorized filename generation
  all_files <- paste0(inputStart, chr_names, inputEnd)

  # Vectorized file checking (much faster than a for-loop)
  # This filters the list to only existing, non-empty files
  infiles <- all_files[file.exists(all_files) & file.info(all_files)$size > 0]
  if (length(infiles) == 0) {
    return(data.frame())
  }
  log_info("Using {length(infiles)} infiles in concatenateAlleleCountFiles. Example: {infiles[1]}")

  # Bulk read using vroom for significant speedup
  # Allele counter files typically have no header or start with '#' comments
  combined <- vroom::vroom(
    infiles,
    delim = "\t",
    col_names = c("CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth"),
    col_types = "ciiiiii",
    comment = "#",
    show_col_types = FALSE
  )
  data.table::setDF(combined)
  return(combined)
}

#' Function to concatenate 1000 Genomes SNP reference files
#' @noRd
concatenateG1000SnpFiles <- function(inputStart, inputEnd, chr_names) {
  # Vectorized filename generation
  filenames <- paste0(inputStart, chr_names, inputEnd)
  names(filenames) <- chr_names

  # Filter for valid files
  existing_files <- filenames[file.exists(filenames) & file.info(filenames)$size > 0]

  if (length(existing_files) == 0) {
    return(data.frame())
  }

  # Bulk read using vroom for speed
  # Reference files have a header
  combined <- vroom::vroom(
    existing_files,
    delim = "\t",
    col_types = vroom::cols(.default = "c"),
    show_col_types = FALSE
  )

  data.table::setDF(combined)
  return(combined)
}
