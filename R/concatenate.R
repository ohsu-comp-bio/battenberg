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
  # Robust filename resolution: try both '1' and 'chr1'
  find_file <- function(prefix, chrom, suffix) {
    f1 <- paste0(prefix, chrom, suffix)
    if (file.exists(f1)) {
      return(f1)
    }
    # Try with/without 'chr'
    if (grepl("^chr", chrom, ignore.case = TRUE)) {
      f2 <- paste0(prefix, gsub("^chr", "", chrom, ignore.case = TRUE), suffix)
    } else {
      f2 <- paste0(prefix, "chr", chrom, suffix)
    }
    if (file.exists(f2)) {
      return(f2)
    }
    return(NULL)
  }

  infiles <- character(0)
  for (cn in chr_names) {
    f <- find_file(inputStart, cn, inputEnd)
    if (!is.null(f) && file.info(f)$size > 0) {
      infiles <- c(infiles, f)
    }
  }

  if (length(infiles) == 0) {
    return(data.frame())
  }
  log_info("Using {length(infiles)} infiles in concatenateAlleleCountFiles. Example: {infiles[1]}")

  # Bulk read using vroom. We remove delim="\t" to allow guessing,
  # which handles both space and tab delimited counts.
  combined <- vroom::vroom(
    infiles,
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
  # Robust filename resolution
  find_file <- function(prefix, chrom, suffix) {
    f1 <- paste0(prefix, chrom, suffix)
    if (file.exists(f1)) {
      return(f1)
    }
    if (grepl("^chr", chrom, ignore.case = TRUE)) {
      f2 <- paste0(prefix, gsub("^chr", "", chrom, ignore.case = TRUE), suffix)
    } else {
      f2 <- paste0(prefix, "chr", chrom, suffix)
    }
    if (file.exists(f2)) {
      return(f2)
    }
    return(NULL)
  }

  existing_files <- character(0)
  for (cn in chr_names) {
    f <- find_file(inputStart, cn, inputEnd)
    if (!is.null(f) && file.info(f)$size > 0) {
      existing_files[cn] <- f
    }
  }

  if (length(existing_files) == 0) {
    return(data.frame())
  }

  # Read files individually to inject chromosome if missing (common in some bundles)
  # using data.table::fread for multi-delimiter robustness
  datalist <- lapply(names(existing_files), function(cn) {
    f <- existing_files[cn]

    # Force colClasses to character for initial read to prevent parsing issues
    d <- data.table::fread(f, sep = "auto", header = "auto", colClasses = "character", data.table = FALSE)

    if (ncol(d) == 3) {
      # File has (POS, A0, A1), we prepend the CHR from filename
      d <- cbind(CHR = cn, d)
    }

    # Ensure consistent column naming to prevent binding issues
    colnames(d)[1:4] <- c("CHR", "POS", "A0", "A1")

    # Standardise structure to exactly 4 columns: CHR, POS, A0, A1
    return(d[, 1:4])
  })

  combined <- data.table::as.data.table(data.table::rbindlist(datalist, use.names = TRUE))
  return(combined)
}
