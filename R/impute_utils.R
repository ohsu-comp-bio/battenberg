#' Read in the reference_info_file.
#'
#' Reads in a file with the following columns:
#'   chromosome : 1-X
#'   impute_legend : Legend file in IMPUTE -l format
#'   genetic_map : Genetic map file in IMPUTE -m format
#'   impute_hap : Phased haplotype file in IMPUTE -h format
#'   start : Start of the chromosome
#'   end : End of the chromosome
#'   is_par : 1 when pseudo autosomal region, 0 when not
#'
#' @param reference_info_file Path to the reference_info_file on disk.
#' @param is_male A boolean describing whether the sample under study is male.
#' @param chrom The name of a chromosome to subset the contents of the reference_info_file with (optional)
#' @return A data.frame with 7 columns: Chromosome, impute_legend, genetic_map, impute_hap, start, end, is_par
#' @author sd11
#' @export
parse_imputeinfofile <- function(reference_info_file, is_male, chrom = NA) {
    if (is.na(reference_info_file) || !file.exists(reference_info_file)) {
        return(data.table::data.table())
    }

    # Use fread for high-speed reading.
    impute_info <- data.table::fread(
        reference_info_file,
        col.names = c(
            "chrom", "impute_legend", "genetic_map",
            "impute_hap", "start", "end", "is_par"
        ),
        stringsAsFactors = FALSE
    )

    expected_cols <- c("chrom", "impute_legend", "genetic_map", "impute_hap", "start", "end", "is_par")
    if (!all(expected_cols %in% names(impute_info))) {
        # If columns are missing, try to assign them if possible, or fail
        if (ncol(impute_info) == length(expected_cols)) {
            names(impute_info) <- expected_cols
        } else {
            log_failure("Reference info file does not have the expected number of columns (7). Found: {ncol(impute_info)}")
        }
    }

    # Filter based on gender
    if (!is.na(is_male) && !is_male) {
        impute_info <- impute_info[impute_info[["chrom"]] != "Y", ]
    }
    # Subset for a particular chromosome
    if (!is.na(chrom)) {
        impute_info <- impute_info[impute_info[["chrom"]] == chrom, ]
    }
    return(impute_info)
}

#' Check reference info file consistency
#' @param reference_info_file Path to the reference_info_file on disk.
#' @author sd11
check_imputeinfofile <- function(reference_info_file, is_male, usebeagle) {
    if (is.na(reference_info_file)) {
        return(invisible(NULL))
    }

    impute_info <- parse_imputeinfofile(reference_info_file, is_male)
    if (nrow(impute_info) == 0) {
        return(invisible(NULL))
    }

    if (usebeagle) {
        # For Beagle input generation, we only strictly need the legend file
        if (any(!file.exists(as.character(impute_info$impute_legend)))) {
            log_failure("Could not find reference legend files, make sure paths in reference_info_file point to the correct location")
        }
    } else {
        if (any(!file.exists(as.character(impute_info$impute_legend)) |
            !file.exists(as.character(impute_info$genetic_map)) |
            !file.exists(as.character(impute_info$impute_hap)))) {
            log_failure("Could not find reference files, make sure paths in reference_info_file point to the correct location")
        }
    }
}

#' Returns the chromosome names that are supported
#' @param chrom_names A vector of chromosome names to use directly (optional)
#' @return A vector containing the supported chromosome names
#' @author sd11
#' @export
get_chrom_names <- function(reference_info_file = NA, is_male = NA, chrom = NA, analysis = "paired", chrom_names = NULL,
                            usebeagle = FALSE, beagleref_dir = NA) {
    if (!is.null(chrom_names)) {
        return(chrom_names)
    }

    if (is.na(reference_info_file)) {
        # If we are using Beagle, we might be able to infer chroms from beagleref_dir
        if (usebeagle && !is.na(beagleref_dir) && dir.exists(beagleref_dir)) {
            vcfs <- list.files(beagleref_dir, pattern = "\\.vcf(\\.gz)?$")
            found_chroms <- gsub(".*chr([0-9XY]+).*", "\\1", vcfs)
            found_chroms <- unique(found_chroms[found_chroms %in% c(as.character(1:22), "X", "Y")])
            if (length(found_chroms) > 0) {
                log_info("Inferred chromosomes from Beagle reference directory: {paste(found_chroms, collapse=', ')}")
                return(sort(found_chroms))
            }
        }
        # Fallback to standard human autosomes if nothing else provided
        log_warning("No reference_info_file or chrom_names provided. Defaulting to 1-22.")
        return(as.character(1:22))
    }

    chrom_names <- unique(parse_imputeinfofile(reference_info_file, is_male, chrom = chrom)$chrom)
    if (analysis == "cell_line" || analysis == "germline") {
        # Both cell line and germline analysis do not yield usable data on X and Y, so remove
        chrom_names <- chrom_names[!chrom_names %in% c("X", "Y")]
    }
    return(chrom_names)
}

#' Concatenate the impute output generated for each of the regions.
#'
#' This function assembles the impute output generated.
#' @param inputfile.prefix Prefix of the input files.
#' @param outputfile Where to store the output.
#' @param is_male Boolean describing whether the sample is male (TRUE) or female (FALSE).
#' @param reference_info_file Path to the reference_info_file on disk.
#' @param region.size An integer describing the region size to be used by impute (optional).
#' @param chrom The name of a chromosome on which this function should run.
#' @author dw9
#' @export
combine_impute_output <- function(inputfile.prefix, outputfile, is_male, reference_info_file, region.size = 5000000, chrom = NA) {
    # Read in the impute file information
    impute_info <- parse_imputeinfofile(reference_info_file, is_male, chrom = chrom)

    # Assemble the start and end points of all regions
    all.boundaries <- array(0, c(0, 2))
    for (r in seq_len(nrow(impute_info))) {
        boundaries <- seq(as.numeric(impute_info[r, ]$start), as.numeric(impute_info[r, ]$end), region.size)
        if (boundaries[length(boundaries)] != impute_info[r, ]$end) {
            boundaries <- c(boundaries, impute_info[r, ]$end)
        }
        all.boundaries <- rbind(all.boundaries, cbind(boundaries[-(length(boundaries))], boundaries[-1]))
    }
    # Concatenate all the regions
    impute.output <- concatenateImputeFiles(inputfile.prefix, all.boundaries)
    data.table::fwrite(
        impute.output,
        file = outputfile,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = " "
    )
}

#' Load centromere coordinates from a reference file
#'
#' @param coord_file Path to the gcCorrect_chromosome_coordinates_hg38.txt or similar file.
#' @return A named list of centromere split points.
#' @keywords internal
load_centromere_splits <- function(coord_file) {
    if (!file.exists(coord_file)) {
        log_failure("Centromere coordinate file not found: {coord_file}")
    }
    coords <- data.table::fread(coord_file, header = TRUE)
    # Map columns (chr, cen.left.base, cen.right.base) to a single split point (mean)
    splits <- list()
    for (i in seq_len(nrow(coords))) {
        chr <- as.character(coords$chr[i])
        # split point is the middle of the centromere range
        splits[[chr]] <- (coords$cen.left.base[i] + coords$cen.right.base[i]) / 2
    }
    return(splits)
}
