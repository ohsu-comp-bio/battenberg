#' Read in the imputeinfofile.
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
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @param is_male A boolean describing whether the sample under study is male.
#' @param chrom The name of a chromosome to subset the contents of the imputeinfofile with (optional)
#' @return A data.frame with 7 columns: Chromosome, impute_legend, genetic_map, impute_hap, start, end, is_par
#' @author sd11
#' @export
parse_imputeinfofile <- function(imputeinfofile, is_male, chrom = NA) {
  # Use fread for high-speed reading.
  impute_info <- data.table::fread(
    imputeinfofile,
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
      log_failure("Impute info file does not have the expected number of columns (7). Found: {ncol(impute_info)}")
    }
  }

  # Filter based on gender
  if (!is.na(is_male) && !is_male) {
    # If female, we exclude Y chromosome regions
    # and we might want to handle PAR specifically if the pipeline requires it.
    # But generally, we just want to ensure we don't return Y.
    impute_info <- impute_info[impute_info[["chrom"]] != "Y", ]
  }
  # Subset for a particular chromosome
  if (!is.na(chrom)) {
    impute_info <- impute_info[impute_info[["chrom"]] == chrom, ]
  }
  return(impute_info)
}

#' Check impute info file consistency
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @author sd11
check_imputeinfofile <- function(imputeinfofile, is_male, usebeagle) {
  impute_info <- parse_imputeinfofile(imputeinfofile, is_male)
  if (usebeagle) {
    if (any(!file.exists(impute_info$impute_legend))) {
      log_failure("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
    }
  } else {
    if (any(!file.exists(impute_info$impute_legend) | !file.exists(impute_info$genetic_map) | !file.exists(impute_info$impute_hap))) {
      log_failure("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
    }
  }
}

#' Returns the chromosome names that are supported
#' @param chrom_names A vector of chromosome names to use directly (optional)
#' @return A vector containing the supported chromosome names
#' @author sd11
#' @export
get_chrom_names <- function(imputeinfofile = NA, is_male = NA, chrom = NA, analysis = "paired", chrom_names = NULL) {
  if (!is.null(chrom_names)) {
    return(chrom_names)
  }

  if (is.na(imputeinfofile)) {
    # Fallback to standard human autosomes if nothing else provided
    log_warning("No imputeinfofile or chrom_names provided. Defaulting to 1-22.")
    return(as.character(1:22))
  }

  chrom_names <- unique(parse_imputeinfofile(imputeinfofile, is_male, chrom = chrom)$chrom)
  if (analysis == "cell_line" || analysis == "germline") {
    # Both cell line and germline analysis do not yield usable data on X and Y, so remove
    chrom_names <- chrom_names[!chrom_names %in% c("X", "Y")]
  }
  return(chrom_names)
}

#' Concatenate the impute output generated for each of the regions.
#'
#' This function assembles the impute output generated.
#' @param inputfile.prefix Prefix of the input files (this is typically the outputfile_prefix option supplied when calling run_impute).
#' @param outputfile Where to store the output.
#' @param is_male Boolean describing whether the sample is male (TRUE) or female (FALSE).
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @param region.size An integer describing the region size to be used by impute (optional).
#' @param chrom The name of a chromosome on which this function should run (names are used, supply X as 'X').
#' @author dw9
#' @export
combine_impute_output <- function(inputfile.prefix, outputfile, is_male, imputeinfofile, region.size = 5000000, chrom = NA) {
  # Read in the impute file information
  impute_info <- parse_imputeinfofile(imputeinfofile, is_male, chrom = chrom)

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


#' @export
convert_impute_input_to_beagle_vcf <- function(impute_input_data, chrom) {
  # Standardize chrom for VCF
  chr_vcf <- if (chrom == "23") "X" else as.character(chrom)

  # Column mapping: Battenberg intermediate format to VCF
  # VCF Columns: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
  coln <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMP001")

  # Battenberg intermediate (impute_input) columns (standardized by read_impute_input):
  # X1: snpID, X2: Chr, X3: Pos, X4: Ref, X5: Alt, X6: HomRef, X7: Het, X8: HomAlt
  vcf <- data.frame(
    CHROM = rep(chr_vcf, nrow(impute_input_data)),
    POS = impute_input_data$X3,
    ID = rep(".", nrow(impute_input_data)),
    REF = impute_input_data$X4,
    ALT = impute_input_data$X5,
    QUAL = rep(".", nrow(impute_input_data)),
    FILTER = rep("PASS", nrow(impute_input_data)),
    INFO = rep(".", nrow(impute_input_data)),
    FORMAT = rep("GT", nrow(impute_input_data)),
    GT = paste(impute_input_data$X6, impute_input_data$X7, impute_input_data$X8, sep = "-"),
    stringsAsFactors = FALSE
  )

  # Convert 1-hot encoding to VCF GT format (0/0, 0/1, 1/1)
  vcf$GT[vcf$GT == "1-0-0"] <- "0/0"
  vcf$GT[vcf$GT == "0-1-0"] <- "0/1"
  vcf$GT[vcf$GT == "0-0-1"] <- "1/1"

  # Filter samples with no valid genotype (0-0-0)
  vcf <- vcf[vcf$GT != "0-0-0", ]

  colnames(vcf) <- coln
  return(vcf)
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

#' Split a VCF into p and q arms
#'
#' @param vcf A data.table containing VCF data.
#' @param chrom Chromosome name.
#' @param pathP Path to write p-arm VCF.
#' @param pathQ Path to write q-arm VCF.
#' @param coord_file Path to chromosome coordinates file.
#' @export
split_and_writevcf_by_arm <- function(vcf, chrom, pathP, pathQ, coord_file) {
  centromere_split <- load_centromere_splits(coord_file)

  # Standardize chrom name for lookup
  lookup_chrom <- if (chrom == "X") "23" else as.character(chrom)

  if (!(lookup_chrom %in% names(centromere_split))) {
    log_warning("Chromosome '{chrom}' not found in centromere table. Phasing as single unit.")
    writevcf_beagle(vcf, pathP)
    return(invisible(NULL))
  }

  split_point <- centromere_split[[lookup_chrom]]
  vcf_p <- vcf[as.numeric(vcf$POS) <= split_point]
  vcf_q <- vcf[as.numeric(vcf$POS) > split_point]

  if (nrow(vcf_p) > 0) writevcf_beagle(vcf_p, pathP)
  if (nrow(vcf_q) > 0) writevcf_beagle(vcf_q, pathQ)
}

#' Merge Beagle output from p and q arms back into IMPUTE format
#'
#' @param vcfP Path to p-arm Beagle VCF.
#' @param vcfQ Path to q-arm Beagle VCF.
#' @param outfile Path to the output IMPUTE format file.
#' @export
writebeagle_as_impute_arms <- function(vcfP = NULL, vcfQ = NULL, outfile) {
  read_vcf <- function(path) {
    if (!is.null(path) && file.exists(path)) {
      # Beagle VCFs are gzipped by default
      return(data.table::fread(path, skip = "#CHROM", header = TRUE))
    }
    return(NULL)
  }

  outP <- read_vcf(vcfP)
  outQ <- read_vcf(vcfQ)

  if (is.null(outP) && is.null(outQ)) {
    log_failure("Neither p-arm nor q-arm Beagle output found for merging.")
  }

  combined <- data.table::rbindlist(list(outP, outQ), use.names = TRUE)

  # Extract GT (Genotype)
  gt_col <- names(combined)[10]
  gt_data <- combined[[gt_col]]
  haplo <- data.table::tstrsplit(gt_data, "[|/]")

  impute_dt <- data.table::data.table(
    V1 = "---",
    V2 = combined$ID,
    V3 = combined$POS,
    V4 = combined$REF,
    V5 = combined$ALT,
    V6 = haplo[[1]],
    V7 = haplo[[2]]
  )

  data.table::fwrite(impute_dt, file = outfile, sep = " ", col.names = FALSE, quote = FALSE)
}

#' helper for writing Beagle VCFs
#' @export
writevcf_beagle <- function(vcf, filepath, vcfversion = "4.2", genomereference = "GRCh38") {
  header <- paste0(
    "##fileformat=VCFv", vcfversion, "\n",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    "##reference=", genomereference, "\n"
  )
  cat(header, file = filepath)
  data.table::fwrite(vcf, file = filepath, sep = "\t", append = TRUE, col.names = TRUE, quote = FALSE)
}

#' Generate Beagle input directly from allele counts
#' @export
generate_beagle_input_from_counts <- function(chrom, tumour_allele_counts_file, normal_allele_counts_file,
                                              output_file, imputeinfofile, is_male, problem_loci_file = NA,
                                              heterozygous_filter = 0.1) {
  # Load reference info
  impute_info <- parse_imputeinfofile(imputeinfofile, is_male, chrom = chrom)

  # Load reference legend (using vroom for speed)
  # Expected columns: id, position, a0, a1, type
  log_info("Reading legend from {impute_info$impute_legend}")
  known_SNPs <- vroom::vroom(unlist(impute_info$impute_legend), delim = " ", col_types = "ciccc", show_col_types = FALSE)
  data.table::setDT(known_SNPs)

  # Filter problem SNPs
  if (!is.na(problem_loci_file) && problem_loci_file != "NA" && file.exists(problem_loci_file)) {
    problem_snps_raw <- data.table::fread(problem_loci_file, header = TRUE, sep = "\t", data.table = FALSE)
    problem_positions <- problem_snps_raw$Pos[problem_snps_raw$Chr == chrom]
    known_SNPs <- known_SNPs[!(known_SNPs$position %in% problem_positions), ]
  }

  # Load allele counts using the package's robust reader (handles headers and #)
  log_info("Reading normal allele counts from {normal_allele_counts_file}")
  snp_normal <- read_alleleFrequencies(normal_allele_counts_file)

  # Intersection based on position
  common_pos <- intersect(known_SNPs$position, snp_normal$POS)

  if (length(common_pos) == 0) {
    # Try with chr prefix if match failed
    if (any(grepl("^chr", snp_normal$CHR))) {
      # This is already handled by read_alleleFrequencies returning numeric or whatever
      # But POS is what matters.
    }
    log_failure("No overlap between reference legend and normal allele counts for chr {chrom}. Check positions and chromosome versions.")
  }

  log_info("Found {length(common_pos)} SNPs overlapping with reference for chr {chrom}")

  # Subset and sort both
  valid_known_snps <- known_SNPs[match(common_pos, known_SNPs$position), ]
  found_normal_data <- snp_normal[match(common_pos, snp_normal$POS), ]

  # Define base columns (A=3, C=4, G=5, T=6 in our table)
  bases <- c("A", "C", "G", "T")

  # Get indices for Ref and Alt (a0 and a1)
  # We use the matched normal counts to determine GT
  ref_base_idx <- match(valid_known_snps$a0, bases)
  alt_base_idx <- match(valid_known_snps$a1, bases)

  # Safely extract counts using matrix indexing for speed
  # Columns 3,4,5,6 correspond to bases
  normal_counts_matrix <- as.matrix(found_normal_data[, 3:6, with = FALSE])

  ref_counts <- as.numeric(vapply(seq_along(ref_base_idx), function(i) {
    if (is.na(ref_base_idx[i])) {
      return(0)
    }
    normal_counts_matrix[i, ref_base_idx[i]]
  }, numeric(1)))

  alt_counts <- as.numeric(vapply(seq_along(alt_base_idx), function(i) {
    if (is.na(alt_base_idx[i])) {
      return(0)
    }
    normal_counts_matrix[i, alt_base_idx[i]]
  }, numeric(1)))

  # Combined depth at the reference alleles
  total_counts <- ref_counts + alt_counts
  keep_mask <- total_counts > 0

  if (sum(keep_mask) == 0) {
    log_failure("No SNPs with coverage in normal for chr {chrom}")
  }

  log_info("Keeping {sum(keep_mask)} SNPs with coverage in normal")

  # Subset one last time
  valid_known_snps <- valid_known_snps[keep_mask]
  ref_counts <- ref_counts[keep_mask]
  alt_counts <- alt_counts[keep_mask]
  total_counts <- total_counts[keep_mask]

  bafs <- alt_counts / total_counts

  # Determine genotypes (GT)
  gt <- rep("0/1", length(bafs))
  gt[bafs <= heterozygous_filter] <- "0/0"
  gt[bafs >= (1.0 - heterozygous_filter)] <- "1/1"

  # Format directly for VCF
  chr_vcf <- if (chrom == "23") "X" else as.character(chrom)

  vcf <- data.table::data.table(
    "#CHROM" = rep(chr_vcf, length(gt)),
    POS = valid_known_snps$position,
    ID = valid_known_snps$id,
    REF = valid_known_snps$a0,
    ALT = valid_known_snps$a1,
    QUAL = ".",
    FILTER = "PASS",
    INFO = ".",
    FORMAT = "GT",
    SAMP001 = gt
  )

  log_info("Writing {nrow(vcf)} SNPs to {output_file}")
  writevcf_beagle(vcf, output_file)
}

#' Construct haplotypes for a chromosome
#'
#' This function takes preprocessed data and performs haplotype reconstruction.
#'
#' @param chrom The chromosome for which to reconstruct haplotypes
#' @param tumourname Identifier of the tumour, used to match data files on disk
#' @param normalname Identifier of the normal, used to match data files on disk
#' @param ismale Boolean, set to TRUE if the sample is male
#' @param imputeinfofile Full path to the imputeinfo reference file
#' @param problemloci Full path to the problematic loci reference file
#' @param impute_exe Path to the impute executable (can be found if its in $PATH)
#' @param min_normal_depth Minimal depth in the matched normal required for a SNP to be used
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param snp6_reference_info_file SNP6 only parameter Default: NA
#' @param heterozygous_filter SNP6 only parameter Default: NA
#' @param usebeagle Should use beagle5 instead of impute2 Default: FALSE
#' @param beaglejar Full path to Beagle java jar file Default: NA
#' @param beagleref Full path to Beagle reference file Default: NA
#' @param beagleplink Full path to Beagle plink file  Default: NA
#' @param beaglemaxmem Integer Beagle max heap size in Gb  Default: 10
#' @param beaglenthreads Integer number of threads used by beagle5 Default:1
#' @param beaglewindow Integer size of the genomic window for beagle5 (cM) Default:40
#' @param beagleoverlap Integer size of the overlap between windows beagle5 Default:4
#' @param javajre Path to the Java JRE executable (default java, i.e. in $PATH)
#' @author sd11, maxime.tarabichi, jdemeul
#' @export
convert_beagle_to_impute <- function(beagle_file, output_file) {
  # Robust VCF reading: Beagle files are often gzipped and might have sparse headers
  if (!file.exists(beagle_file)) {
    log_failure("Beagle VCF file not found: {beagle_file}")
  }

  # First try reading with skip="#CHROM"
  vcf <- tryCatch(
    {
      data.table::fread(beagle_file, skip = "#CHROM", header = TRUE)
    },
    error = function(e) {
      # Fallback: if #CHROM is missing, try reading without skip if the file is tiny
      if (file.info(beagle_file)$size < 500) {
        return(data.table::data.table())
      }
      stop(e)
    }
  )

  # If we have no data, return empty table
  if (nrow(vcf) == 0) {
    log_info("Beagle VCF is empty. Writing empty output.")
    data.table::fwrite(data.table::data.table(), file = output_file, sep = " ", col.names = FALSE)
    return(NULL)
  }

  # Identify the GT data column (standard VCF col 10)
  if (ncol(vcf) < 10) {
    log_warning("Beagle VCF file {beagle_file} has fewer than 10 columns. Writing empty output.")
    data.table::fwrite(data.table::data.table(), file = output_file, sep = " ", col.names = FALSE)
    return(NULL)
  }

  gt_data <- vcf[[10]]
  gt_only <- data.table::tstrsplit(gt_data, ":")[[1]]
  haplo <- data.table::tstrsplit(gt_only, "[|/]")

  if (length(haplo) < 2) {
    log_warning("Could not parse genotypes from Beagle VCF {beagle_file}. Writing empty output.")
    data.table::fwrite(data.table::data.table(), file = output_file, sep = " ", col.names = FALSE)
    return(NULL)
  }

  # Construct IMPUTE2 format
  impute_dt <- data.table::data.table(
    V1 = "---",
    V2 = vcf[["ID"]],
    V3 = as.integer(as.numeric(vcf[["POS"]])),
    V4 = vcf[["REF"]],
    V5 = vcf[["ALT"]],
    V6 = haplo[[1]],
    V7 = haplo[[2]]
  )

  log_info("Extracted {nrow(impute_dt)} SNPs from Beagle VCF")

  # Write out space-separated, no header (as expected by GetChromosomeBAFs read logic 'header=FALSE')
  data.table::fwrite(impute_dt, file = output_file, sep = " ", col.names = FALSE, quote = FALSE)
}

#' @param impute_results_dir Directory containing the impute/beagle output files
#' @author sd11, maxime.tarabichi, jdemeul
#' @export
run_haplotyping <- function(
  chrom, tumourname, normalname,
  ismale, problemloci,
  impute_results_dir, min_normal_depth, chrom_names,
  imputeinfofile = NA,
  externalhaplotypeprefix = NA,
  use_previous_imputation = FALSE,
  snp6_reference_info_file = NA,
  heterozygous_filter = NA,
  beagle_input_dir = NA,
  allele_frequencies_dir = NA,
  chrom_coord_file = NA
) {
  # determine if we are using beagle based on beagle_input_dir
  usebeagle <- !is.na(beagle_input_dir)

  if (usebeagle) {
    # Check if we already have the IMPUTE-converted file locally first (prevents redundant conversion)
    local_haplo <- paste0(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt")
    if (file.exists(local_haplo)) {
      haplotype_file <- local_haplo
    } else {
      # Construct path to Beagle VCF
      # If beagle_input_dir is provided, look there.
      beagle_search_dir <- if (!is.na(beagle_input_dir)) beagle_input_dir else impute_results_dir
      haplotype_file <- local_haplo

      beagle_vcf_p <- file.path(beagle_search_dir, paste0(tumourname, "_beagle5_output_chr", chrom, "_P.vcf.gz"))
    beagle_vcf_q <- file.path(beagle_search_dir, paste0(tumourname, "_beagle5_output_chr", chrom, "_Q.vcf.gz"))

    if (file.exists(beagle_vcf_p) || file.exists(beagle_vcf_q)) {
      log_info("Merging Beagle arm-specific outputs for chr {chrom}")
      writebeagle_as_impute_arms(
        vcfP = if (file.exists(beagle_vcf_p)) beagle_vcf_p else NULL,
        vcfQ = if (file.exists(beagle_vcf_q)) beagle_vcf_q else NULL,
        outfile = haplotype_file
      )
    } else {
      # Fallback to single file patterns
      beagle_patterns <- c(
        paste0(tumourname, "_beagle5_output_chr", chrom, ".txt.vcf.gz"),
        paste0(tumourname, "_beagle5_output_chr", chrom, ".txt.vcf"),
        paste0(tumourname, "_beagle_output_chr", chrom, ".vcf.gz"),
        paste0(tumourname, "_beagle_output_chr", chrom, ".vcf")
      )

      beagle_vcf <- NA
      for (pat in beagle_patterns) {
        temp_path <- file.path(beagle_search_dir, pat)
        if (file.exists(temp_path)) {
          beagle_vcf <- temp_path
          break
        }
      }

      if (is.na(beagle_vcf)) {
        log_failure("Expected Beagle VCF file not found in {beagle_search_dir} (single file or arm-specific).")
      }

      log_info("Converting Beagle VCF to IMPUTE format: {beagle_vcf} -> {haplotype_file}")
      convert_beagle_to_impute(beagle_vcf, haplotype_file)
    }
  } else {
    # Non-Beagle (Standard Impute2) mode
    # Local first check
    local_haplo <- paste0(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt")
    if (file.exists(local_haplo)) {
      haplotype_file <- local_haplo
    } else if (!is.na(impute_results_dir)) {
      haplotype_file <- file.path(impute_results_dir, local_haplo)
      if (!file.exists(haplotype_file)) {
        log_failure("Expected haplotype file not found: {haplotype_file}")
      }
    } else {
      log_failure("No haplotype file found and no impute_results_dir provided.")
    }
  }


  # If an allele counts file exists we assume this is a WGS sample and run the corresponding step, otherwise it must be SNP6
  if (is.na(allele_frequencies_dir)) {
    log_failure("allele_frequencies_dir must be provided to run_haplotyping")
  }
  # Use robust find_file logic for allele frequencies
  find_ac_file <- function(dir, sample, chrom) {
    p1 <- file.path(dir, paste0(sample, "_alleleFrequencies_chr", chrom, ".txt"))
    if (file.exists(p1)) return(p1)
    norm_c <- gsub("chr", "", as.character(chrom), ignore.case = TRUE)
    p2 <- file.path(dir, paste0(sample, "_alleleFrequencies_chr", norm_c, ".txt"))
    if (file.exists(p2)) return(p2)
    p3 <- file.path(dir, paste0(sample, "_alleleFrequencies_", norm_c, ".txt"))
    if (file.exists(p3)) return(p3)
    return(NULL)
  }

  allelefrequenciesfile <- find_ac_file(allele_frequencies_dir, tumourname, chrom)

  if (file.exists(allelefrequenciesfile)) {
    # WGS - Transform the impute output into haplotyped BAFs

    # if present, input external haplotype blocks
    if (!is.na(externalhaplotypeprefix) && file.exists(paste0(externalhaplotypeprefix, chrom, ".vcf"))) {
      log_info("Adding in the external haplotype blocks")

      # output BAFs to plot pre-external haplotyping
      GetChromosomeBAFs(
        chrom = chrom,
        SNP_file = allelefrequenciesfile,
        haplotypeFile = haplotype_file,
        samplename = tumourname,
        outfile = paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt", sep = ""),
        chr_names = chrom_names,
        minCounts = min_normal_depth
      )

      # Plot what we have before external haplotyping is incorporated
      plot_haplotype_data(
        haplotyped_baf_file = paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt", sep = ""),
        image_file_name = paste(tumourname, "_chr", chrom, "_heterozygousData_noExt.png", sep = ""),
        samplename = tumourname,
        chrom = chrom
      )

      input_known_haplotypes(
        chrom = chrom,
        chrom_names = chrom_names,
        imputedHaplotypeFile = haplotype_file,
        externalHaplotypeFile = paste0(externalhaplotypeprefix, chrom, ".vcf")
      )
    }

    GetChromosomeBAFs(
      chrom = chrom,
      SNP_file = allelefrequenciesfile,
      haplotypeFile = haplotype_file,
      samplename = tumourname,
      outfile = paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
      chr_names = chrom_names,
      minCounts = min_normal_depth
    )
  } else {
    log_info("SNP6 get BAFs")
    # SNP6 - Transform the impute output into haplotyped BAFs
    GetChromosomeBAFs_SNP6(
      chrom = chrom,
      alleleFreqFile = paste(tumourname, "_impute_input_chr", chrom, "_withAlleleFreq.csv", sep = ""),
      haplotypeFile = haplotype_file,
      samplename = tumourname,
      outputfile = paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
      chr_names = chrom_names
    )
  }

  # Plot what we have until this point
  plot_haplotype_data(
    haplotyped_baf_file = paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
    image_file_name = paste(tumourname, "_chr", chrom, "_heterozygousData.png", sep = ""),
    samplename = tumourname,
    chrom = chrom
  )
}

#' Construct haplotypes for a chromosome - germline WGS version
#'
#' This function takes preprocessed data and performs haplotype reconstruction.
#'
#' @param chrom The chromosome for which to reconstruct haplotypes
#' @param germlinename Identifier of the germline sample, used to match data files on disk
#' @param normalname Identifier of the reconstructed normal, used to match data files on disk
#' @param ismale Boolean, set to TRUE if the sample is male
#' @param imputeinfofile Full path to the imputeinfo reference file
#' @param problemloci Full path to the problematic loci reference file
#' @param impute_exe Path to the impute executable (can be found if its in $PATH)
#' @param min_normal_depth Minimal depth in the matched normal required for a SNP to be used
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param snp6_reference_info_file SNP6 only parameter Default: NA
#' @param heterozygous_filter SNP6 only parameter Default: NA
#' @param usebeagle Should use beagle5 instead of impute2 Default: FALSE
#' @param beaglejar Full path to Beagle java jar file Default: NA
#' @param beagleref Full path to Beagle reference file Default: NA
#' @param beagleplink Full path to Beagle plink file  Default: NA
#' @param beaglemaxmem Integer Beagle max heap size in Gb  Default: 10
#' @param beaglenthreads Integer number of threads used by beagle5 Default:1
#' @param beaglewindow Integer size of the genomic window for beagle5 (cM) Default:40
#' @param beagleoverlap Integer size of the overlap between windows beagle5 Default:4
#' @param javajre Path to the Java JRE executable (default java, i.e. in $PATH)
#' @author sd11, maxime.tarabichi, jdemeul, Naser Ansari-Pour (BDI, Oxford)
#' @export

#' @param usebeagle Logical, if TRUE expects Beagle VCF output and converts to IMPUTE format.
#' @author sd11, maxime.tarabichi, jdemeul, Naser Ansari-Pour (BDI, Oxford)
#' @export
run_haplotyping_germline <- function(
  chrom, germlinename, normalname, ismale, problemloci,
  impute_results_dir, min_normal_depth, chrom_names,
  imputeinfofile = NA,
  externalhaplotypeprefix = NA,
  use_previous_imputation = FALSE,
  snp6_reference_info_file = NA, heterozygous_filter = NA,
  beagle_input_dir = NA,
  allele_frequencies_dir = NA,
  chrom_coord_file = NA
) {
  # determine if we are using beagle based on beagle_input_dir
  usebeagle <- !is.na(beagle_input_dir)

  # Point to the existing haplotype file in the external directory
  if (usebeagle) {
    # If beagle_input_dir is provided, look there.
    beagle_search_dir <- if (!is.na(beagle_input_dir)) beagle_input_dir else impute_results_dir

    # Try multiple common naming patterns for Beagle VCFs
    beagle_patterns <- c(
      paste0(germlinename, "_beagle5_output_chr", chrom, ".txt.vcf.gz"),
      paste0(germlinename, "_beagle5_output_chr", chrom, ".txt.vcf"),
      paste0(germlinename, "_beagle_output_chr", chrom, ".vcf.gz"),
      paste0(germlinename, "_beagle_output_chr", chrom, ".vcf")
    )

    beagle_vcf <- NA
    for (pat in beagle_patterns) {
      temp_path <- file.path(beagle_search_dir, pat)
      if (file.exists(temp_path)) {
        beagle_vcf <- temp_path
        break
      }
    }

    if (is.na(beagle_vcf)) {
      log_failure("Expected Beagle VCF file not found in {beagle_search_dir}. Tried patterns: {paste(beagle_patterns, collapse=', ')}")
    }

    haplotype_file <- paste0(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt")
    log_info("Converting Beagle VCF to IMPUTE format: {beagle_vcf} -> {haplotype_file}")
    convert_beagle_to_impute(beagle_vcf, haplotype_file)
  } else {
    haplotype_file <- file.path(impute_results_dir, paste0(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"))
    if (!file.exists(haplotype_file)) {
      log_failure("Expected haplotype file not found: {haplotype_file}")
    }
  }

  if (is.na(allele_frequencies_dir)) {
    log_failure("allele_frequencies_dir must be provided to run_haplotyping_germline")
  }
  allelefrequenciesfile <- file.path(allele_frequencies_dir, paste0(germlinename, "_alleleFrequencies_chr", chrom, ".txt"))

  if (file.exists(allelefrequenciesfile)) {
    # WGS - Transform the impute output into haplotyped BAFs

    # if present, input external haplotype blocks
    if (!is.na(externalhaplotypeprefix) && file.exists(paste0(externalhaplotypeprefix, chrom, ".vcf"))) {
      log_info("Adding in the external haplotype blocks")

      # output BAFs to plot pre-external haplotyping
      GetChromosomeBAFs(
        chrom = chrom,
        SNP_file = allelefrequenciesfile,
        haplotypeFile = haplotype_file,
        samplename = germlinename,
        outfile = paste(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt", sep = ""),
        chr_names = chrom_names,
        minCounts = min_normal_depth
      )

      # Plot what we have before external haplotyping is incorporated
      plot_haplotype_data(
        haplotyped_baf_file = paste(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt", sep = ""),
        image_file_name = paste(germlinename, "_chr", chrom, "_heterozygousData_noExt.png", sep = ""),
        samplename = germlinename,
        chrom = chrom
      )

      input_known_haplotypes(
        chrom = chrom,
        chrom_names = chrom_names,
        imputedHaplotypeFile = haplotype_file,
        externalHaplotypeFile = paste0(externalhaplotypeprefix, chrom, ".vcf")
      )
    }

    GetChromosomeBAFs(
      chrom = chrom,
      SNP_file = allelefrequenciesfile,
      haplotypeFile = haplotype_file,
      samplename = germlinename,
      outfile = paste(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
      chr_names = chrom_names,
      minCounts = min_normal_depth
    )
  } else {
    log_failure("Germline calling is only on WGS data - SNParray data not sufficiently dense")
  }

  # Plot what we have until this point
  plot_haplotype_data(
    haplotyped_baf_file = paste(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
    image_file_name = paste(germlinename, "_chr", chrom, "_heterozygousData.png", sep = ""),
    samplename = germlinename,
    chrom = chrom
  )
}
