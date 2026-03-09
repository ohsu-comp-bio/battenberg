#' Obtain BAF and LogR from the allele counts (Memory Optimized)
#' @export
getBAFsAndLogRs <- function(tumourAlleleCountsFile.prefix, normalAlleleCountsFile.prefix, figuresFile.prefix, BAFnormalFile, BAFmutantFile, logRnormalFile, logRmutantFile, combinedAlleleCountsFile, chr_names, g1000file.prefix, minCounts = NA, samplename = "sample1", seed = as.integer(Sys.time())) {
  set.seed(seed)

  # Initialize files (delete if already exists to avoid double-appending)
  out_files <- c(BAFnormalFile, BAFmutantFile, logRnormalFile, logRmutantFile, combinedAlleleCountsFile)
  for (f in out_files) if (file.exists(f)) file.remove(f)

  # Containers for thinned plotting data (to prevent graphical OOM)
  plot_data_list <- list()
  total_snps_processed <- 0

  for (chrom in chr_names) {
    log_info("Processing chromosome {chrom}...")

    # Load data for THIS chromosome only
    input_data <- concatenateAlleleCountFiles(tumourAlleleCountsFile.prefix, ".txt", chrom)
    normal_input_data <- concatenateAlleleCountFiles(normalAlleleCountsFile.prefix, ".txt", chrom)
    allele_data <- concatenateG1000SnpFiles(g1000file.prefix, ".txt", chrom)

    log_info("  - Raw SNPs: Tumour={nrow(input_data)}, Normal={nrow(normal_input_data)}, G1000={nrow(allele_data)}")

    if (nrow(input_data) == 0 || nrow(normal_input_data) == 0 || nrow(allele_data) == 0) {
      log_warning("  - Missing data for chromosome {chrom}. Skipping.")
      next
    }

    # Convert to data.table
    data.table::setDT(input_data)
    data.table::setDT(normal_input_data)
    data.table::setDT(allele_data)

    # Standardize
    input_data[[1]] <- gsub("chr", "", as.character(input_data[[1]]))
    normal_input_data[[1]] <- gsub("chr", "", as.character(normal_input_data[[1]]))
    allele_data[[1]] <- gsub("chr", "", as.character(allele_data[[1]]))

    names(allele_data)[1:4] <- c("CHR", "POS", "A0", "A1")
    names(normal_input_data)[1:7] <- c("CHR", "POS", "nCountA", "nCountC", "nCountG", "nCountT", "nDepth")
    names(input_data)[1:7] <- c("CHR", "POS", "tCountA", "tCountC", "tCountG", "tCountT", "tDepth")

    # Ensure types match for join
    input_data[, `:=`(CHR = as.character(CHR), POS = as.integer(POS))]
    normal_input_data[, `:=`(CHR = as.character(CHR), POS = as.integer(POS))]
    allele_data[, `:=`(CHR = as.character(CHR), POS = as.integer(POS))]

    # Fast Join logic
    data.table::setkey(input_data, CHR, POS)
    data.table::setkey(normal_input_data, CHR, POS)
    data.table::setkey(allele_data, CHR, POS)

    # Join
    joined <- normal_input_data[input_data, nomatch = 0]
    joined <- allele_data[joined, nomatch = 0]

    log_info("  - Synced SNPs: {nrow(joined)}")

    if (nrow(joined) == 0) {
      log_warning("  - Zero overlap for chromosome {chrom}. Check reference compatibility.")
      next
    }

    # cleanup temp objects
    rm(input_data, normal_input_data, allele_data)

    # Matrix extraction
    norm_m <- as.matrix(joined[, .(nCountA, nCountC, nCountG, nCountT)])
    mut_m <- as.matrix(joined[, .(tCountA, tCountC, tCountG, tCountT)])

    len <- nrow(joined)
    idx_matrix <- cbind(seq_len(len), as.integer(joined$A0))
    idx_matrix2 <- cbind(seq_len(len), as.integer(joined$A1))

    normCount1 <- norm_m[idx_matrix]
    normCount2 <- norm_m[idx_matrix2]
    mutCount1 <- mut_m[idx_matrix]
    mutCount2 <- mut_m[idx_matrix2]

    totalNormal <- normCount1 + normCount2
    totalMutant <- mutCount1 + mutCount2

    rm(norm_m, mut_m)

    # Apply coverage filters
    valid_indices <- seq_len(len)
    if (!is.na(minCounts)) {
      valid_indices <- which(totalNormal >= minCounts & totalMutant >= 1)
      totalNormal <- totalNormal[valid_indices]
      totalMutant <- totalMutant[valid_indices]
      normCount1 <- normCount1[valid_indices]
      normCount2 <- normCount2[valid_indices]
      mutCount1 <- mutCount1[valid_indices]
      mutCount2 <- mutCount2[valid_indices]
    }

    n <- length(valid_indices)
    log_info("  - Final Filtered SNPs: {n}")

    if (n == 0) {
      log_warning("  - No SNPs passed coverage filters for {chrom}.")
      next
    }

    # BAF/LogR Calc
    selector <- round(stats::runif(n))
    is_zero <- selector == 0
    is_one <- !is_zero

    normalBAF <- numeric(n)
    mutantBAF <- numeric(n)
    normalBAF[is_zero] <- normCount1[is_zero] / totalNormal[is_zero]
    normalBAF[is_one] <- normCount2[is_one] / totalNormal[is_one]
    mutantBAF[is_zero] <- mutCount1[is_zero] / totalMutant[is_zero]
    mutantBAF[is_one] <- mutCount2[is_one] / totalMutant[is_one]

    mutantLogR_raw <- totalMutant / totalNormal
    # Mean shift will be approximate per chromosome here, but we can fix the global mean shift later
    # Actually, original code used log2(ratio / mean(all_ratios))
    # For now, let's keep the raw ratio and we'll normalize at the very end of this loop?
    # No, let's calculate the log2(ratio) and keep the global mean shift in mind.
    # Actually, we should probably calculate the global mean first...
    # But that requires loading all ratios.
    # Let's just use log2(ratio) and we'll shift the file afterwards.
    tumorLogR_unshifted <- log2(mutantLogR_raw)

    CHR_final <- joined$CHR[valid_indices]
    POS_final <- joined$POS[valid_indices]

    # Write results appending to disk
    baseDT <- data.table::data.table(Chromosome = CHR_final, Position = POS_final)

    # Normal BAF
    baseDT[[samplename]] <- normalBAF
    data.table::fwrite(baseDT, file = BAFnormalFile, sep = "\t", append = TRUE, col.names = !file.exists(BAFnormalFile))

    # Mutant BAF
    baseDT[[samplename]] <- mutantBAF
    data.table::fwrite(baseDT, file = BAFmutantFile, sep = "\t", append = TRUE, col.names = !file.exists(BAFmutantFile))

    # Normal LogR
    baseDT[[samplename]] <- integer(n)
    data.table::fwrite(baseDT, file = logRnormalFile, sep = "\t", append = TRUE, col.names = !file.exists(logRnormalFile))

    # Mutant LogR
    baseDT[[samplename]] <- tumorLogR_unshifted
    data.table::fwrite(baseDT, file = logRmutantFile, sep = "\t", append = TRUE, col.names = !file.exists(logRmutantFile))

    # Combined counts
    baseDT[[samplename]] <- NULL
    combinedDT <- cbind(baseDT, data.table::data.table(
      mutCountT1 = mutCount1, mutCountT2 = mutCount2,
      mutCountN1 = normCount1, mutCountN2 = normCount2
    ))
    data.table::fwrite(combinedDT, file = combinedAlleleCountsFile, sep = "\t", append = TRUE, col.names = !file.exists(combinedAlleleCountsFile))

    # Thinned plotting data: keep 1 in every 25 SNPs
    thin_idx <- seq(1, n, by = 25)
    plot_data_list[[chrom]] <- data.table::data.table(
      Chromosome = CHR_final[thin_idx],
      Position = POS_final[thin_idx],
      Tumor_LogR = tumorLogR_unshifted[thin_idx],
      Tumor_BAF = mutantBAF[thin_idx],
      Germline_BAF = normalBAF[thin_idx]
    )

    total_snps_processed <- total_snps_processed + n
    rm(joined, baseDT, combinedDT, normalBAF, mutantBAF, tumorLogR_unshifted)
    gc()
  }

  log_info("Sync complete. Total SNPs processed across all chromosomes: {total_snps_processed}")

  # GLOBAL MEAN SHIFT for LogR (Battenberg requires center at 0)
  log_info("Performing global LogR mean shift...")
  # We read the LogR column to calculate the global mean.
  # vroom is faster for column selection on large files.
  global_mean <- mean(vroom::vroom(logRmutantFile, col_select = 3, show_col_types = FALSE)[[1]], na.rm = TRUE)
  log_info("Global LogR Mean: {global_mean}. Shifting values...")

  # Read full file, shift, write. (This is high RAM but only for 2 columns Chrom/Pos + 1 Float)
  # 28M rows * 3 cols * 8 bytes ≈ 672 MB. Totally safe.
  full_logr <- data.table::fread(logRmutantFile)
  full_logr[[3]] <- full_logr[[3]] - global_mean
  data.table::fwrite(full_logr, file = logRmutantFile, sep = "\t")
  rm(full_logr)
  gc()

  # CONSTRUCT PLOTTING OBJECT (FROM THINNED DATA)
  log_info("Constructing thinned ASCAT plot...")
  plot_data <- data.table::rbindlist(plot_data_list)
  # Standardize Chromosome names for ASCAT factor sorting
  ch <- lapply(chr_names, function(x) {
    # Match robustly (handling both '1' and 'chr1' in the data)
    normalized_data_chrs <- gsub("chr", "", as.character(plot_data$Chromosome))
    normalized_target_chr <- gsub("chr", "", as.character(x))
    tmp <- which(normalized_data_chrs == normalized_target_chr)

    if (length(tmp) == 0) {
      return(numeric(0))
    }
    return(tmp[1]:tmp[length(tmp)])
  })

  ascat_bc <- list(
    Tumor_LogR = data.frame(plot_data$Tumor_LogR - global_mean),
    Tumor_BAF = data.frame(plot_data$Tumor_BAF),
    Germline_LogR = data.frame(integer(nrow(plot_data))),
    Germline_BAF = data.frame(plot_data$Germline_BAF),
    Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL,
    Tumor_counts = NULL, Germline_counts = NULL,
    SNPpos = data.frame(Chromosome = plot_data$Chromosome, Position = plot_data$Position, stringsAsFactors = FALSE),
    chrs = chr_names,
    samples = samplename,
    chrom = split_genome(plot_data[, 1:2]),
    ch = ch
  )
  ASCAT::ascat.plotRawData(ascat_bc)
}

#' Prepare data for impute
#'
#' @param chrom The chromosome for which impute input should be generated.
#' @param tumour_allele_counts_file Output from the allele counter on the matched tumour for this chromosome.
#' @param normal_allele_counts_file Output from the allele counter on the matched normal for this chromosome.
#' @param output_file File where the impute input for this chromosome will be written.
#' @param imputeinfofile Info file with impute reference information.
#' @param is_male Boolean denoting whether this sample is male (TRUE), or female (FALSE).
#' @param problem_loci_file A file containing genomic locations that must be discarded (optional).
#' @param use_loci_file A file containing genomic locations that must be included (optional).
#' @param heterozygous_filter The cutoff where a SNP will be considered as heterozygous (default 0.1).
#' @author dw9, sd11
#' @export
generate_impute_input_wgs <- function(
  chrom, tumour_allele_counts_file, normal_allele_counts_file,
  output_file, imputeinfofile, is_male, problem_loci_file = NA,
  use_loci_file = NA, heterozygous_filter = 0.1
) {
  # Read in the reference file paths for the specified chrom
  impute_info <- parse_imputeinfofile(imputeinfofile, is_male, chrom = chrom)
  chrom_name <- chrom

  # Efficiently load and combine known SNP legend files
  # Replaces the for-loop/rbind pattern which is very slow in R
  # Efficiently load known SNP legend files using vroom
  known_SNPs <- vroom::vroom(
    unlist(impute_info$impute_legend),
    delim = " ",
    show_col_types = FALSE
  )
  data.table::setDF(known_SNPs)

  # Filter out 'problem' SNPs (BAF streaks)
  if (!is.na(problem_loci_file) && problem_loci_file != "NA") {
    problem_snps_raw <- data.table::fread(problem_loci_file, header = TRUE, sep = "auto", data.table = FALSE)
    problem_positions <- problem_snps_raw$Pos[problem_snps_raw$Chr == chrom_name]
    known_SNPs <- known_SNPs[!(known_SNPs$position %in% problem_positions), ]
  }

  # Filter for 'good' SNPs (e.g., SNP6 positions)
  if (!is.na(use_loci_file) && use_loci_file != "NA") {
    good_snps_raw <- data.table::fread(use_loci_file, header = TRUE, sep = "auto", data.table = FALSE)
    good_positions <- good_snps_raw$pos[good_snps_raw$chr == chrom_name]
    known_SNPs <- known_SNPs[known_SNPs$position %in% good_positions, ]
  }

  # Load allele counts using fread (ignoring comments)
  # Tumour and Normal are combined column-wise to match legacy indexing
  snp_tumour <- data.table::fread(tumour_allele_counts_file, sep = "auto", header = FALSE, data.table = FALSE)
  snp_normal <- data.table::fread(normal_allele_counts_file, sep = "auto", header = FALSE, data.table = FALSE)

  # Combined data: [Tumour Cols 1-6] [Normal Cols 7-12]
  snp_combined <- cbind(snp_tumour, snp_normal)

  # Match known SNPs to the allele counter positions
  indices <- match(known_SNPs$position, snp_combined[, 2])
  mask <- !is.na(indices)
  found_snp_data <- snp_combined[indices[mask], ]
  valid_known_snps <- known_SNPs[mask, ]

  # Calculate BAF for the NORMAL sample to determine genotypes
  # Logic: Alt / (Alt + Ref).
  # Ref column index: match allele in col 3 + normal offset (ncol) + 2
  # Alt column index: match allele in col 4 + normal offset (ncol) + 2
  nucleotides <- c("A", "C", "G", "T")
  norm_col_count <- ncol(snp_normal)

  ref_cols <- match(valid_known_snps[, 3], nucleotides) + norm_col_count + 2
  alt_cols <- match(valid_known_snps[, 4], nucleotides) + norm_col_count + 2

  # Matrix indexing for high-speed extraction of specific allele counts
  row_idx <- seq_len(nrow(found_snp_data))
  alt_counts <- as.numeric(found_snp_data[cbind(row_idx, alt_cols)])
  ref_counts <- as.numeric(found_snp_data[cbind(row_idx, ref_cols)])

  bafs <- alt_counts / (alt_counts + ref_counts)
  bafs[is.nan(bafs)] <- 0

  # Determine genotypes for IMPUTE2 (1-hot encoded: HomRef, Het, HomAlt)
  min_baf <- min(heterozygous_filter, 1.0 - heterozygous_filter)
  max_baf <- max(heterozygous_filter, 1.0 - heterozygous_filter)

  genotypes <- matrix(0, nrow = nrow(found_snp_data), ncol = 3)
  genotypes[bafs <= min_baf, 1] <- 1
  genotypes[bafs > min_baf & bafs < max_baf, 2] <- 1
  genotypes[bafs >= max_baf, 3] <- 1

  # Create final output table
  # Format: [snpID] [Chr] [Pos] [Ref] [Alt] [G1] [G2] [G3]
  snp_names <- paste0("snp", seq_len(nrow(genotypes)))
  out_data <- cbind(snp_names, valid_known_snps[, 1:4], genotypes)

  # Write main output
  data.table::fwrite(out_data, file = output_file, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Legacy check: Write sample_g.txt if chrom_name is NA (usually for non-standard chrom processing)
  if (is.na(chrom_name)) {
    sample_g_file <- file.path(dirname(output_file), "sample_g.txt")
    sample_g_data <- data.frame(
      ID_1 = c(0, "INDIVI1"),
      ID_2 = c(0, "INDIVI1"),
      missing = c(0, 0),
      sex = c("D", 2)
    )
    data.table::fwrite(sample_g_data, file = sample_g_file, sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}

#' Function to correct LogR for waivyness that correlates with GC content
#' @param Tumour_LogR_file String pointing to the tumour LogR output
#' @param outfile String pointing to where the GC corrected LogR should be written
#' @param correlations_outfile File where correlations are to be saved
#' @param gc_content_file_prefix String pointing to where GC windows for this reference genome can be
#' found. These files should be split per chromosome and this prefix must contain the full path until
#' chr in its name. The .txt extension is automatically added.
#' @param replic_timing_file_prefix Like the gc_content_file_prefix, containing replication timing info (supply NULL if no replication timing correction is to be applied)
#' @param chrom_names A vector containing chromosome names to be considered
#' @param recalc_corr_afterwards Set to TRUE to recalculate correlations after correction
#' @author jdemeul, sd11
#' @export
gc_correct_wgs <- function(Tumour_LogR_file, outfile, correlations_outfile, gc_content_file_prefix, replic_timing_file_prefix, chrom_names) {
  if (is.null(gc_content_file_prefix)) log_failure("GC content reference files must be supplied")

  log_info("Starting two-pass memory-optimized GC correction...")

  # Helper to identify reference file properties (names, index presence)
  get_ref_info <- function(f) {
    if (!file.exists(f)) {
      return(NULL)
    }
    # Use suppressWarnings ONLY once to peek at the format
    h_orig <- suppressWarnings(names(data.table::fread(f, nrows = 0)))
    d_check <- suppressWarnings(data.table::fread(f, nrows = 5, header = FALSE))
    has_idx <- ncol(d_check) > length(h_orig)

    h_clean <- h_orig
    if ("chr" %in% h_clean) h_clean[h_clean == "chr"] <- "Chromosome"
    if ("pos" %in% h_clean) h_clean[h_clean == "pos"] <- "Position"
    wins <- setdiff(h_clean, c("Chromosome", "Position"))

    return(list(has_index = has_idx, orig_names = h_orig, clean_names = h_clean, win_cols = wins))
  }

  # Helper to load reference files robustly without causing fread warnings
  load_ref_dt <- function(f, info) {
    if (info$has_index) {
      dt <- data.table::fread(f, skip = 1, header = FALSE, col.names = c("V1_idx", info$clean_names))
      return(dt[, -1, with = FALSE])
    } else {
      # Use col.names even if no index to ensure standardized names (Chromosome/Position)
      dt <- data.table::fread(f, header = TRUE, col.names = info$clean_names)
      return(dt)
    }
  }

  # Peeking at the first GC file
  first_gc_file <- paste0(gc_content_file_prefix, chrom_names[1], ".txt.gz")
  if (!file.exists(first_gc_file)) log_failure("GC reference file not found: {first_gc_file}")
  gc_info <- get_ref_info(first_gc_file)
  win_cols <- gc_info$win_cols
  log_info("GC Reference Windows: {paste(win_cols, collapse=', ')}")

  # Accumulators for cross-genome correlation statistics
  N_vec <- setNames(numeric(length(win_cols)), win_cols)
  SX_vec <- setNames(numeric(length(win_cols)), win_cols)
  SXX_vec <- setNames(numeric(length(win_cols)), win_cols)
  SXY_vec <- setNames(numeric(length(win_cols)), win_cols)
  SY <- 0
  SYY <- 0
  Total_N <- 0

  has_replic <- !is.null(replic_timing_file_prefix) && !is.na(replic_timing_file_prefix)
  rep_info <- NULL
  rep_win_cols <- NULL
  if (has_replic) {
    first_rep_file <- paste0(replic_timing_file_prefix, chrom_names[1], ".txt.gz")
    rep_info <- get_ref_info(first_rep_file)
    if (!is.null(rep_info)) {
      rep_win_cols <- rep_info$win_cols
      RN_vec <- setNames(numeric(length(rep_win_cols)), rep_win_cols)
      RSX_vec <- setNames(numeric(length(rep_win_cols)), rep_win_cols)
      RSXX_vec <- setNames(numeric(length(rep_win_cols)), rep_win_cols)
      RSXY_vec <- setNames(numeric(length(rep_win_cols)), rep_win_cols)
    } else {
      has_replic <- FALSE
    }
  }

  log_info("Pass 1: Identifying best GC windows via online correlation accumulation...")
  all_logr <- data.table::fread(Tumour_LogR_file) # High but manageable RAM usage
  all_logr[, `:=`(Chromosome = gsub("chr", "", as.character(Chromosome)), Position = as.integer(Position))]
  data.table::setkey(all_logr, Chromosome, Position)

  for (cn in chrom_names) {
    log_info("  - Pass 1: Processing {cn}...")
    gc_f <- paste0(gc_content_file_prefix, cn, ".txt.gz")
    if (!file.exists(gc_f)) next
    dt_gc <- load_ref_dt(gc_f, gc_info)
    dt_gc[, `:=`(Chromosome = gsub("chr", "", as.character(Chromosome)), Position = as.integer(Position))]
    sub_logr <- all_logr[gsub("chr", "", as.character(cn))]

    if (nrow(sub_logr) == 0) next

    data.table::setkey(dt_gc, Position)
    data.table::setkey(sub_logr, Position)
    m <- dt_gc[sub_logr, nomatch = 0]
    log_info("    - Joined with GC: {nrow(m)} SNPs")
    if (nrow(m) == 0) next

    y <- as.numeric(m[[ncol(m)]])
    SY <- SY + sum(y, na.rm = TRUE)
    SYY <- SYY + sum(y^2, na.rm = TRUE)
    Total_N <- Total_N + length(y)

    for (w in win_cols) {
      if (!w %in% names(m)) next
      x <- as.numeric(m[[w]])
      valid <- !is.na(x) & !is.na(y)
      N_vec[w] <- N_vec[w] + sum(valid)
      SX_vec[w] <- SX_vec[w] + sum(x[valid])
      SXX_vec[w] <- SXX_vec[w] + sum(x[valid]^2)
      SXY_vec[w] <- SXY_vec[w] + sum(x[valid] * y[valid])
    }

    if (has_replic) {
      rep_f <- paste0(replic_timing_file_prefix, cn, ".txt.gz")
      if (file.exists(rep_f)) {
        dt_rep <- load_ref_dt(rep_f, rep_info)
        dt_rep[, `:=`(Chromosome = gsub("chr", "", as.character(Chromosome)), Position = as.integer(Position))]
        data.table::setkey(dt_rep, Position)
        mr <- dt_rep[m, nomatch = 0]
        log_info("    - Joined with Replication: {nrow(mr)} SNPs")
        if (nrow(mr) > 0) {
          yr <- as.numeric(mr[[ncol(mr)]])
          for (rw in rep_win_cols) {
            if (!rw %in% names(mr)) next
            rx <- as.numeric(mr[[rw]])
            v <- !is.na(rx) & !is.na(yr)
            RN_vec[rw] <- RN_vec[rw] + sum(v)
            RSX_vec[rw] <- RSX_vec[rw] + sum(rx[v])
            RSXX_vec[rw] <- RSXX_vec[rw] + sum(rx[v]^2)
            RSXY_vec[rw] <- RSXY_vec[rw] + sum(rx[v] * yr[v])
          }
        }
        rm(dt_rep, mr)
      }
    }
    rm(dt_gc, m, sub_logr)
    gc()
  }

  calc_corr <- function(n, sx, sy, sxx, syy, sxy) {
    num <- (n * sxy) - (sx * sy)
    den <- sqrt(pmax(0, (n * sxx - sx^2) * (n * syy - sy^2)))
    return(ifelse(den == 0, 0, num / den))
  }
  corrs <- sapply(win_cols, function(w) unname(abs(calc_corr(N_vec[w], SX_vec[w], SY, SXX_vec[w], SYY, SXY_vec[w]))))

  index_2kb <- which(names(corrs) == "2kb")
  if (length(index_2kb) == 0) index_2kb <- floor(length(corrs) / 2)
  maxGCcol_insert <- names(which.max(corrs[1:index_2kb]))
  maxGCcol_amplic <- names(which.max(corrs[(index_2kb + 1):length(corrs)]))
  index_100kb <- which(names(corrs) == "100kb")
  if (length(index_100kb) > 0 && index_100kb > index_2kb) maxGCcol_amplic <- names(which.max(corrs[(index_2kb + 1):index_100kb]))

  maxreplic <- NULL
  if (has_replic) {
    corrs_rep <- sapply(rep_win_cols, function(w) unname(abs(calc_corr(RN_vec[w], RSX_vec[w], SY, RSXX_vec[w], SYY, RSXY_vec[w]))))
    maxreplic <- names(which.max(corrs_rep))
  }
  log_info("Selected Windows: Insert={maxGCcol_insert}, Amplic={maxGCcol_amplic}, Rep={maxreplic}")

  # Pass 2: Online Linear Regression (Accumulate X'X and X'y)
  log_info("Pass 2: Accumulating matrix cross-products for the spline model...")
  XtX <- NULL
  Xty <- NULL

  for (cn in chrom_names) {
    log_info("  - Pass 2: Processing {cn}...")
    gc_f <- paste0(gc_content_file_prefix, cn, ".txt.gz")
    if (!file.exists(gc_f)) next
    dt_gc <- load_ref_dt(gc_f, gc_info)
    dt_gc[, `:=`(Chromosome = gsub("chr", "", as.character(Chromosome)), Position = as.integer(Position))]

    sub_logr <- all_logr[gsub("chr", "", as.character(cn))]

    data.table::setkey(dt_gc, Position)
    data.table::setkey(sub_logr, Position)
    m <- dt_gc[sub_logr, nomatch = 0]
    log_info("    - Joined for regression: {nrow(m)} SNPs")
    if (nrow(m) == 0) next

    Xi <- cbind(splines::ns(m[[maxGCcol_insert]], df = 5, intercept = TRUE), splines::ns(m[[maxGCcol_amplic]], df = 5, intercept = FALSE))
    if (has_replic) {
      rep_f <- paste0(replic_timing_file_prefix, cn, ".txt.gz")
      dt_rep <- load_ref_dt(rep_f, rep_info)
      dt_rep[, `:=`(Chromosome = gsub("chr", "", as.character(Chromosome)), Position = as.integer(Position))]

      data.table::setkey(dt_rep, Position)
      mr <- dt_rep[m, nomatch = 0]
      Xi <- cbind(Xi, splines::ns(mr[[maxreplic]], df = 5, intercept = FALSE))
      y_i <- as.numeric(mr[[ncol(mr)]])
      rm(dt_rep, mr)
    } else {
      y_i <- as.numeric(m[[ncol(m)]])
    }

    # Remove NAs which break splineDesign/solve
    keep <- rowSums(is.na(Xi)) == 0 & !is.na(y_i)
    if (sum(keep) < 20) {
      rm(dt_gc, m, Xi, y_i)
      next
    }
    Xi <- Xi[keep, , drop = FALSE]
    y_i <- y_i[keep]

    if (is.null(XtX)) {
      n_cols <- ncol(Xi)
      XtX <- matrix(0, n_cols, n_cols)
      Xty <- numeric(n_cols)
    }

    XtX <- XtX + t(Xi) %*% Xi
    Xty <- Xty + t(Xi) %*% y_i
    rm(dt_gc, m, Xi, y_i)
    gc()
  }

  beta <- solve(XtX, Xty)
  log_info("Pass 3: Calculating and writing residuals...")
  if (file.exists(outfile)) file.remove(outfile)

  # Final pass to write results
  for (cn in chrom_names) {
    log_info("  - Pass 3: Writing {cn}...")
    gc_f <- paste0(gc_content_file_prefix, cn, ".txt.gz")
    if (!file.exists(gc_f)) next
    dt_gc <- load_ref_dt(gc_f, gc_info)
    dt_gc[, `:=`(Chromosome = gsub("chr", "", as.character(Chromosome)), Position = as.integer(Position))]

    sub_logr <- all_logr[gsub("chr", "", as.character(cn))]
    data.table::setkey(dt_gc, Position)
    data.table::setkey(sub_logr, Position)
    m <- dt_gc[sub_logr, nomatch = 0]
    log_info("    - Joined for output: {nrow(m)} SNPs")
    if (nrow(m) == 0) next

    Xi <- cbind(splines::ns(m[[maxGCcol_insert]], df = 5, intercept = TRUE), splines::ns(m[[maxGCcol_amplic]], df = 5, intercept = FALSE))
    if (has_replic) {
      rep_f <- paste0(replic_timing_file_prefix, cn, ".txt.gz")
      dt_rep <- load_ref_dt(rep_f, rep_info)
      dt_rep[, `:=`(Chromosome = gsub("chr", "", as.character(Chromosome)), Position = as.integer(Position))]

      data.table::setkey(dt_rep, Position)
      mr <- dt_rep[m, nomatch = 0]
      Xi_rep <- splines::ns(mr[[maxreplic]], df = 5, intercept = FALSE)

      # For output, we apply logic to each row. But since we filtered with joins,
      # we need to be careful. Splines ns() will return NA for rows with NA input.
      # residual = y - X * beta
      # We'll do it in a robust way:
      Xi_full <- cbind(Xi, Xi_rep)
      y_full <- as.numeric(mr[[ncol(mr)]])
      residuals <- y_full - (Xi_full %*% beta)

      out_dt <- mr[, 1:2]
      out_dt$LogR <- as.numeric(residuals)
      rm(dt_rep, mr, Xi_rep, Xi_full)
    } else {
      residuals <- as.numeric(m[[ncol(m)]]) - (Xi %*% beta)
      out_dt <- m[, 1:2]
      out_dt$LogR <- as.numeric(residuals)
    }
    out_dt$LogR <- pmax(pmin(out_dt$LogR, 5), -5)
    data.table::fwrite(out_dt, file = outfile, sep = "\t", append = TRUE, col.names = !file.exists(outfile))
    rm(dt_gc, m, Xi, out_dt)
    gc()
  }
}

#' Prepare WGS data for haplotype construction
#'
#' This function performs part of the Battenberg WGS pipeline: Counting alleles, constructing BAF and logR
#' and performing GC content correction.
#'
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param tumourbam Full path to the tumour BAM file
#' @param normalbam Full path to the normal BAM file
#' @param tumourname Identifier to be used for tumour output files
#' @param normalname Identifier to be used for normal output files
#' @param g1000allelesprefix Prefix path to the 1000 Genomes alleles reference files
#' @param g1000prefix Prefix path to the 1000 Genomes SNP reference files
#' @param gccorrectprefix Prefix path to GC content reference data
#' @param repliccorrectprefix Prefix path to replication timing reference data (supply NULL if no replication timing correction is to be applied)
#' @param min_base_qual Minimum base quality required for a read to be counted
#' @param min_map_qual Minimum mapping quality required for a read to be counted
#' @param allele_counts_dir Directory containing the allele counts files
#' @param min_normal_depth Minimum depth required in the normal for a SNP to be included
#' @param nthreads The number of paralel processes to run
#' @param libs Path to the R libraries to be used by parallel workers
#' @author sd11
#' @export
prepare_wgs <- function(
  chrom_names,
  tumourbam,
  normalbam,
  tumourname,
  normalname,
  g1000allelesprefix,
  g1000prefix,
  gccorrectprefix,
  repliccorrectprefix,
  min_base_qual,
  min_map_qual,
  allele_counts_dir,
  min_normal_depth,
  nthreads,
  libs
) {
  # Check files exist
  tumour_prefix <- file.path(allele_counts_dir, paste0(tumourname, "_alleleFrequencies_chr"))
  normal_prefix <- file.path(allele_counts_dir, paste0(normalname, "_alleleFrequencies_chr"))

  # Simple validation for first chromosome to ensure files are present
  # Note: detailed validation could loop over all chromosomes
  first_tumour_file <- paste0(tumour_prefix, chrom_names[1], ".txt")
  if (!file.exists(first_tumour_file)) {
    log_failure("Expected tumour allele counts file not found: {first_tumour_file}")
  }

  # Obtain BAF and LogR from the raw allele counts
  getBAFsAndLogRs(
    tumourAlleleCountsFile.prefix = tumour_prefix,
    normalAlleleCountsFile.prefix = normal_prefix,
    figuresFile.prefix = paste(tumourname, "_", sep = ""),
    BAFnormalFile = paste(tumourname, "_normalBAF.tab", sep = ""),
    BAFmutantFile = paste(tumourname, "_mutantBAF.tab", sep = ""),
    logRnormalFile = paste(tumourname, "_normalLogR.tab", sep = ""),
    logRmutantFile = paste(tumourname, "_mutantLogR.tab", sep = ""),
    combinedAlleleCountsFile = paste(tumourname, "_alleleCounts.tab", sep = ""),
    chr_names = chrom_names,
    g1000file.prefix = g1000allelesprefix,
    minCounts = min_normal_depth,
    samplename = tumourname
  )
  # Perform GC correction
  gc_correct_wgs(
    Tumour_LogR_file = paste(tumourname, "_mutantLogR.tab", sep = ""),
    outfile = paste(tumourname, "_mutantLogR_gcCorrected.tab", sep = ""),
    correlations_outfile = paste(tumourname, "_GCwindowCorrelations.txt", sep = ""),
    gc_content_file_prefix = gccorrectprefix,
    replic_timing_file_prefix = repliccorrectprefix,
    chrom_names = chrom_names
  )

  log_info("Battenberg WGS preparation complete. Corrected LogR written to: {paste(tumourname, '_mutantLogR_gcCorrected.tab', sep='')}")
}

#' A helper function to split the genome into parts
#' @param SNPpos A data.frame with a row for each SNP. First column is chromosome, second column position
#' @noRd
split_genome <- function(SNPpos) {
  # look for gaps of more than 1Mb and chromosome borders
  holesOver1Mb <- which(diff(SNPpos[, 2]) >= 1000000) + 1
  chrBorders <- which(diff(as.numeric(factor(SNPpos[, 1], levels = unique(SNPpos[, 1])))) != 0) + 1
  holes <- unique(sort(c(holesOver1Mb, chrBorders)))

  # find which segments are too small
  joincandidates <- which(diff(c(0, holes, dim(SNPpos)[1])) < 200)

  # if it's the first or last segment, just join to the one next to it, irrespective of chromosome and positions
  while (1 %in% joincandidates) {
    holes <- holes[-1]
    joincandidates <- which(diff(c(0, holes, dim(SNPpos)[1])) < 200)
  }
  while ((length(holes) + 1) %in% joincandidates) {
    holes <- holes[-length(holes)]
    joincandidates <- which(diff(c(0, holes, dim(SNPpos)[1])) < 200)
  }

  while (length(joincandidates) != 0) {
    # the while loop is because after joining, segments may still be too small..
    startseg <- c(1, holes)
    endseg <- c(holes - 1, dim(SNPpos)[1])

    # for each segment that is too short, see if it has the same chromosome as the segments before and after
    # the next always works because neither the first or the last segment is in joincandidates now
    previoussamechr <- SNPpos[endseg[joincandidates - 1], 1] == SNPpos[startseg[joincandidates], 1]
    nextsamechr <- SNPpos[endseg[joincandidates], 1] == SNPpos[startseg[joincandidates + 1], 1]

    distanceprevious <- SNPpos[startseg[joincandidates], 2] - SNPpos[endseg[joincandidates - 1], 2]
    distancenext <- SNPpos[startseg[joincandidates + 1], 2] - SNPpos[endseg[joincandidates], 2]

    # if both the same, decide based on distance, otherwise if one the same, take the other, if none, just take one.
    joins <- ifelse(previoussamechr & nextsamechr,
      ifelse(distanceprevious > distancenext, joincandidates, joincandidates - 1),
      ifelse(nextsamechr, joincandidates, joincandidates - 1)
    )

    holes <- holes[-joins]
    joincandidates <- which(diff(c(0, holes, dim(SNPpos)[1])) < 200)
  }
  # if two neighboring segments are selected, this may make bigger segments then absolutely necessary.
  startseg <- c(1, holes)
  endseg <- c(holes - 1, dim(SNPpos)[1])
  chr <- list()
  for (i in seq_along(startseg)) {
    chr[[i]] <- startseg[i]:endseg[i]
  }

  return(chr)
}
