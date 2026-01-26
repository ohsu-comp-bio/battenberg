####################################################################################################
#' ASCAT like function to obtain a clonal copy number profile
#'
#' This function takes an initial optimum rho/psi pair and uses
#' an internal distance metric to calculate a score for each rho/psi pair allowed.
#' The solution with the best score is then taken to obtain a global copy number
#' profile. This function performs both a grid search and tries to find a reference
#' segment, but the grid search result is always used for now.
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param baf (unsegmented) B Allele Frequency, in genomic sequence (all probes),
#' with probe IDs
#' @param lrrsegmented log R, segmented, in genomic sequence (all probes), with
#' probe IDs
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only
#' probes heterozygous in germline), with probe IDs
#' @param chromosomes a list containing c vectors, where c is the number of
#' chromosomes and every vector contains all probe numbers per chromosome
#' @param segBAF_table Segmented BAF data.frame from \code{get_segment_info}
#' @param input_optimum_pair A list containing fields for rho, psi and ploidy,
#' as is output from \code{runASCAT}
#' @param dist_choice The distance metric to be used internally to penalise a copy
#' number solution
#' @param distancepng if NA: distance is plotted, if filename is given, the plot
#' is written to a .png file (Default NA)
#' @param copynumberprofilespng if NA: possible copy number profiles are plotted,
#' if filename is given, the plot is written to a .png file (Default NA)
#' @param nonroundedprofilepng if NA: copy number profile before rounding is
#' plotted (total copy number as well as the copy number of the minor allele), if
#' filename is given, the plot is written to a .png file (Default NA)
#' @param gamma_param technology parameter, compaction of Log R profiles (expected
#' decrease in case of deletion in diploid sample, 100 "\%" aberrant cells; 1 in
#' ideal case, 0.55 of Illumina 109K arrays) (Default 0.55)
#' @param read_depth TODO: unused parameter that should be removed
#' @param uninformative_baf_threshold The threshold beyond which BAF becomes
#' uninformative
#' @param allow100percent A boolean whether to allow a 100"\%" cellularity
#' solution
#' @param reliabilityFile String to where fit reliabilty information should be
#' written. This file contains backtransformed BAF and LogR values for segments
#' using the fitted copy number profile (Default NA)
#' @param psi_min_initial Minimum psi value to be considered (Default: 1.0)
#' @param psi_max_initial Maximum psi value to be considered (Default: 5.4)
#' @param rho_min_initial Minimum rho value to be considered (Default: 0.1)
#' @param rho_max_initial Maximum rho value to be considered (Default: 1.05)
#' @param chr_names A vector with chromosome names used for plotting
#' @param nthreads The number of paralel processes to run
#' @return A list with fields output_optimum_pair, output_optimum_pair_without_ref,
#' distance, distance_without_ref, minimise and is_ref_better
#' @export
run_clonal_ASCAT <- function(
  lrr, baf, lrrsegmented,
  bafsegmented, chromosomes,
  segBAF_table, input_optimum_pair,
  dist_choice, distancepng = NA,
  copynumberprofilespng = NA,
  nonroundedprofilepng = NA,
  gamma_param, read_depth,
  uninformative_baf_threshold,
  allow100percent,
  reliabilityFile = NA,
  psi_min_initial = 1.0,
  psi_max_initial = 5.4,
  rho_min_initial = 0.1,
  rho_max_initial = 1.05,
  chr_names,
  nthreads = 1
) {
  siglevel_BAF <- 0.05
  maxdist_BAF <- 0.01

  # DCW 160314 - much more lenient logR thresholds (allow anything!)
  #  # TODO: This parameter is pushed down to is_segment_clonal but not used there (maybe not used at all?)
  siglevel_LogR <- -0.01
  maxdist_LogR <- 1

  initial_bounds <- list(psi_min = psi_min_initial, psi_max = psi_max_initial, rho_min = rho_min_initial, rho_max = rho_max_initial)

  new_bounds <- get_new_bounds(input_optimum_pair, initial_bounds)


  ch <- chromosomes
  b <- bafsegmented
  r <- lrrsegmented[names(bafsegmented)]

  # CRITICAL FIX: Subset LRR using PROBE NAMES (names of bafsegmented)
  # segBAF_table rownames are numeric indices (1..N) which causes mismatch with named lrrsegmented vector
  s <- get_segment_info(lrrsegmented[names(bafsegmented)], segBAF_table)
  log_debug("get_segment_info returned: {nrow(s)} rows, {ncol(s)} columns")
  if (nrow(s) > 0) {
    log_debug("get_segment_info head: {paste(head(s, 1), collapse=', ')}")
  } else {
    log_debug("get_segment_info returned empty matrix")
  }

  if (is.null(s) || nrow(s) == 0) {
    log_failure("No valid segments found in run_clonal_ASCAT. Cannot proceed with clonal copy number fitting.")
  }

  # Make sure no segment of length 1 remains
  s <- s[s[, "length"] > 1, , drop = FALSE]
  log_debug("After filtering length > 1: {nrow(s)} rows")
  if (nrow(s) == 0) {
    log_failure("No segments with length > 1 found in run_clonal_ASCAT.")
  }

  dist_matrix_info <- create_distance_matrix_clonal(
    s, dist_choice, gamma_param, read_depth, siglevel_BAF, maxdist_BAF,
    siglevel_LogR, maxdist_LogR, uninformative_baf_threshold, new_bounds,
    nthreads = nthreads
  ) # kjd 10-2-2013

  d <- dist_matrix_info$distance_matrix
  if (all(is.na(d)) || all(is.infinite(d))) {
    log_failure("Distance matrix is entirely NA or Inf in run_clonal_ASCAT. No valid copy number solution possible.")
  }
  minimise <- dist_matrix_info$minimise

  # DCW 210314
  if (minimise) {
    best.distance <- min(d)
  } else {
    best.distance <- max(d)
  }

  ref_seg_matrix <- dist_matrix_info$ref_seg_matrix

  ref_major <- dist_matrix_info$ref_major
  ref_minor <- dist_matrix_info$ref_minor

  #########################################################

  ret <- find_centroid_of_global_minima(
    d, ref_seg_matrix, ref_major,
    ref_minor, s, dist_choice, minimise,
    new_bounds, distancepng, gamma_param,
    siglevel_BAF, maxdist_BAF, siglevel_LogR,
    maxdist_LogR, allow100percent,
    uninformative_baf_threshold, read_depth
  )
  optima_info_without_ref <- ret$optima_info_without_ref
  optima_info <- ret$optima_info

  nropt <- optima_info$nropt
  psi_opt1 <- optima_info$psi_opt1
  rho_opt1 <- optima_info$rho_opt1
  ploidy_opt1 <- optima_info$ploidy_opt1
  goodness_of_fit_opt1 <- optima_info$goodness_of_fit_opt1

  distance.from.ref.seg <- goodness_of_fit_opt1

  is_ref_better <- FALSE
  if (is.na(rho_opt1)) {
    log_info("reference segment did not provide a possible solution")
  } else if (psi_opt1 >= psi_min_initial && psi_opt1 <= psi_max_initial &&
    rho_opt1 >= rho_min_initial && rho_opt1 <= rho_max_initial &&
    ((minimise && distance.from.ref.seg < best.distance) ||
      (!minimise && distance.from.ref.seg > best.distance))) {
    is_ref_better <- T
    log_info("reference segment gives better results than grid search")
  } else {
    log_info("reference segment gives no better results than grid search. \\
             Reverting to grid search solution")
  }

  psi_without_ref <- optima_info_without_ref$psi_opt1
  rho_without_ref <- optima_info_without_ref$rho_opt1
  ploidy_without_ref <- optima_info_without_ref$ploidy_opt1
  goodness_of_fit_without_ref <- optima_info_without_ref$goodness_of_fit_opt1

  #########################################################

  if (nropt > 0) {
    if (is_ref_better) {
      rho <- rho_opt1
      psi <- psi_opt1
      ploidy <- ploidy_opt1
      goodness_of_fit <- goodness_of_fit_opt1
    } else {
      rho <- rho_without_ref
      psi <- psi_without_ref
      ploidy <- ploidy_without_ref
      goodness_of_fit <- goodness_of_fit_without_ref
    }
    nAfull <- (rho - 1 - (b - 1) * 2^(r / gamma_param) *
      ((1 - rho) * 2 + rho * psi)) / rho
    nBfull <- (rho - 1 + b * 2^(r / gamma_param) *
      ((1 - rho) * 2 + rho * psi)) / rho
    nA <- pmax(round(nAfull), 0)
    nB <- pmax(round(nBfull), 0)

    rBacktransform <- gamma_param *
      log((rho * (nA + nB) + (1 - rho) * 2) / ((1 - rho) * 2 + rho * psi), 2)
    bBacktransform <- (1 - rho + rho * nB) / (2 - 2 * rho + rho * (nA + nB))
    rDiff <- 1 - abs(rBacktransform - r) / abs(r)
    rConf <- ifelse(abs(rBacktransform) > 0.15,
      pmin(100, pmax(0, 100 * rDiff)), NA
    )
    bDiff <- 1 - abs(bBacktransform - b) / abs(b - 0.5)
    bConf <- ifelse(bBacktransform != 0.5,
      pmin(100, pmax(0, ifelse(b == 0.5, 100, 100 * bDiff))), NA
    )
    # DCW 150711 - get deviations from expected values
    if (!is.na(reliabilityFile)) {
      data.table::fwrite(
        data.frame(
          segmentedBAF = b, backTransformedBAF = bBacktransform,
          confidenceBAF = bConf, segmentedR = r,
          backTransformedR = rBacktransform, confidenceR = rConf,
          nA = nA, nB = nB, nAfull = nAfull, nBfull = nBfull
        ),
        reliabilityFile,
        sep = ",", row.names = FALSE
      )
    }

    # SMART DOWNSAMPLING for performance
    log_info("Applying chromosome-aware smart downsampling to plotting data...")
    target_total <- 500000
    total_probes <- length(lrr)
    lrr_list <- vector("list", length(ch))
    baf_list <- vector("list", length(ch))
    nA_list <- vector("list", length(ch))
    nB_list <- vector("list", length(ch))
    nAfull_list <- vector("list", length(ch))
    nBfull_list <- vector("list", length(ch))
    ch_ds <- vector("list", length(ch))
    curr_pos <- 1

    for (i in seq_along(ch)) {
      idx <- ch[[i]]
      if (length(idx) == 0) next
      chr_target <- max(500, round(target_total * length(idx) / total_probes))
      keep_rel <- bt_downsample_indices(lrr[idx], chr_target)
      keep_abs <- idx[keep_rel]

      lrr_list[[i]] <- lrr[keep_abs]
      baf_list[[i]] <- bafsegmented[keep_abs]
      nA_list[[i]] <- nA[keep_abs]
      nB_list[[i]] <- nB[keep_abs]
      nAfull_list[[i]] <- nAfull[keep_abs]
      nBfull_list[[i]] <- nBfull[keep_abs]

      new_len <- length(keep_abs)
      ch_ds[[i]] <- seq(curr_pos, length.out = new_len)
      curr_pos <- curr_pos + new_len
    }

    lrr_ds <- unlist(lrr_list)
    bafsegmented_ds <- unlist(baf_list)
    nA_ds <- unlist(nA_list)
    nB_ds <- unlist(nB_list)
    nAfull_ds <- unlist(nAfull_list)
    nBfull_ds <- unlist(nBfull_list)
    if (!is.null(names(ch))) names(ch_ds) <- names(ch)

    # Make plots in parallel if requested
    plot_tasks <- list()

    if (!is.na(copynumberprofilespng)) {
      plot_tasks[["profile"]] <- function() {
        grDevices::png(
          filename = copynumberprofilespng,
          width = 2000, height = 500,
          res = 200, type = "cairo"
        )
        ASCAT::ascat.plotAscatProfile(
          n1all = nA_ds, n2all = nB_ds,
          heteroprobes = TRUE,
          ploidy = ploidy, rho = rho,
          goodnessOfFit = goodness_of_fit * 100,
          nonaberrant = FALSE,
          ch = ch_ds, lrr = lrr_ds,
          bafsegmented = bafsegmented_ds,
          chrs = chr_names
        )
        grDevices::dev.off()
      }
    }

    if (!is.na(nonroundedprofilepng)) {
      plot_tasks[["nonrounded"]] <- function() {
        grDevices::png(
          filename = nonroundedprofilepng,
          width = 2000, height = 500,
          res = 200, type = "cairo"
        )
        ASCAT::ascat.plotNonRounded(
          ploidy = ploidy, rho = rho,
          goodnessOfFit = goodness_of_fit * 100,
          nonaberrant = FALSE, nAfull = nAfull_ds,
          nBfull = nBfull_ds, bafsegmented = bafsegmented_ds,
          ch = ch_ds, lrr = lrr_ds, chrs = chr_names
        )
        grDevices::dev.off()
      }
    }

    if (length(plot_tasks) > 0) {
      log_info("Generating {length(plot_tasks)} genome-wide plots sequentially to ensure container stability...")
      lapply(plot_tasks, function(f) f())
    }
  }

  # Recalculate the psi_t for this rho using only clonal segments
  psi_t <- recalc_psi_t(
    psi_without_ref, rho_without_ref, gamma_param, r, segBAF_table,
    siglevel_BAF, maxdist_BAF,
    include_subcl_segments = FALSE
  )

  # If there aren't any clonally fit segments, the above yields NA. In this case, revert to the original grid search psi_t
  if (is.na(psi_t)) {
    log_info("Recalculated psi_t was NA, reverting to grid search solution. This occurs when no segment could be fit with a clonal state, check sample for contamination")
    psi_t <- psi_without_ref
  }

  output_optimum_pair <- list(psi = psi_opt1, rho = rho_opt1, ploidy = ploidy_opt1)
  # output_optimum_pair_without_ref = list(psi = psi_without_ref, rho = rho_without_ref, ploidy = ploidy_without_ref)
  # Use the recalculated psi_t from the clonal segments as our final estimate
  # of psi_t which is data driven with rho fixed
  output_optimum_pair_without_ref <- list(
    psi = psi_t, rho = rho_without_ref, ploidy = ploidy_without_ref
  )
  return(
    list(
      output_optimum_pair = output_optimum_pair,
      output_optimum_pair_without_ref = output_optimum_pair_without_ref,
      distance = goodness_of_fit_opt1,
      distance_without_ref = goodness_of_fit_without_ref,
      minimise = minimise,
      is_ref_better = is_ref_better,
      dist_matrix_info = dist_matrix_info
    )
  )
}

#' Function extends the ASCAT \code{make_segments} function to make segments
#' of constant BAF and LogR. This function returns a matrix with for each
#' segment the LogR, BAF, the length of the segment (twice), and the mean and
#' standard deviation of the BAF values
#' @noRd
get_segment_info <- function(segLogR, segBAF_table) {
  # Column names for robust access
  col_names <- names(segBAF_table)
  baf_col <- if ("BAFseg" %in% col_names) "BAFseg" else 5
  phased_col <- if ("BAFphased" %in% col_names) "BAFphased" else 4

  b_raw_full <- segBAF_table[[baf_col]]
  b_phased_full <- segBAF_table[[phased_col]]

  # 1. Consensus filtering: ensures r, b, and phased BAF are aligned and non-NA
  # This matches the internal logic of make_segments for consistency
  valid_mask <- !is.na(segLogR) & !is.na(b_raw_full)
  r <- segLogR[valid_mask]
  b <- b_raw_full[valid_mask]
  bp <- b_phased_full[valid_mask]

  if (length(r) == 0) {
    return(matrix(0, 0, 7))
  }

  # 2. Identify contiguous segments using Run-Length Encoding ID
  # This is the defining logic of a segment: contiguous regions with same values
  seg_id <- data.table::rleid(r, b)

  # 3. Aggregate stats per segment
  # This avoids any indexing mismatch errors and handles the matrix creation in one pass
  dt <- data.table::data.table(r = r, b = b, bp = bp, seg_id = seg_id)

  # We need 7 columns: r, b, length, length.1, size, mean, sd
  # Battenberg legacy format repeats length/size columns
  stats_dt <- dt[, .(
    r = .subset2(r, 1),
    b = .subset2(b, 1),
    len1 = .N,
    len2 = .N,
    size = .N,
    mean_bp = mean(bp, na.rm = TRUE),
    sd_bp = sd(bp, na.rm = TRUE)
  ), by = seg_id]

  # 4. Return as matrix with exact column names expected by Battenberg
  res <- as.matrix(stats_dt[, .(r, b, len1, len2, size, mean_bp, sd_bp)])
  colnames(res) <- c("r", "b", "length", "length.1", "size", "mean", "sd")

  return(res)
}


#' Optimized Segment Maker - Returns 3 columns like ASCAT original
#' @noRd
make_segments <- function(r, b) {
  m <- matrix(ncol = 2, nrow = length(b))
  m[, 1] <- r
  m[, 2] <- b
  m <- as.matrix(na.omit(m))

  if (nrow(m) == 0) {
    return(matrix(nrow = 0, ncol = 3, dimnames = list(NULL, c("r", "b", "length"))))
  }

  pcf_segments <- matrix(ncol = 3, nrow = dim(m)[1])
  colnames(pcf_segments) <- c("r", "b", "length")

  index <- 0
  previousb <- -1
  previousr <- 1E10

  for (i in seq_len(dim(m)[1])) {
    # Use a small tolerance for floating point comparisons to ensure segmented values collapse correctly
    if (abs(m[i, 2] - previousb) > 1e-10 || abs(m[i, 1] - previousr) > 1e-10) {
      index <- index + 1
      count <- 1
      pcf_segments[index, "r"] <- m[i, 1]
      pcf_segments[index, "b"] <- m[i, 2]
    } else {
      count <- count + 1
    }
    pcf_segments[index, "length"] <- count
    previousb <- m[i, 2]
    previousr <- m[i, 1]
  }

  # Clean up the matrix to remove unused pre-allocated rows
  pcf_segments <- pcf_segments[seq_len(index), , drop = FALSE]
  return(pcf_segments)
}
