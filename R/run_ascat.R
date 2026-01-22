#' A modified ASCAT main function to fit Battenberg
#'
#' This function returns an initial rho and psi estimate for a clonal copy number fit. It uses an internal distance metric to create a distance matrix.
#' Using that matrix it will search for a rho and psi combination that yields the least heavy penalty.
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param baf (unsegmented) B Allele Frequency, in genomic sequence (all probes), with probe IDs
#' @param lrrsegmented log R, segmented, in genomic sequence (all probes), with probe IDs
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param chromosomes a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#' @param dist_choice The distance metric to be used internally to penalise a copy number solution
#' @param distancepng if NA: distance is plotted, if filename is given, the plot is written to a .png file (Default NA)
#' @param copynumberprofilespng if NA: possible copy number profiles are plotted, if filename is given, the plot is written to a .png file (Default NA)
#' @param nonroundedprofilepng if NA: copy number profile before rounding is plotted (total copy number as well as the copy number of the minor allele), if filename is given, the plot is written to a .png file (Default NA)
#' @param cnaStatusFile File where the copy number profile status is written to. This contains either the message "No suitable copy number solution found" or "X copy number solutions found" (Default copynumber_solution_status.txt)
#' @param gamma technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 "\%" aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays) (Default 0.55)
#' @param allow100percent A boolean whether to allow a 100"\%" cellularity solution
#' @param reliabilityFile String to where fit reliabilty information should be written. This file contains backtransformed BAF and LogR values for segments using the fitted copy number profile (Default NA)
#' @param min_ploidy The minimum ploidy to consider (Default 1.6)
#' @param max_ploidy The maximum ploidy to consider (Default 4.8)
#' @param min_rho The minimum cellularity to consider (Default 0.1)
#' @param max_rho The maximum cellularity to consider (Default 1.0)
#' @param min_goodness The minimum goodness of fit for a solution to have to be considered (Default 63)
#' @param uninformative_baf_threshold The threshold beyond which BAF becomes uninformative (Default 0.51)
#' @param chr_names A vector with chromosome names used for plotting
#' @param analysis A String representing the type of analysis to be run, this determines whether the distance figure is produced (Default paired)
#' @param nthreads The number of paralel processes to run
#' @return A list with fields psi, rho and ploidy
#' @export
# the limit on rho is lenient and may lead to spurious solutions
runASCAT <- function(
  lrr, baf, lrrsegmented,
  bafsegmented, chromosomes,
  dist_choice, distancepng = NA,
  copynumberprofilespng = NA,
  nonroundedprofilepng = NA,
  cnaStatusFile = "copynumber_solution_status.txt",
  gamma = 0.55, allow100percent,
  reliabilityFile = NA, min_ploidy = 1.6,
  max_ploidy = 4.8, min_rho = 0.1,
  max_rho = 1.0, min_goodness = 63,
  uninformative_baf_threshold = 0.51,
  chr_names, analysis = "paired",
  nthreads = 1
) {
  # Setup inputs and segments
  ch <- chromosomes
  b <- bafsegmented
  r <- lrrsegmented[names(bafsegmented)]

  # Adapt the rho/psi boundaries
  dist_min_psi <- max(min_ploidy - 0.6, 0)
  dist_max_psi <- max_ploidy + 0.6
  dist_min_rho <- max(min_rho - 0.1, 0.05)
  dist_max_rho <- max_rho + 0.1

  s <- make_segments(r, b)
  dist_matrix_info <- create_distance_matrix(
    s, dist_choice, gamma,
    uninformative_baf_threshold = uninformative_baf_threshold,
    min_psi = dist_min_psi,
    max_psi = dist_max_psi,
    min_rho = dist_min_rho,
    max_rho = dist_max_rho,
    nthreads = nthreads
  )
  d <- dist_matrix_info$distance_matrix
  if (all(is.na(d)) || all(is.infinite(d))) {
    log_failure("Distance matrix is entirely NA or Inf in runASCAT. No valid copy number solution possible.")
  }
  minimise <- dist_matrix_info$minimise

  # Calculate theoretical max distance for goodness of fit
  TheoretMaxdist <- sum(rep(0.25, dim(s)[1]) * s[, "length"], na.rm = TRUE)
  total_len <- sum(s[, "length"])

  # Ensure we are always searching for a minimum
  if (!minimise) d <- -d

  # VECTORIZED LOCAL MINIMA SEARCH (Pixel-perfect replacement for 7x7 loop)
  nr <- nrow(d)
  nc <- ncol(d)
  is_local_min <- matrix(TRUE, nrow = nr, ncol = nc)

  # Constrain search to the interior to match 4:(dim-3) logic
  row_range <- 4:(nr - 3)
  col_range <- 4:(nc - 3)

  # Check every neighbor in the 7x7 window (48 neighbors)
  for (dx in -3:3) {
    for (dy in -3:3) {
      if (dx == 0 && dy == 0) next
      is_local_min[row_range, col_range] <- is_local_min[row_range, col_range] &
        (d[row_range, col_range] < d[row_range + dx, col_range + dy])
    }
  }

  # Zero out the margins to match original loop boundaries
  is_local_min[-row_range, ] <- FALSE
  is_local_min[, -col_range] <- FALSE

  # Extraction helper to process candidates
  evaluate_candidates <- function(indices, current_d) {
    if (nrow(indices) == 0) {
      return(NULL)
    }

    # Pre-calculate segment masks for efficiency
    is_not_balanced <- s[, "b"] != 0.5
    weight_unbalanced <- sum(s[, "length"] * is_not_balanced)

    results <- apply(indices, 1, function(idx) {
      i <- idx[1]
      j <- idx[2]
      m <- current_d[i, j]
      psi <- as.numeric(rownames(current_d)[i])
      rho <- as.numeric(colnames(current_d)[j])

      # Copy number algebra
      common_term <- 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)
      nA <- (rho - 1 - (s[, "b"] - 1) * common_term) / rho
      nB <- (rho - 1 + s[, "b"] * common_term) / rho

      ploidy <- sum((nA + nB) * s[, "length"]) / total_len

      # Biological viability checks
      is_nA_zero <- round(nA) == 0
      is_nB_zero <- round(nB) == 0
      percentzero <- (sum(is_nA_zero * s[, "length"]) + sum(is_nB_zero * s[, "length"])) / total_len
      perczeroAbb <- (sum(is_nA_zero * s[, "length"] * is_not_balanced) + sum(is_nB_zero * s[, "length"] * is_not_balanced)) / weight_unbalanced
      if (is.na(perczeroAbb)) perczeroAbb <- 0

      # Goodness of fit calculation
      fit <- if (minimise) (1 - m / TheoretMaxdist) * 100 else -m / TheoretMaxdist * 100

      # Return data if it meets primary constraints (percentzero checks applied later if allow100percent is used)
      return(list(m = m, i = i, j = j, ploidy = ploidy, fit = fit, pz = percentzero, pza = perczeroAbb, rho = rho, psi = psi))
    })
    return(results)
  }

  # First pass: find optima meeting the percentzero conditions
  opt_indices <- which(is_local_min, arr.ind = TRUE)
  candidates <- evaluate_candidates(opt_indices, d)

  # Debug stats container
  debug_stats <- list(
    ploidy_bounds = 0,
    rho_bounds = 0,
    low_goodness = 0,
    zero_constraint = 0
  )

  # Filtering based on standard Battenberg criteria with logging
  valid_optima <- list()
  if (!is.null(candidates)) {
    valid_optima <- Filter(function(x) {
      if (x$ploidy < min_ploidy || x$ploidy > max_ploidy) {
        debug_stats$ploidy_bounds <<- debug_stats$ploidy_bounds + 1
        return(FALSE)
      }
      if (x$rho < min_rho) {
        debug_stats$rho_bounds <<- debug_stats$rho_bounds + 1
        return(FALSE)
      }
      if (x$fit < min_goodness) {
        debug_stats$low_goodness <<- debug_stats$low_goodness + 1
        return(FALSE)
      }
      if (!(x$pz > 0.01 || x$pza > 0.1)) {
        debug_stats$zero_constraint <<- debug_stats$zero_constraint + 1
        return(FALSE)
      }
      return(TRUE)
    }, candidates)
  }

  # Second pass: If allow100percent is TRUE and no solutions found, relax constraints
  if (allow100percent && length(valid_optima) == 0) {
    # Penalize cellularity > 1 as per original code
    cold_idx <- which(as.numeric(colnames(d)) > 1)
    d[, cold_idx] <- 1e20

    # Reset debug stats for second pass (optional, or keep cumulative)
    # Re-evaluate all local minima with relaxed biological constraints
    valid_optima <- Filter(function(x) {
      return(x$ploidy > min_ploidy && x$ploidy < max_ploidy &&
        x$rho >= min_rho && x$fit >= min_goodness)
    }, candidates)
  }

  # Process the winning solution
  nropt <- length(valid_optima)
  psi_opt1_plot <- vector(mode = "numeric")
  rho_opt1_plot <- vector(mode = "numeric")

  if (nropt > 0) {
    data.table::fwrite(
      list(paste(nropt, " copy number solutions found", sep = "")),
      file = cnaStatusFile, quote = FALSE, col.names = FALSE, row.names = FALSE
    )

    # Find the global minimum among the local optima
    all_m <- sapply(valid_optima, function(x) x$m)
    optlim <- min(all_m)

    # Extract ties for plotting and set the final result
    for (opt in valid_optima) {
      if (opt$m == optlim) {
        psi_opt1 <- opt$psi
        rho_opt1 <- min(opt$rho, 1)
        ploidy_opt1 <- opt$ploidy
        goodness_of_fit_opt1 <- opt$fit

        psi_opt1_plot <- c(psi_opt1_plot, psi_opt1)
        rho_opt1_plot <- c(rho_opt1_plot, rho_opt1)
      }
    }

    log_info("DEBUG: After filtering, {nropt} valid solutions remain")
    log_info("DEBUG: Selected solution: rho={round(rho_opt1, 3)}, psi={round(psi_opt1, 3)}, ploidy={round(ploidy_opt1, 3)}, goodness={round(goodness_of_fit_opt1, 2)}%")
  } else {
    writeLines("no copy number solutions found", con = cnaStatusFile)
    log_info("No suitable copy number solution found.")
    log_info("Debug Rejection Stats: PloidyBounds={debug_stats$ploidy_bounds}, RhoBounds={debug_stats$rho_bounds}, LowGoodness={debug_stats$low_goodness}, ZeroConstraint={debug_stats$zero_constraint}")
    psi <- ploidy <- rho <- NA
    psi_opt1_plot <- rho_opt1_plot <- -1
  }

  # Plotting Sunrise (if paired) - Delayed to run in parallel with other plots
  # (Logic moved to plotting section below)

  # Final calculations for the best solution
  if (nropt > 0) {
    rho <- rho_opt1
    psi <- psi_opt1
    ploidy <- ploidy_opt1

    # Full genomic fit
    # Full genomic fit
    # Optimized Back-transformation with data.table chunking
    # This matches the enhanced version's logic for speed and memory efficiency
    log_info("Starting back-transformation (Chunked execution, threads={nthreads})...")

    indices <- seq_along(r)
    num_chunks <- max(1, nthreads)
    chunks <- parallel::splitIndices(length(indices), num_chunks)

    results <- parallel::mclapply(chunks, function(idx) {
      b_sub <- b[idx]
      r_sub <- r[idx]

      # Calculate mult locally
      mult_sub <- 2^(r_sub / gamma) * ((1 - rho) * 2 + rho * psi)

      nAfull_sub <- (rho - 1 - (b_sub - 1) * mult_sub) / rho
      nBfull_sub <- (rho - 1 + b_sub * mult_sub) / rho
      nA_sub <- pmax(round(nAfull_sub), 0)
      nB_sub <- pmax(round(nBfull_sub), 0)

      rBT_sub <- gamma * log(
        (rho * (nA_sub + nB_sub) + (1 - rho) * 2) / ((1 - rho) * 2 + rho * psi),
        2
      )
      bBT_sub <- (1 - rho + rho * nB_sub) / (2 - 2 * rho + rho * (nA_sub + nB_sub))

      return(data.table::data.table(
        segmentedBAF = b_sub, backTransformedBAF = bBT_sub, segmentedR = r_sub,
        backTransformedR = rBT_sub, nA = nA_sub, nB = nB_sub, nAfull = nAfull_sub,
        nBfull = nBfull_sub
      ))
    }, mc.cores = nthreads)

    log_info("Aggregating results...")
    final_dt <- data.table::rbindlist(results)

    # Extract variables for standard plotting/usage downstream
    nA <- final_dt$nA
    nB <- final_dt$nB
    nAfull <- final_dt$nAfull
    nBfull <- final_dt$nBfull
    rBT <- final_dt$backTransformedR
    bBT <- final_dt$backTransformedBAF

    if (!is.na(reliabilityFile)) {
      # Use threaded writing
      data.table::fwrite(
        list(
          segmentedBAF = b, backTransformedBAF = bBT, segmentedR = r,
          backTransformedR = rBT, nA = nA, nB = nB, nAfull = nAfull,
          nBfull = nBfull
        ),
        reliabilityFile,
        sep = ",", row.names = FALSE,
        nThread = nthreads
      )
    }

    # Generate Profile Plots in Parallel
    plot_tasks <- list()

    if (analysis == "paired" && !is.na(distancepng)) {
      plot_tasks[["sunrise"]] <- function() {
        # Recalculate res based on original logic (1000/7 approx 142.8)
        grDevices::png(filename = distancepng, width = 1000, height = 1000, res = 1000 / 7, type = "cairo")
        ASCAT::ascat.plotSunrise(-d, psi_opt1_plot, rho_opt1_plot, minimise)
        grDevices::dev.off()
      }
    }

    if (!is.na(copynumberprofilespng)) {
      plot_tasks[["profile"]] <- function() {
        grDevices::png(
          filename = copynumberprofilespng,
          width = 2000, height = 500,
          res = 200, type = "cairo"
        )
        ASCAT::ascat.plotAscatProfile(
          n1all = nA, n2all = nB, heteroprobes = TRUE,
          ploidy = ploidy_opt1, rho = rho_opt1,
          goodnessOfFit = goodness_of_fit_opt1,
          nonaberrant = FALSE, ch = ch,
          lrr = lrr, bafsegmented = bafsegmented,
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
          ploidy = ploidy_opt1, rho = rho_opt1,
          goodnessOfFit = goodness_of_fit_opt1,
          nonaberrant = FALSE, nAfull = nAfull,
          nBfull = nBfull, bafsegmented = bafsegmented,
          ch = ch, lrr = lrr, chrs = chr_names
        )
        grDevices::dev.off()
      }
    }

    if (length(plot_tasks) > 0) {
      if (nthreads > 1 && length(plot_tasks) > 1) {
        log_info("Generating plots in parallel (threads={min(nthreads, length(plot_tasks))})...")
        parallel::mclapply(plot_tasks, function(f) f(), mc.cores = min(nthreads, length(plot_tasks)))
      } else {
        lapply(plot_tasks, function(f) f())
      }
    }
  }

  return(list(psi = psi, rho = rho, ploidy = ploidy))
}
