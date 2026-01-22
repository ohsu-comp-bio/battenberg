#' Key optimizations:
#' 1. Early termination after first good solution (like original)
#' 2. Vectorized distance calculations
#' 3. Optimized constraint checking
#' 4. Smart search ordering (best regions first)
#' 5. Reduced memory allocations
#' @export
runASCAT_enhanced <- function(
  lrr, baf, lrrsegmented, bafsegmented, chromosomes, dist_choice,
  distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA,
  cnaStatusFile = "copynumber_solution_status.txt", gamma = 0.55,
  allow100percent, reliabilityFile = NA, min_ploidy = 1.6, max_ploidy = 4.8,
  min_rho = 0.1, max_rho = 1.0, min_goodness = 63,
  uninformative_baf_threshold = 0.51, chr_names, analysis = "paired",
  smart_ordering = TRUE, early_termination = TRUE, verbose = TRUE, nthreads = 1
) {
  start_time <- Sys.time()
  log_info("BATTENBERG ASCAT ENHANCED - VERSION CHECK: FAILSAVE & SD-FIX APPLIED !!!")

  # 0. Input Validation
  if (missing(lrr) || missing(baf) || missing(lrrsegmented) || missing(bafsegmented)) {
    log_failure("Missing required input arguments for runASCAT_enhanced")
    stop("Missing input arguments")
  }

  if (!is.numeric(lrr) || length(lrr) == 0) log_failure("Invalid lrr: must be numeric and non-empty")
  if (!is.numeric(baf) || length(baf) == 0) log_failure("Invalid baf: must be numeric and non-empty")
  if (!is.numeric(lrrsegmented) || length(lrrsegmented) == 0) log_failure("Invalid lrrsegmented: must be numeric and non-empty")
  if (!is.numeric(bafsegmented) || length(bafsegmented) == 0) log_failure("Invalid bafsegmented: must be numeric and non-empty")

  if (length(lrrsegmented) != length(bafsegmented)) {
    log_failure("Length mismatch: lrrsegmented ({length(lrrsegmented)}) != bafsegmented ({length(bafsegmented)})")
    stop("Input length mismatch")
  }


  # 1. Setup Data Processing
  ch <- chromosomes
  b <- bafsegmented
  # Use direct assignment - names(bafsegmented) is often NULL which empties r
  logR_segmented <- lrrsegmented
  if (length(logR_segmented) != length(b)) {
    log_failure("Length mismatch in runASCAT_enhanced: LRR {length(logR_segmented)} vs BAF {length(b)}")
  }

  dist_min_psi <- max(min_ploidy - 0.6, 0)
  dist_max_psi <- max_ploidy + 0.6
  dist_min_rho <- max(min_rho - 0.03, 0.05)
  dist_max_rho <- max_rho + 0.03

  # 2. Create Segments & Distance Matrix
  s <- make_segments(logR_segmented, b)

  log_info("Number of segments created: {nrow(s)}")
  # ADD THESE DEBUG LINES:
  log_info("DEBUG: Segment matrix dimensions: {nrow(s)} x {ncol(s)}")
  log_info("DEBUG: Column names: {paste(colnames(s), collapse=', ')}")
  log_info("DEBUG: First few rows:")
  print(head(s, 10))
  log_info("DEBUG: Sum of segment lengths: {sum(s[,'length'])}")
  log_info("DEBUG: Length of input b: {length(b)}")


  log_info("Number of segments created: {nrow(s)}")

  if (nrow(s) == 0) {
    log_failure("No valid segments created in runASCAT_enhanced. Cannot proceed with grid search.")
  }

  dist_matrix_info <- create_distance_matrix(s, dist_choice, gamma,
    uninformative_baf_threshold = uninformative_baf_threshold,
    min_psi = dist_min_psi, max_psi = dist_max_psi,
    min_rho = dist_min_rho, max_rho = dist_max_rho,
    nthreads = nthreads
  )
  d <- dist_matrix_info$distance_matrix

  # Theoretical maximum distance (weighted by length)
  TheoretMaxdist <- collapse::fsum(rep(0.25, nrow(s)) * s[, "length"],
    na.rm = TRUE
  )

  minimise <- dist_matrix_info$minimise

  log_debug("Distance matrix dimensions: {nrow(d)} x: {ncol(d)}")
  log_debug("Theoretical Max Distance: {round(TheoretMaxdist, 4)}")

  log_info("DEBUG: Distance matrix stats BEFORE negation: min={min(d, na.rm=TRUE)}, max={max(d, na.rm=TRUE)}, mean={mean(d, na.rm=TRUE)}")

  if (!minimise) d <- -d

  # 3. Pre-compute Search Parameters
  rho_values <- as.numeric(colnames(d))
  psi_values <- as.numeric(rownames(d))
  s_length <- s[, "length"]
  s_b <- s[, "b"]
  s_r <- s[, "r"]
  total_length <- collapse::fsum(s_length)

  # Pre-compute masks for calculate_solution_fast
  baf_mask <- s_b != 0.5
  denom_abb <- collapse::fsum(s_length[baf_mask])

  # 3.1 Vectorized Local Minima Detection
  # We use the original Battenberg logic: a point is a local minimum if it is STRICTLY LESS
  # than all other points in its 7x7 neighborhood.
  nr <- nrow(d)
  nc <- ncol(d)
  is_local_min <- matrix(TRUE, nrow = nr, ncol = nc)

  # Constrain to interior 4:(dim-3)
  row_range <- 4:(nr - 3)
  col_range <- 4:(nc - 3)

  # Fill with FALSE for safety, only interior can be TRUE
  is_local_min[, ] <- FALSE
  is_local_min[row_range, col_range] <- TRUE

  # Check neighbors
  for (dx in -3:3) {
    for (dy in -3:3) {
      if (dx == 0 && dy == 0) next
      # Use is_local_min & (...) and handle NAs by treating them as larger than any value
      # This ensures NAs don't invalidate the whole mask
      neighbor_vals <- d[row_range + dx, col_range + dy]
      neighbor_vals[is.na(neighbor_vals)] <- Inf # NAs are not minima

      comparison <- (d[row_range, col_range] < neighbor_vals)
      comparison[is.na(comparison)] <- FALSE
      is_local_min[row_range, col_range] <- is_local_min[row_range, col_range] & comparison
    }
  }

  # Get search matrix (i, j)
  search_order <- create_smart_search_order(d, smart_ordering, verbose)
  total_points_in_grid <- nrow(search_order)

  # Failsafe: If no strict local minima found, we MUST check the full grid
  # as per the fallback logic in the original runASCAT.
  if (sum(is_local_min, na.rm = TRUE) == 0 && total_points_in_grid > 0) {
    if (verbose) log_info("No strict local minima found. Activating FULL GRID search...")
    is_local_min[row_range, col_range] <- TRUE
  }


  # 4. Main Search Loop
  nropt <- 0
  optima <- list()
  localmin_vals <- numeric()
  points_checked <- 0

  # Debug stats
  debug_stats <- list(
    pre_check_bounds = 0,
    ploidy_bounds = 0,
    low_goodness = 0,
    zero_constraint = 0,
    max_goodness = -1
  )

  if (total_points_in_grid > 0) {
    # Pre-calculate max possible goodness
    min_dist <- min(d, na.rm = TRUE)
    max_poss_goodness <- if (minimise) (1 - min_dist / TheoretMaxdist) * 100 else -min_dist / TheoretMaxdist * 100
    log_info("DEBUG START SEARCH: Min Dist={min_dist}, Max Possible Goodness={round(max_poss_goodness, 2)}% (Threshold: {min_goodness}%)")

    if (verbose) log_info("Starting grid search over {total_points_in_grid} points...")
    for (idx in seq_len(total_points_in_grid)) {
      i <- search_order[idx, 1]
      j <- search_order[idx, 2]

      # Use the pre-computed mask
      if (!is_local_min[i, j]) next

      m <- d[i, j]
      points_checked <- points_checked + 1

      solution <- calculate_solution_fast(
        psi_values[i], rho_values[j], s_b, s_r, s_length, total_length,
        gamma, min_ploidy, max_ploidy, min_rho, max_rho,
        min_goodness, m, TheoretMaxdist, minimise, allow100percent,
        baf_mask = baf_mask, denom_abb = denom_abb
      )

      if (solution$valid) {
        nropt <- nropt + 1
        # Store as vector for consistency with original optima extraction
        optima[[nropt]] <- c(m, i, j, solution$ploidy, solution$goodness)
        localmin_vals[nropt] <- m

        if (verbose) {
          log_info("Found solution {nropt} at point {points_checked}: rho={round(rho_values[j], 3)}, psi={round(psi_values[i], 3)}")
        }

        if (early_termination && solution$goodness >= (min_goodness + 5)) {
          if (verbose) log_info("Early termination triggered: Good solution found.")
          break
        }
      } else {
        # Track rejection reason
        reject_reason <- solution$reason
        if (!is.null(reject_reason)) {
          debug_stats[[reject_reason]] <- debug_stats[[reject_reason]] + 1
        }
        if (!is.null(solution$goodness) && solution$goodness > debug_stats$max_goodness) {
          debug_stats$max_goodness <- solution$goodness
        }
      }

      # Correctly report progress inside the loop
      if (verbose && (points_checked %% 1000 == 0 || points_checked == total_points_in_grid)) {
        pct_val <- round(points_checked / total_points_in_grid * 100, 1)
        log_info("Progress: {points_checked}/{total_points_in_grid} ({pct_val}%) points checked")
      }
    }
  }


  # 5. Handle 100% Aberrant Fallback
  if (allow100percent && nropt == 0) {
    log_info("DEBUG FIRST PASS FAILED: Rejected: pre_check={debug_stats$pre_check_bounds}, ploidy_bounds={debug_stats$ploidy_bounds}, low_goodness={debug_stats$low_goodness}, zero_constraint={debug_stats$zero_constraint}")
    log_info("DEBUG FIRST PASS FAILED: Max Goodness found: {round(debug_stats$max_goodness, 2)}")

    if (verbose) log_info("Trying 100% aberrant solutions...")
    d_mod <- d
    d_mod[, rho_values <= 1] <- 1e20
    search_order_100 <- create_smart_search_order(d_mod, smart_ordering, FALSE)

    # Pre-compute local minima for d_mod (interior only)
    is_local_min_mod <- matrix(FALSE, nrow = nr, ncol = nc)
    if (nr >= 7 && nc >= 7) {
      is_local_min_mod[row_range, col_range] <- TRUE
      for (dx in -3:3) {
        for (dy in -3:3) {
          if (dx == 0 && dy == 0) next
          is_local_min_mod[row_range, col_range] <- is_local_min_mod[row_range, col_range] &
            (d_mod[row_range, col_range] <= d_mod[row_range + dx, col_range + dy])
        }
      }
    }

    if (nrow(search_order_100) > 0) {
      for (idx in seq_len(nrow(search_order_100))) {
        i <- search_order_100[idx, 1]
        j <- search_order_100[idx, 2]

        if (!is_local_min_mod[i, j]) next

        m <- d_mod[i, j]
        solution <- calculate_solution_fast(
          psi_values[i], rho_values[j], s_b, s_r, s_length, total_length, gamma,
          min_ploidy, max_ploidy, min_rho, max_rho,
          min_goodness, m, TheoretMaxdist, minimise, allow100percent,
          baf_mask = baf_mask, denom_abb = denom_abb,
          skip_zero_check = TRUE # RELAX CONSTRAINTS FOR FALLBACK
        )
        if (solution$valid) {
          nropt <- 1
          optima[[1]] <- c(m, i, j, solution$ploidy, solution$goodness)
          localmin_vals[1] <- m
          break # Stop after finding first valid solution in fallback mode
        }
      }
    }
  }


  optimization_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Select Best Solution & Collect Sunrise Plot Data
  if (nropt > 0) {
    data.table::fwrite(list(paste0(nropt, " copy number solutions found")), cnaStatusFile)

    optlim <- sort(localmin_vals)[1]
    psi_opt1_plot <- numeric()
    rho_opt1_plot <- numeric()

    # Original logic: collect all solutions that share the global minimum distance
    for (idx in seq_along(optima)) {
      if (optima[[idx]][1] == optlim) {
        psi_opt1 <- psi_values[optima[[idx]][2]]
        rho_opt1 <- min(rho_values[optima[[idx]][3]], 1.0)
        ploidy_opt1 <- optima[[idx]][4]
        goodness_of_fit_opt1 <- optima[[idx]][5]

        psi_opt1_plot <- c(psi_opt1_plot, psi_opt1)
        rho_opt1_plot <- c(rho_opt1_plot, rho_opt1)
      }
    }
  } else {
    data.table::fwrite(list("no copy number solutions found"), cnaStatusFile)

    log_failure("ASCAT Optimization failed. Rejected: pre_check={debug_stats$pre_check_bounds}, ploidy_bounds={debug_stats$ploidy_bounds}, low_goodness={debug_stats$low_goodness}, zero_constraint={debug_stats$zero_constraint}")
    if (debug_stats$max_goodness > -1) {
      log_failure("Best rejected candidate had goodness: {round(debug_stats$max_goodness, 2)} (threshold: {min_goodness}). If this is high, check ploidy/zero constraints.")
    }

    return(list(
      psi = NA, rho = NA, ploidy = NA,
      convergence_info = list(
        converged = FALSE, n_solutions_found = 0,
        optimization_time = optimization_time, points_checked = points_checked,
        search_efficiency = points_checked / total_points_in_grid
      )
    ))
  }

  # Use the extracted "best" values for the final vectors
  rho <- rho_opt1
  psi <- psi_opt1
  ploidy <- ploidy_opt1

  # 7. Final Back-transformation
  log_info("Debug Backtransform: rho={rho}, psi={psi}, length(logR_segmented)={length(logR_segmented)}, class={class(logR_segmented)}, gamma={gamma}")
  if (!is.numeric(logR_segmented)) {
    log_failure("CRITICAL: logR_segmented corrupted. Value: {paste(head(logR_segmented), collapse=', ')}")
  }

  # Always use chunked execution to manage memory and provide consistent logging
  # Even with nthreads=1, this prevents massive single-step allocations
  log_info("Starting back-transformation (Chunked execution, threads={nthreads})...")

  indices <- seq_along(logR_segmented)
  # Ensure at least 1 chunk
  num_chunks <- max(1, nthreads)
  chunks <- parallel::splitIndices(length(indices), num_chunks)

  results <- parallel::mclapply(chunks, function(idx) {
    # Extract subset
    r_sub <- logR_segmented[idx]
    b_sub <- b[idx]

    # Calculate mult locally to save memory
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

    # Reliability
    # Handle potentially empty r_sub
    if (length(r_sub) > 0) {
      rDiff <- 1 - abs(rBT_sub - r_sub) / abs(r_sub)
      rConf_sub <- ifelse(abs(rBT_sub) > 0.15, pmin(100, pmax(0, 100 * rDiff)), NA)

      bDiff <- 1 - abs(bBT_sub - b_sub) / abs(b_sub - 0.5)
      bConf_sub <- ifelse(bBT_sub != 0.5,
        pmin(100, pmax(0, ifelse(b_sub == 0.5, 100, 100 * bDiff))), NA
      )
    } else {
      rConf_sub <- numeric(0)
      bConf_sub <- numeric(0)
    }

    # Return as a data.table chunk for fast rbindlist
    return(data.table::data.table(
      segmentedBAF = b_sub,
      backTransformedBAF = bBT_sub,
      confidenceBAF = bConf_sub,
      segmentedR = r_sub,
      backTransformedR = rBT_sub,
      confidenceR = rConf_sub,
      nA = nA_sub,
      nB = nB_sub,
      nAfull = nAfull_sub,
      nBfull = nBfull_sub
    ))
  }, mc.cores = nthreads)

  # Fast aggregation
  log_info("Aggregating results...")
  start_agg <- Sys.time()
  final_dt <- data.table::rbindlist(results)
  log_info(paste("Aggregation complete in", round(difftime(Sys.time(), start_agg, units = "secs"), 2), "seconds"))

  if (!is.na(reliabilityFile)) {
    # Optimization: Write the prepared data.table directly
    log_info(paste("Writing reliability file to", reliabilityFile, "..."))
    start_write <- Sys.time()
    # Use threaded writing if available
    data.table::fwrite(
      final_dt,
      reliabilityFile,
      sep = ",", row.names = FALSE,
      nThread = nthreads
    )
    log_info(paste("Writing complete in", round(difftime(Sys.time(), start_write, units = "secs"), 2), "seconds"))
  }

  # Extract vectors for plotting (plotting functions expect these variable names)
  nA <- final_dt$nA
  nB <- final_dt$nB
  nAfull <- final_dt$nAfull
  nBfull <- final_dt$nBfull
  # Ensure these are numeric vectors
  if (is.null(nA)) log_failure("Critical: nA missing from results")

  # 8. Plotting
  # Define plotting tasks as closures
  plot_tasks <- list()

  # SMART DOWNSAMPLING for performance
  # Target ~100k points across the whole genome
  # We downsample each chromosome to preserve original indexing mapping in 'ch'
  log_info("Applying chromosome-aware smart downsampling to plotting data...")

  # helper to find min/max indices in a vector segment
  get_keep_indices <- function(v, target) {
    n <- length(v)
    if (n <= target) {
      return(seq_along(v))
    }
    bin_size <- ceiling(n / (target / 2))
    dt_ds <- data.table::data.table(val = as.numeric(v), id = seq_along(v))
    dt_ds[, bin := ceiling(id / bin_size)]
    keep <- dt_ds[, .(id_min = id[which.min(val)], id_max = id[which.max(val)]), by = bin]
    return(sort(unique(c(keep$id_min, keep$id_max))))
  }

  target_total <- 100000
  total_probes <- length(lrr)

  # Accumulate in lists to avoid O(N^2) overhead
  lrr_list <- vector("list", length(ch))
  baf_list <- vector("list", length(ch))
  nA_list <- vector("list", length(ch))
  nB_list <- vector("list", length(ch))
  nAfull_list <- vector("list", length(ch))
  nBfull_list <- vector("list", length(ch))
  ch_ds <- vector("list", length(ch))

  curr_pos <- 1
  start_ds <- Sys.time()

  for (i in seq_along(ch)) {
    idx <- ch[[i]]
    if (length(idx) == 0) next

    # Proportionate target for this chromosome
    chr_target <- max(500, round(target_total * length(idx) / total_probes))

    # Relies on data.table for speed
    keep_rel <- get_keep_indices(lrr[idx], chr_target)
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

  # Flatten lists once
  lrr_ds <- unlist(lrr_list)
  bafsegmented_ds <- unlist(baf_list)
  nA_ds <- unlist(nA_list)
  nB_ds <- unlist(nB_list)
  nAfull_ds <- unlist(nAfull_list)
  nBfull_ds <- unlist(nBfull_list)

  log_info("Downsampling complete in {round(difftime(Sys.time(), start_ds, units='secs'), 2)} seconds. Reduced to {length(lrr_ds)} points.")
  # Preserve names for plotter consistency if they exist
  if (!is.null(names(ch))) names(ch_ds) <- names(ch)

  if (analysis == "paired" && !is.na(distancepng)) {
    plot_tasks[["sunrise"]] <- function() {
      log_info("SUNRISE: Starting calculation for {distancepng}...")
      log_info("SUNRISE DEBUG: d matrix stats - min={min(d, na.rm=TRUE)}, max={max(d, na.rm=TRUE)}, NA_count={sum(is.na(d))}")
      log_info("SUNRISE DEBUG: psi_opt1_plot length={length(psi_opt1_plot)}, rho_opt1_plot length={length(rho_opt1_plot)}")
      if (length(psi_opt1_plot) > 0) {
        log_info("SUNRISE DEBUG: first sol: rho={rho_opt1_plot[1]}, psi={psi_opt1_plot[1]}")
      }
      t1 <- Sys.time()
      grDevices::png(filename = distancepng, width = 1000, height = 1000, res = 150, type = "cairo")
      ASCAT::ascat.plotSunrise(-d, psi_opt1_plot, rho_opt1_plot, minimise)
      grDevices::dev.off()
      t2 <- Sys.time()
      log_info("SUNRISE: Finished in {round(difftime(t2, t1, units='secs'), 2)}s")
    }
  }

  if (!is.na(copynumberprofilespng)) {
    plot_tasks[["profile"]] <- function() {
      log_info("PROFILE: Starting genome-wide plot (probes={length(lrr_ds)})...")
      t1 <- Sys.time()
      grDevices::png(filename = copynumberprofilespng, width = 2000, height = 500, res = 200, type = "cairo")
      ASCAT::ascat.plotAscatProfile(
        n1all = nA_ds, n2all = nB_ds, heteroprobes = TRUE, ploidy = ploidy,
        rho = rho, goodnessOfFit = goodness_of_fit_opt1, nonaberrant = FALSE,
        ch = ch_ds, lrr = lrr_ds, bafsegmented = bafsegmented_ds, chrs = chr_names
      )
      grDevices::dev.off()
      t2 <- Sys.time()
      log_info("PROFILE: Finished in {round(difftime(t2, t1, units='secs'), 2)}s")
    }
  }

  if (!is.na(nonroundedprofilepng)) {
    plot_tasks[["nonrounded"]] <- function() {
      log_info("NONROUNDED: Starting genome-wide plot (probes={length(lrr_ds)})...")
      t1 <- Sys.time()
      grDevices::png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200, type = "cairo")
      ASCAT::ascat.plotNonRounded(
        ploidy = ploidy, rho = rho, goodnessOfFit = goodness_of_fit_opt1,
        nonaberrant = FALSE, nAfull = nAfull_ds, nBfull = nBfull_ds,
        bafsegmented = bafsegmented_ds, ch = ch_ds, lrr = lrr_ds, chrs = chr_names
      )
      grDevices::dev.off()
      t2 <- Sys.time()
      log_info("NONROUNDED: Finished in {round(difftime(t2, t1, units='secs'), 2)}s")
    }
  }

  if (length(plot_tasks) > 0) {
    if (nthreads > 1 && length(plot_tasks) > 1 && .Platform$OS.type != "windows") {
      n_workers <- min(nthreads, length(plot_tasks))
      log_info("Generating {length(plot_tasks)} plots in parallel (FORK, threads={n_workers})...")

      # Use mclapply for high-performance forking
      # This is much faster than PSOCK as it avoids copying the downsampled data
      parallel::mclapply(plot_tasks, function(f) f(), mc.cores = n_workers)
    } else {
      log_info("Generating {length(plot_tasks)} plots sequentially...")
      lapply(plot_tasks, function(f) f())
    }
    log_info("All plotting tasks completed.")
  }

  return(list(
    psi = psi, rho = rho, ploidy = ploidy,
    convergence_info = list(
      converged = TRUE,
      n_solutions_found = nropt,
      optimization_time = optimization_time,
      points_checked = points_checked,
      search_efficiency = points_checked / total_points_in_grid
    )
  ))
}

#' Create search order for grid search
#' CRITICAL: Must use row-major order (i then j) to match original battenberg/R/grid_search.R:258
#' @noRd
create_smart_search_order <- function(d, smart_ordering, verbose) {
  nr <- nrow(d)
  nc <- ncol(d)

  # CRITICAL: Use row-major order (i then j) to match original
  # Original: for(i in 4:(nr-3)) for(j in 4:(nc-3))
  # This is psi-first order, NOT rho-first from which()
  search_points <- list()

  if (nr >= 7 && nc >= 7) {
    # Match original border exclusion
    for (i in 4:(nr - 3)) {
      for (j in 4:(nc - 3)) {
        if (is.finite(d[i, j])) {
          search_points[[length(search_points) + 1]] <- list(i = i, j = j, distance = d[i, j])
        }
      }
    }
  } else {
    # Small matrix fallback
    idx_mat <- which(is.finite(d), arr.ind = TRUE)
    for (k in seq_len(nrow(idx_mat))) {
      i <- idx_mat[k, 1]
      j <- idx_mat[k, 2]
      search_points[[length(search_points) + 1]] <- list(i = i, j = j, distance = d[i, j])
    }
  }

  if (length(search_points) == 0) {
    return(matrix(0, 0, 2))
  }

  if (smart_ordering) {
    distances <- sapply(search_points, function(p) p$distance)
    search_points <- search_points[order(distances)]
  }

  # Convert to matrix
  result <- matrix(0, nrow = length(search_points), ncol = 2)
  for (k in seq_along(search_points)) {
    result[k, 1] <- search_points[[k]]$i
    result[k, 2] <- search_points[[k]]$j
  }

  return(result)
}

#' Fast solution calculation (vectorized and optimized)
calculate_solution_fast <- function(
  psi, rho, s_b, s_r, s_length, total_length, gamma,
  min_ploidy, max_ploidy, min_rho, max_rho,
  min_goodness, distance_value, TheoretMaxdist, minimise,
  allow100percent, baf_mask, denom_abb, skip_zero_check = FALSE
) {
  # Guard against rho = 0 to prevent Inf
  safe_rho <- pmax(rho, 1e-6)

  # Constraint pre-check
  if (psi < min_ploidy || psi > max_ploidy || rho < min_rho || rho > max_rho) {
    return(list(valid = FALSE, reason = "pre_check_bounds"))
  }

  # Vectorized calculation
  multiplier <- 2^(s_r / gamma) * ((1 - safe_rho) * 2 + safe_rho * psi)
  nA <- (safe_rho - 1 - (s_b - 1) * multiplier) / safe_rho
  nB <- (safe_rho - 1 + s_b * multiplier) / safe_rho

  # Ploidy check
  ploidy <- collapse::fsum((nA + nB) * s_length) / total_length
  if (is.na(ploidy) || ploidy < min_ploidy || ploidy > max_ploidy) {
    return(list(valid = FALSE, reason = "ploidy_bounds", ploidy = ploidy))
  }

  # Goodness check
  goodness_of_fit <- if (minimise) {
    (1 - distance_value / TheoretMaxdist) * 100
  } else {
    -distance_value / TheoretMaxdist * 100
  }
  if (is.na(goodness_of_fit) || goodness_of_fit < min_goodness) {
    return(list(valid = FALSE, reason = "low_goodness", goodness = goodness_of_fit))
  }

  if (!skip_zero_check && !allow100percent) {
    nA_r <- round(nA)
    nB_r <- round(nB)
    # Edge case: sum(s_length[logical]) can be 0 if no indices match
    percentzero <- (collapse::fsum(s_length[which(nA_r == 0)]) +
      collapse::fsum(s_length[which(nB_r == 0)])) / total_length

    perczeroAbb <- 0
    if (denom_abb > 0) {
      # Use which() to avoid NA issues in logical indexing
      # Use which() to avoid NA issues in logical indexing
      perczeroAbb <- (collapse::fsum(s_length[which(baf_mask & nA_r == 0)]) +
        collapse::fsum(s_length[which(baf_mask & nB_r == 0)])) /
        denom_abb
    }
    # Ensure we don't have NAs or empty results in our proportions
    if (length(percentzero) == 0 || is.na(percentzero)) percentzero <- 0
    if (length(perczeroAbb) == 0 || is.na(perczeroAbb)) perczeroAbb <- 0

    if (!isTRUE(percentzero > 0.01 || perczeroAbb > 0.1)) {
      return(list(valid = FALSE, reason = "zero_constraint", goodness = goodness_of_fit))
    }
  }

  return(list(valid = TRUE, psi = psi, rho = min(rho, 1.0), ploidy = ploidy, goodness = goodness_of_fit))
}
