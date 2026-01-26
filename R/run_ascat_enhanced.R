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
  min_rho = 0.1, max_rho = 1.0, min_goodness = 0.63,
  uninformative_baf_threshold = 0.51, chr_names, analysis = "paired",
  smart_ordering = TRUE, early_termination = FALSE, verbose = TRUE,
  n_neighbors_search = NULL, psi_step = 0.05, rho_step = 0.01,
  local_min_window_size = 7, nthreads = 1
) {
  start_time <- Sys.time()

  # 0. Input Validation
  if (missing(lrr) || missing(baf) || missing(lrrsegmented) || missing(bafsegmented)) {
    log_failure("Missing required input arguments for runASCAT_enhanced")
    stop("Missing input arguments")
  }

  # Validate new parameters
  if (!is.numeric(local_min_window_size) || local_min_window_size < 3 || local_min_window_size %% 2 == 0) {
    log_failure("local_min_window_size must be an odd integer >= 3, got: {local_min_window_size}")
    stop("Invalid local_min_window_size")
  }

  if (!is.null(n_neighbors_search)) {
    if (!is.numeric(n_neighbors_search) || (!is.infinite(n_neighbors_search) && n_neighbors_search < 1)) {
      log_failure("n_neighbors_search must be NULL, a positive integer, or Inf, got: {n_neighbors_search}")
      stop("Invalid n_neighbors_search")
    }
  }

  # 1. Setup Data Processing
  ch <- chromosomes
  b <- bafsegmented
  # CRITICAL FIX: Match original logic - subset LRR to match BAF (heterozygous probes)
  # The refactor used the full LRR vector which caused segment misalignment and garbage results
  if (!is.null(names(bafsegmented))) {
    logR_segmented <- lrrsegmented[names(bafsegmented)]
  } else {
    # Fallback if names are missing (should not happen in standard pipeline)
    log_info("names(bafsegmented) is NULL. Assuming lrrsegmented and bafsegmented are already aligned or this will fail.")
    logR_segmented <- lrrsegmented
  }

  if (length(logR_segmented) != length(b)) {
    log_failure("Length mismatch in runASCAT_enhanced: LRR {length(logR_segmented)} vs BAF {length(b)}")
  }

  dist_min_psi <- max(min_ploidy - 0.6, 0)
  dist_max_psi <- max_ploidy + 0.6
  dist_min_rho <- max(min_rho - 0.03, 0.05)
  dist_max_rho <- max_rho + 0.03

  # 2. Create Segments & Distance Matrix
  # Use internal tolerance-based make_segments to handle floating point jitter
  s <- make_segments_internal(logR_segmented, b)
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

  log_info("Distance matrix stats: min={min(d, na.rm=TRUE)}, max={max(d, na.rm=TRUE)}, mean={mean(d, na.rm=TRUE)}")
  # We handle minimization/maximization explicitly in the search functions.

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
  nr <- nrow(d)
  nc <- ncol(d)
  is_local_min <- matrix(TRUE, nrow = nr, ncol = nc)
  half_window <- (local_min_window_size - 1) / 2
  row_range <- (half_window + 1):(nr - half_window)
  col_range <- (half_window + 1):(nc - half_window)
  is_local_min[, ] <- FALSE

  if (!is.null(n_neighbors_search)) {
    if (verbose) {
      if (is.infinite(n_neighbors_search)) {
        log_info("Search Mode: Exhaustive search (all grid points)")
      } else {
        log_info("Search Mode: Top {n_neighbors_search} neighbors by distance")
      }
    }
    is_local_min[row_range, col_range] <- TRUE
  } else {
    if (verbose) log_info("Search Mode: Local minima only (window size: {local_min_window_size})")
    is_local_min[row_range, col_range] <- TRUE
    for (dx in -half_window:half_window) {
      for (dy in -half_window:half_window) {
        if (dx == 0 && dy == 0) next
        neighbor_vals <- d[row_range + dx, col_range + dy]

        if (minimise) {
          neighbor_vals[is.na(neighbor_vals)] <- Inf
          comparison <- (d[row_range, col_range] < neighbor_vals)
        } else {
          neighbor_vals[is.na(neighbor_vals)] <- -Inf
          comparison <- (d[row_range, col_range] > neighbor_vals)
        }

        comparison[is.na(comparison)] <- FALSE
        is_local_min[row_range, col_range] <- is_local_min[row_range, col_range] & comparison
      }
    }
  }

  search_order <- create_smart_search_order(
    d, smart_ordering, verbose, minimise,
    local_min_window_size = local_min_window_size,
    skip_local_min = !is.null(n_neighbors_search)
  )
  total_points_in_grid <- nrow(search_order)

  # Apply top N filtering if specified
  if (!is.null(n_neighbors_search) && !is.infinite(n_neighbors_search)) {
    if (total_points_in_grid > n_neighbors_search) {
      # Search order is already sorted by distance (best first)
      # Just take the top N
      search_order <- search_order[1:n_neighbors_search, , drop = FALSE]
      total_points_in_grid <- n_neighbors_search
      if (verbose) log_info("Limited search to top {n_neighbors_search} neighbors")
    }
  }

  # Log how many local minima detected by each method
  num_vectorized_minima <- sum(is_local_min, na.rm = TRUE)
  log_info("Vectorized detection found {num_vectorized_minima} local minima")
  log_info("Smart search order returns {total_points_in_grid} points")

  # Check specific grid point (psi=4.45, rho=0.74) if it exists
  target_psi <- 4.45
  target_rho <- 0.74
  psi_idx <- which.min(abs(psi_values - target_psi))
  rho_idx <- which.min(abs(rho_values - target_rho))
  if (length(psi_idx) > 0 && length(rho_idx) > 0) {
    actual_psi <- psi_values[psi_idx]
    actual_rho <- rho_values[rho_idx]


    # Show window to see why it's not a local min
    if (psi_idx >= (half_window + 1) && psi_idx <= (nr - half_window) &&
      rho_idx >= (half_window + 1) && rho_idx <= (nc - half_window)) {
      window_vals <- d[
        (psi_idx - half_window):(psi_idx + half_window),
        (rho_idx - half_window):(rho_idx + half_window)
      ]
      center_val <- d[psi_idx, rho_idx]
      min_neighbor <- min(window_vals[window_vals != center_val], na.rm = TRUE)
    }
  }

  # Failsafe: If no strict local minima found, we MUST check the full grid
  # as per the fallback logic in the original runASCAT.
  if (sum(is_local_min, na.rm = TRUE) == 0 && total_points_in_grid > 0 && is.null(n_neighbors_search)) {
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
    max_poss_goodness <- if (minimise) {
      min_dist <- min(d, na.rm = TRUE)
      (1 - min_dist / TheoretMaxdist)
    } else {
      max_sim <- max(d, na.rm = TRUE)
      max_sim / TheoretMaxdist
    }
    log_info("Start Search: Optimal Grid Value={if(minimise) min(d, na.rm=TRUE) else max(d, na.rm=TRUE)}, Max Possible Goodness={round(max_poss_goodness * 100, 2)}% (Threshold: {round(min_goodness * 100, 2)}%)")

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
        min_goodness, m, TheoretMaxdist, minimise,
        allow100percent = FALSE, # FIRST PASS ALWAYS REQUIRES LOH/DELETIONS
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

        if (early_termination && solution$goodness >= (min_goodness + 0.05)) {
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

    # 5. Handle 100% Aberrant Fallback
    if (allow100percent && nropt == 0) {
      log_info("DEBUG FIRST PASS FAILED: Rejected: pre_check={debug_stats$pre_check_bounds}, ploidy_bounds={debug_stats$ploidy_bounds}, low_goodness={debug_stats$low_goodness}, zero_constraint={debug_stats$zero_constraint}")
      log_info("DEBUG FIRST PASS FAILED: Max Goodness found: {round(debug_stats$max_goodness, 2)}")

      if (verbose) log_info("Trying 100% aberrant solutions...")
      d_mod <- d
      if (minimise) {
        d_mod[, rho_values > 1] <- 1e20 # Bad for distance
      } else {
        d_mod[, rho_values > 1] <- -1e20 # Bad for similarity
      }

      # CONSISTENCY FIX: Use the same search strategy as first pass
      # If n_neighbors_search was specified, use it for fallback too
      if (!is.null(n_neighbors_search)) {
        # Use top-N search (same as first pass)
        search_order_100 <- create_smart_search_order(d_mod, smart_ordering, FALSE, minimise,
          local_min_window_size = local_min_window_size,
          skip_local_min = TRUE
        )

        # Apply top-N filtering if needed
        fallback_search_limit <- if (!is.infinite(n_neighbors_search)) {
          min(n_neighbors_search, nrow(search_order_100))
        } else {
          nrow(search_order_100)
        }

        if (nrow(search_order_100) > fallback_search_limit) {
          search_order_100 <- search_order_100[1:fallback_search_limit, , drop = FALSE]
        }

        if (verbose) log_info("100% fallback: Searching top {nrow(search_order_100)} points (same as first pass)")

        # Search all points in the order (no local min filtering)
        if (nrow(search_order_100) > 0) {
          for (idx in seq_len(nrow(search_order_100))) {
            i <- search_order_100[idx, 1]
            j <- search_order_100[idx, 2]

            m <- d_mod[i, j]
            solution <- calculate_solution_fast(
              psi_values[i], rho_values[j], s_b, s_r, s_length, total_length, gamma,
              min_ploidy, max_ploidy, min_rho, max_rho,
              min_goodness, m, TheoretMaxdist, minimise, allow100percent,
              baf_mask = baf_mask, denom_abb = denom_abb,
              skip_zero_check = TRUE # RELAX CONSTRAINTS FOR FALLBACK
            )
            if (solution$valid) {
              nropt <- nropt + 1
              optima[[nropt]] <- c(m, i, j, solution$ploidy, solution$goodness)
              localmin_vals[nropt] <- m
            }
          }
        }
      } else {
        # Original local minima search (when n_neighbors_search is NULL)
        search_order_100 <- create_smart_search_order(d_mod, smart_ordering, FALSE, minimise)

        # Pre-compute local minima for d_mod (interior only)
        is_local_min_mod <- matrix(FALSE, nrow = nr, ncol = nc)
        if (nr >= local_min_window_size && nc >= local_min_window_size) {
          is_local_min_mod[row_range, col_range] <- TRUE
          for (dx in -half_window:half_window) {
            for (dy in -half_window:half_window) {
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
              nropt <- nropt + 1
              optima[[nropt]] <- c(m, i, j, solution$ploidy, solution$goodness)
              localmin_vals[nropt] <- m
            }
          }
        }
      }
    }


    optimization_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    # Select Best Solution & Collect Sunrise Plot Data
    if (nropt > 0) {
      data.table::fwrite(list(paste0(nropt, " copy number solutions found")), cnaStatusFile)

      # IMPLMENTATION OF ORIGINAL CENTROID LOGIC
      # Original Battenberg does NOT just take the best goodness.
      # It calculates the "geometric center" of all valid solutions and picks the one closest to it.

      # 1. Extract Grid Coordinates
      grid_x_vect <- sapply(optima, function(z) psi_values[z[2]]) # Psi
      grid_y_vect <- sapply(optima, function(z) rho_values[z[3]]) # Rho

      # 2. Calculate Centroid (Median of means? Original code says: mean(median(grid_x_vect)))
      # This seems redundant (mean of a scalar median is just the median), but we follow it exactly.
      centre_x <- mean(stats::median(grid_x_vect))
      centre_y <- mean(stats::median(grid_y_vect))
      centre <- c(centre_x, centre_y)

      # 3. Find optimum closest to centroid
      best_idx <- 1
      min_sq_dist <- Inf

      # Function to calculate Euclidean distance squared
      calc_sq_dist <- function(p1, p2) {
        sum((p1 - p2)^2)
      }

      for (i in seq_along(optima)) {
        grid_point <- c(psi_values[optima[[i]][2]], rho_values[optima[[i]][3]])
        sq_dist <- calc_sq_dist(grid_point, centre)

        if (sq_dist <= min_sq_dist) {
          min_sq_dist <- sq_dist
          best_idx <- i
        }
      }

      # 4. Extract Winner
      psi_opt1 <- psi_values[optima[[best_idx]][2]]
      rho_opt1 <- min(rho_values[optima[[best_idx]][3]], 1.0)
      ploidy_opt1 <- optima[[best_idx]][4]
      goodness_of_fit_opt1 <- optima[[best_idx]][5] # This is now the clonal genomic proportion (0-1)

      # 5. Collect Plotting Data (All points passing filters)
      psi_opt1_plot <- grid_x_vect
      rho_opt1_plot <- grid_y_vect
    } else {
      data.table::fwrite(list("no copy number solutions found"), cnaStatusFile)

      log_info("ASCAT Optimization failed. Rejected: pre_check={debug_stats$pre_check_bounds}, ploidy_bounds={debug_stats$ploidy_bounds}, low_goodness={debug_stats$low_goodness}, zero_constraint={debug_stats$zero_constraint}")
      if (debug_stats$max_goodness > -1) {
        log_info("Best rejected candidate had goodness: {round(debug_stats$max_goodness * 100, 2)}% (threshold: {round(min_goodness * 100, 2)}%). If this is high, check ploidy/zero constraints.")
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
    log_info("Backtransform: rho={rho}, psi={psi}, length(logR_segmented)={length(logR_segmented)}, class={class(logR_segmented)}, gamma={gamma}")
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

    results <- bt_mclapply(chunks, function(idx) {
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
    log_info("Applying chromosome-aware smart downsampling to plotting data...")

    target_total <- 500000
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
        log_info("Sunrise Plot: Starting calculation for {distancepng}...")
        log_info("Sunrise Plot: d matrix stats - min={min(d, na.rm=TRUE)}, max={max(d, na.rm=TRUE)}, NA_count={sum(is.na(d))}")
        log_info("Sunrise Plot: psi_opt1_plot length={length(psi_opt1_plot)}, rho_opt1_plot length={length(rho_opt1_plot)}")
        if (length(psi_opt1_plot) > 0) {
          log_info("Sunrise Plot: first sol: rho={rho_opt1_plot[1]}, psi={psi_opt1_plot[1]}")
        }

        # Construct bounds for the plot
        psi_values <- as.numeric(rownames(d))
        rho_values <- as.numeric(colnames(d))
        new_bounds <- list(
          psi_min = min(psi_values),
          psi_max = max(psi_values),
          rho_min = min(rho_values),
          rho_max = max(rho_values)
        )

        t1 <- Sys.time()
        tryCatch(
          {
            grDevices::png(filename = distancepng, width = 1000, height = 1000, res = 150, type = "cairo")
            # Use internal plotting function instead of ASCAT::ascat.plotSunrise which is unstable
            clonal_findcentroid_plot(minimise, dist_choice, d, psi_opt1_plot, rho_opt1_plot, new_bounds)
            grDevices::dev.off()
          },
          error = function(e) {
            log_failure("CRITICAL ERROR: Failed to create Sunrise plot at {distancepng}. Error: {e$message}")
            stop(paste("Serious Plotting Error:", e$message))
          }
        )
        t2 <- Sys.time()
        log_info("Sunrise: Finished in {round(difftime(t2, t1, units='secs'), 2)}s")
      }
    }

    if (!is.na(copynumberprofilespng)) {
      plot_tasks[["profile"]] <- function() {
        log_info("Profile Plot: Starting genome-wide plot (probes={length(lrr_ds)})...")
        t1 <- Sys.time()
        grDevices::png(filename = copynumberprofilespng, width = 2000, height = 500, res = 200, type = "cairo")
        ASCAT::ascat.plotAscatProfile(
          n1all = nA_ds, n2all = nB_ds, heteroprobes = TRUE, ploidy = ploidy,
          rho = rho, goodnessOfFit = goodness_of_fit_opt1 * 100, nonaberrant = FALSE,
          ch = ch_ds, lrr = lrr_ds, bafsegmented = bafsegmented_ds, chrs = chr_names
        )
        grDevices::dev.off()
        t2 <- Sys.time()
        log_info("Profile Plot: Finished in {round(difftime(t2, t1, units='secs'), 2)}s")
      }
    }

    if (!is.na(nonroundedprofilepng)) {
      plot_tasks[["nonrounded"]] <- function() {
        log_info("Nonrounded Plot: Starting genome-wide plot (probes={length(lrr_ds)})...")
        t1 <- Sys.time()
        grDevices::png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200, type = "cairo")
        ASCAT::ascat.plotNonRounded(
          ploidy = ploidy, rho = rho, goodnessOfFit = goodness_of_fit_opt1 * 100,
          nonaberrant = FALSE, nAfull = nAfull_ds, nBfull = nBfull_ds,
          bafsegmented = bafsegmented_ds, ch = ch_ds, lrr = lrr_ds, chrs = chr_names
        )
        grDevices::dev.off()
        t2 <- Sys.time()
        log_info("Nonrounded Plot: Finished in {round(difftime(t2, t1, units='secs'), 2)}s")
      }
    }

    if (length(plot_tasks) > 0) {
      log_info("Generating {length(plot_tasks)} genome-wide plots sequentially to ensure container stability...")
      lapply(plot_tasks, function(f) f())
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
}


#' Create search order for grid search
create_smart_search_order <- function(d, smart_ordering, verbose, minimise, local_min_window_size = 7, skip_local_min = FALSE) {
  nr <- nrow(d)
  nc <- ncol(d)
  search_points <- list()
  half_window <- (local_min_window_size - 1) / 2

  if (skip_local_min) {
    for (i in (half_window + 1):(nr - half_window)) {
      for (j in (half_window + 1):(nc - half_window)) {
        m <- d[i, j]
        if (is.finite(m)) {
          search_points[[length(search_points) + 1]] <- list(i = i, j = j, distance = m)
        }
      }
    }
  } else if (nr >= local_min_window_size && nc >= local_min_window_size) {
    for (i in (half_window + 1):(nr - half_window)) {
      for (j in (half_window + 1):(nc - half_window)) {
        m <- d[i, j]
        if (is.finite(m)) {
          seld <- d[(i - half_window):(i + half_window), (j - half_window):(j + half_window)]
          center_idx <- half_window + 1

          if (minimise) {
            # Find local minima
            seld[center_idx, center_idx] <- max(seld, na.rm = TRUE) + 1
            if (min(seld, na.rm = TRUE) > m) search_points[[length(search_points) + 1]] <- list(i = i, j = j, distance = m)
          } else {
            # Find local maxima
            seld[center_idx, center_idx] <- min(seld, na.rm = TRUE) - 1
            if (max(seld, na.rm = TRUE) < m) search_points[[length(search_points) + 1]] <- list(i = i, j = j, distance = m)
          }
        }
      }
    }
  } else {
    idx_mat <- which(is.finite(d), arr.ind = TRUE)
    for (k in seq_len(nrow(idx_mat))) {
      search_points[[length(search_points) + 1]] <- list(i = idx_mat[k, 1], j = idx_mat[k, 2], distance = d[idx_mat[k, 1], idx_mat[k, 2]])
    }
  }

  if (length(search_points) == 0) {
    return(matrix(0, 0, 2))
  }
  if (smart_ordering) {
    distances <- sapply(search_points, function(p) p$distance)
    if (minimise) {
      search_points <- search_points[order(distances)]
    } else {
      search_points <- search_points[order(distances, decreasing = TRUE)]
    }
  }
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

  # Ploidy calculation
  ploidy <- collapse::fsum((nA + nB) * s_length) / total_length

  # Goodness check (cap at 1.0 to prevent overflow)
  goodness_of_fit <- pmin(1.0, if (minimise) {
    (1 - distance_value / TheoretMaxdist)
  } else {
    distance_value / TheoretMaxdist
  })

  if (is.na(goodness_of_fit) || goodness_of_fit < min_goodness) {
    return(list(valid = FALSE, reason = "low_goodness", goodness = goodness_of_fit, ploidy = ploidy))
  }

  if (is.na(ploidy) || ploidy < min_ploidy || ploidy > max_ploidy) {
    return(list(valid = FALSE, reason = "ploidy_bounds", ploidy = ploidy, goodness = goodness_of_fit))
  }

  if (!skip_zero_check && !allow100percent) {
    # Battenberg heuristic: valid solutions usually have at least some segments with CN=0
    # (Loss of Heterozygosity or deletion). Solutions with NO losses are often mathematical
    # artifacts of high-ploidy fits.
    # However, if allow100percent is TRUE, we relax this as the sample might actually have no losses.
    nA_r <- round(nA)
    nB_r <- round(nB)
    # Edge case: sum(s_length[logical]) can be 0 if no indices match
    percentzero <- (collapse::fsum(s_length[which(nA_r == 0)]) +
      collapse::fsum(s_length[which(nB_r == 0)])) / total_length

    perczeroAbb <- 0
    if (denom_abb > 0) {
      # Use which() to avoid NA issues in logical indexing
      perczeroAbb <- (collapse::fsum(s_length[which(baf_mask & nA_r == 0)]) +
        collapse::fsum(s_length[which(baf_mask & nB_r == 0)])) /
        denom_abb
    }
    # Ensure we don't have NAs or empty results in our proportions
    if (length(percentzero) == 0 || is.na(percentzero)) percentzero <- 0
    if (length(perczeroAbb) == 0 || is.na(perczeroAbb)) perczeroAbb <- 0

    if (!isTRUE(percentzero > 0.01 || perczeroAbb > 0.1)) {
      if (goodness_of_fit > 0.40) { # Only log high-goodness rejections (decimal scale)
        log_debug("Rejecting high-goodness candidate (no losses): rho={round(rho,3)}, psi={round(psi,3)}, goodness={round(goodness_of_fit,2)}, pz={round(percentzero,4)}, pza={round(perczeroAbb,4)}")
      }
      return(list(valid = FALSE, reason = "zero_constraint", goodness = goodness_of_fit))
    } else {
      if (goodness_of_fit > 0.40) {
        log_debug("Accepting candidate: rho={round(rho,3)}, psi={round(psi,3)}, goodness={round(goodness_of_fit,2)}, pz={round(percentzero,4)}, pza={round(perczeroAbb,4)}")
      }
    }
  }

  return(list(valid = TRUE, psi = psi, rho = min(rho, 1.0), ploidy = ploidy, goodness = goodness_of_fit))
}

#' robust make_segments with tolerance
#' @noRd
make_segments_internal <- function(r, b) {
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
