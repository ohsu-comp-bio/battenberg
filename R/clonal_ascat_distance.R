#' This function computes various "distances", which are used as penalties for a copy number solution.
#' This function is called when searching for a clonal copy number solution.
#' @noRd
calc_distance <- function(segs, dist_choice, rho, psi, gamma_param, uninformative_baf_threshold = 0.51) {
  s <- segs

  # common nA/nB logic
  mult <- 2^(s[, "r"] / gamma_param) * ((1 - rho) * 2 + rho * psi)
  nA <- (rho - 1 - (s[, "b"] - 1) * mult) / rho
  nB <- (rho - 1 + s[, "b"] * mult) / rho


  if (dist_choice == 0) { # original ASCAT distance
    sum_nA <- sum(nA, na.rm = TRUE)
    sum_nB <- sum(nB, na.rm = TRUE)
    if (sum_nA < sum_nB) {
      nMinor <- nA
    } else {
      nMinor <- nB
    }

    # Correctly identify uninformative BAF (near 0.5)
    # Original logic using <= threshold is dangerous for unmirrored BAF (0..1)
    # We want to downweight ONLY values close to 0.5
    # Fallback to a tight window (0.49-0.51) if threshold is weird, or just trust the threshold logic
    # Assuming uninformative_baf_threshold is e.g. 0.51 (meaning deviations < 0.01 from 0.5 are noisy)
    # Let's use a robust check: uninformative if distance to 0.5 is small
    # Correctly identify uninformative BAF (near 0.5)
    # Original logic: weight <- ifelse(s[, "b"] <= uninformative_baf_threshold, 0.05, 1)
    # NOTE: Since inputs are Minor Allele (<0.5) and threshold is 0.51, this effectively weights ALL segments as 0.05.
    # While potentially counter-intuitive, this matches the Original Battenberg behavior exactly.
    weight <- ifelse(s[, "b"] <= uninformative_baf_threshold, 0.05, 1)

    dist_value <- sum(abs(nMinor - pmax(round(nMinor), 0))^2 * s[, "length"] * weight, na.rm = TRUE)
    minimise <- TRUE
  } else if (dist_choice == 1) { # new similarity measure suggested by DW 7-3-2014
    sum_nA <- sum(nA, na.rm = TRUE)
    sum_nB <- sum(nB, na.rm = TRUE)
    if (sum_nA < sum_nB) {
      nMinor <- nA
    } else {
      nMinor <- nB
    }
    dist_value <- sum((pmax(0, 0.5 - abs(nMinor - pmax(round(nMinor), 0))))^2 * s[, "length"], na.rm = TRUE)
    minimise <- FALSE
  } else if (dist_choice == 2) { # adapted DW's 7-3-2014 measure by SD 8-8-2014
    sum_nA <- sum(nA, na.rm = TRUE)
    sum_nB <- sum(nB, na.rm = TRUE)
    if (sum_nA < sum_nB) {
      nMinor <- nA
      nMajor <- nB
    } else {
      nMinor <- nB
      nMajor <- nA
    }

    dist_value <- 0.5 * sum((pmax(0, 0.5 - abs(nMinor - pmax(round(nMinor), 0)))^2 + (pmax(0, 0.5 - abs(nMajor - pmax(round(nMajor), 0)))^2)) * s[, "length"], na.rm = TRUE)
    minimise <- FALSE
  } else if (dist_choice == 3) { # adapted DW's 7-3-2014 measure by SD 8-8-2014 with homozygous deletion penalty
    sum_nA <- sum(nA, na.rm = TRUE)
    sum_nB <- sum(nB, na.rm = TRUE)
    if (sum_nA < sum_nB) {
      nMinor <- nA
      nMajor <- nB
    } else {
      nMinor <- nB
      nMajor <- nA
    }

    segs_penalty <- (pmax(0, 0.5 - abs(nMinor - pmax(round(nMinor), 0))))^2 + (pmax(0, 0.5 - abs(nMajor - pmax(round(nMajor), 0))))^2
    hom_del <- nMinor < 0.5 & nMajor < 0.5 & nMinor >= 0 & nMajor >= 0
    segs_penalty[which(hom_del)] <- segs_penalty[which(hom_del)] * 4
    dist_value <- 0.5 * sum(segs_penalty * (s[, "length"] * ifelse(hom_del, 2, 1)), na.rm = TRUE)
    minimise <- FALSE
  }

  return(list(distance_value = dist_value, minimise = minimise))
}

#' Internal optimized grid search distance matrix calculator
#' @noRd
create_distance_matrix <- function(s, dist_choice, gamma_param, uninformative_baf_threshold = 0.51,
                                   min_rho = 0.1, max_rho = 1, min_psi = 1, max_psi = 5.4, nthreads = 1) {
  psi_pos <- seq(min_psi, max_psi, 0.05)
  rho_pos <- seq(min_rho, max_rho, 0.01)

  d <- matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) <- psi_pos
  colnames(d) <- rho_pos

  if (nthreads > 1 && .Platform$OS.type != "windows") {
    grid <- expand.grid(psi_idx = seq_along(psi_pos), rho_idx = seq_along(rho_pos))
    results <- parallel::mclapply(seq_len(nrow(grid)), function(idx) {
      tryCatch(
        {
          i <- grid$psi_idx[idx]
          j <- grid$rho_idx[idx]
          distance_info <- calc_distance(s, dist_choice, rho_pos[j], psi_pos[i], gamma_param, uninformative_baf_threshold = uninformative_baf_threshold)
          return(list(i = i, j = j, val = distance_info$distance_value, minimise = distance_info$minimise))
        },
        error = function(e) {
          return(e)
        }
      )
    }, mc.cores = nthreads)

    valid_results <- results[sapply(results, function(x) is.list(x) && !inherits(x, "error"))]

    if (length(valid_results) < length(results)) {
      warning("Some parallel distance calculations failed.")
    }

    for (res in valid_results) {
      d[res$i, res$j] <- res$val
    }

    if (length(valid_results) > 0) {
      minimise <- valid_results[[1]]$minimise
    } else {
      # Fallback if all failed or empty grid (unlikely)
      # Calculate once synchronously to determine minimise or catch error
      minimise <- TRUE
    }
  } else {
    for (i in seq_along(psi_pos)) {
      psi <- psi_pos[i]
      for (j in seq_along(rho_pos)) {
        rho <- rho_pos[j]
        distance_info <- calc_distance(s, dist_choice, rho, psi, gamma_param, uninformative_baf_threshold = uninformative_baf_threshold)
        d[i, j] <- distance_info$distance_value
      }
    }
    minimise <- distance_info$minimise
  }
  return(list(distance_matrix = d, minimise = minimise))
}

#' Calculate distance matrix for clonal ASCAT
#' @export
create_distance_matrix_clonal <- function(
  s, dist_choice, gamma_param, read_depth, siglevel_BAF, maxdist_BAF,
  siglevel_LogR, maxdist_LogR, uninformative_baf_threshold, new_bounds,
  nthreads = 1
) {
  psi_min <- new_bounds$psi_min
  psi_max <- new_bounds$psi_max
  rho_min <- new_bounds$rho_min
  rho_max <- new_bounds$rho_max

  psi_range <- psi_max - psi_min
  rho_range <- rho_max - rho_min

  delta_psi <- psi_range / 100
  delta_rho <- rho_range / 100

  psi_pos <- seq(psi_min, psi_max, delta_psi)
  rho_pos <- seq(rho_min, rho_max, delta_rho)

  grid <- expand.grid(psi = psi_pos, rho = rho_pos)

  # Pre-calculate informative segments once for the entire grid search
  lenient_threshold <- pmin(uninformative_baf_threshold, 0.505)
  informative_idx <- which(!is.na(s[, "b"]) & pmax(s[, "b"], 1 - s[, "b"]) > lenient_threshold)
  log_debug("Clonal distance check: {length(informative_idx)}/{nrow(s)} segments informative at >{lenient_threshold} threshold")

  run_grid_point <- function(idx) {
    psi <- grid$psi[idx]
    rho <- grid$rho[idx]

    res <- calc_distance_clonal(
      s, dist_choice, rho, psi, gamma_param, read_depth,
      siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR,
      uninformative_baf_threshold,
      informative_idx = informative_idx
    )
    return(res)
  }

  if (nthreads > 1 && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(seq_len(nrow(grid)), run_grid_point, mc.cores = nthreads)
  } else {
    results <- lapply(seq_len(nrow(grid)), run_grid_point)
  }

  # Extract values
  # Handle both list and atomic vector results from mclapply (e.g. error strings)
  get_val <- function(res, field) {
    if (is.list(res) && field %in% names(res)) res[[field]] else NA
  }

  d_vals <- sapply(results, get_val, "distance_value")
  ref_vals <- sapply(results, get_val, "max_clonal_segment")
  maj_vals <- sapply(results, get_val, "ref_maj")
  min_vals <- sapply(results, get_val, "ref_min")

  dist_mat <- matrix(d_vals, nrow = length(psi_pos), ncol = length(rho_pos))
  ref_seg_mat <- matrix(ref_vals, nrow = length(psi_pos), ncol = length(rho_pos))
  ref_major_mat <- matrix(maj_vals, nrow = length(psi_pos), ncol = length(rho_pos))
  ref_minor_mat <- matrix(min_vals, nrow = length(psi_pos), ncol = length(rho_pos))

  rownames(dist_mat) <- psi_pos
  colnames(dist_mat) <- rho_pos

  # Safety check: If mclapply failed, provide a fallback for 'minimise'
  minimise <- if (length(results) > 0 && is.list(results[[1]])) results[[1]]$minimise else TRUE

  return(list(
    distance_matrix = dist_mat,
    minimise = minimise,
    ref_seg_matrix = ref_seg_mat,
    ref_major = ref_major_mat,
    ref_minor = ref_minor_mat
  ))
}

#' Internal function to calculate distance for a single rho/psi
#' @noRd
calc_distance_clonal <- function(
  s, dist_choice, rho, psi, gamma_param, read_depth,
  siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR,
  uninformative_baf_threshold,
  informative_idx = NULL
) {
  # Initialize accumulators
  genome_size <- 0
  clonal_genome_size <- 0
  seg_count <- 0
  n_included_segments <- 0
  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  sum_ln_lratio <- 0

  max_clonal_segment <- 0
  ref_maj <- NA
  ref_min <- NA

  if (is.null(informative_idx)) {
    lenient_threshold <- pmin(uninformative_baf_threshold, 0.505)
    informative_idx <- which(!is.na(s[, "b"]) & pmax(s[, "b"], 1 - s[, "b"]) > lenient_threshold)
  }

  if (length(informative_idx) == 0) {
    # If no segments informative, we still need to return a structure
    return(list(distance_value = 0, minimise = FALSE, max_clonal_segment = 0, ref_maj = NA, ref_min = NA))
  }

  # Vectorized calculation over informative segments
  subset_s <- s[informative_idx, , drop = FALSE]
  segment_info <- is_segment_clonal(
    LogR = subset_s[, "r"],
    BAF_req = subset_s[, "b"],
    BAF_length = subset_s[, "length"],
    BAF_size = subset_s[, "size"],
    BAF_mean = subset_s[, "mean"],
    BAF_sd = subset_s[, "sd"],
    read_depth = read_depth,
    rho = rho,
    psi = psi,
    gamma_param = gamma_param,
    siglevel_BAF = siglevel_BAF,
    maxdist_BAF = maxdist_BAF,
    siglevel_LogR = siglevel_LogR,
    maxdist_LogR = maxdist_LogR
  )

  is_clonal <- segment_info$is_clonal
  nMaj <- segment_info$nMaj
  nMin <- segment_info$nMin
  is_balanced <- segment_info$balanced

  # Genomic stats (Vectorized)
  segment_sizes <- subset_s[, "length"]
  genome_size <- sum(segment_sizes)
  seg_count <- length(segment_sizes)
  clonal_genome_size <- sum(segment_sizes[is_clonal])

  # Reference segment selection (Vectorized)
  max_clonal_segment <- 0
  ref_maj <- NA
  ref_min <- NA
  is_ref_candidate <- is_clonal & !is_balanced

  if (any(is_ref_candidate)) {
    candidates_sizes <- segment_sizes[is_ref_candidate]
    idx_in_candidates <- which.max(candidates_sizes)
    # Map back to original indices
    max_clonal_segment <- informative_idx[is_ref_candidate][idx_in_candidates]
    ref_maj <- nMaj[is_ref_candidate][idx_in_candidates]
    ref_min <- nMin[is_ref_candidate][idx_in_candidates]
  }

  # Standard error (Batch operation)
  tvars <- calc_batch_standardised_errors(subset_s, rho, psi, gamma_param)
  # Standard errors include segments where size > 0 and sd != 0
  is_valid_se <- subset_s[, "size"] > 0 & !is.na(subset_s[, "sd"]) & subset_s[, "sd"] != 0
  n_included_segments <- sum(is_valid_se)
  sum1 <- sum(tvars[is_valid_se]^2)

  # Distance sums (Vectorized)
  baf_diff_sq <- (subset_s[, "b"] - subset_s[, "mean"])^2
  sum2 <- sum(baf_diff_sq)
  sum3 <- sum(subset_s[, "length"] * baf_diff_sq)

  # Log Likelihood Ratio (Batch operation)
  ln_lratios <- calc_batch_ln_likelihood_ratios(subset_s, read_depth, rho, psi, gamma_param)
  sum_ln_lratio <- sum(ln_lratios)

  # Calculate final distance values
  clonal_proportion <- if (genome_size > 0) clonal_genome_size / genome_size else 0
  dist1 <- if (n_included_segments > 0) sum1 / n_included_segments else 0
  dist2 <- if (seg_count > 0) sum2 / seg_count else 0
  dist3 <- if (genome_size > 0) sum3 / genome_size else 0

  if (dist_choice == 0) {
    dist_value <- clonal_proportion
    minimise <- FALSE
  } else if (dist_choice == 1) {
    dist_value <- dist1
    minimise <- TRUE
  } else if (dist_choice == 2) {
    dist_value <- dist2
    minimise <- TRUE
  } else if (dist_choice == 3) {
    dist_value <- dist3
    minimise <- TRUE
  } else if (dist_choice == 4) {
    dist_value <- sum_ln_lratio
    minimise <- FALSE
  } else {
    # Default fallback
    dist_value <- clonal_proportion
    minimise <- FALSE
  }

  return(list(
    distance_value = dist_value,
    minimise = minimise,
    max_clonal_segment = max_clonal_segment,
    ref_maj = ref_maj,
    ref_min = ref_min
  ))
}
