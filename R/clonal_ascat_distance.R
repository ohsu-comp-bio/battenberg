#' This function computes various "distances", which are used as penalties for a copy number solution.
#' This function is called when searching for a clonal copy number solution.
#' @noRd
calc_distance <- function(segs, dist_choice, rho, psi, gamma_param, uninformative_baf_threshold = 0.51) {
  s <- segs

  # common nA/nB logic
  mult <- 2^(s[, "r"] / gamma_param) * ((1 - rho) * 2 + rho * psi)
  nA <- (rho - 1 - (s[, "b"] - 1) * mult) / rho
  nB <- (rho - 1 + s[, "b"] * mult) / rho

  # CRITICAL FIX #6: Clamp negative copy numbers to 0.01 (from original battenberg/R/clonal_ascat.R:474-480)
  # At low rho values (e.g., 0.07), the formulas can produce negative copy numbers (e.g., nB=-5.37)
  # which create extreme distances ~10x larger than theoretical max, causing optimization to fail
  # The original code clamps these to 0.01 to prevent this mathematical breakdown
  nA[nA < 0 | is.na(nA)] <- 0.01
  nB[nB < 0 | is.na(nB)] <- 0.01

  if (dist_choice == 0) { # original ASCAT distance
    sum_nA <- sum(nA, na.rm = TRUE)
    sum_nB <- sum(nB, na.rm = TRUE)
    if (sum_nA < sum_nB) {
      nMinor <- nA
    } else {
      nMinor <- nB
    }
    dist_value <- sum(abs(nMinor - pmax(round(nMinor), 0))^2 * s[, "length"] * ifelse(s[, "b"] <= uninformative_baf_threshold, 0.05, 1), na.rm = TRUE)
    minimise <- TRUE
  } else if (dist_choice == 1) { # new similarity measure suggested by DW 7-3-2014
    sum_nA <- sum(nA, na.rm = TRUE)
    sum_nB <- sum(nB, na.rm = TRUE)
    if (sum_nA < sum_nB) {
      nMinor <- nA
    } else {
      nMinor <- nB
    }
    dist_value <- sum((0.5 - abs(nMinor - pmax(round(nMinor), 0)))^2 * s[, "length"], na.rm = TRUE)
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
    dist_value <- 0.5 * sum(((0.5 - abs(nMinor - pmax(round(nMinor), 0)))^2 + (0.5 - abs(nMajor - pmax(round(nMajor), 0)))^2) * s[, "length"], na.rm = TRUE)
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
    segs_penalty <- (0.5 - abs(nMinor - pmax(round(nMinor), 0)))^2 + (0.5 - abs(nMajor - pmax(round(nMajor), 0)))^2
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
      i <- grid$psi_idx[idx]
      j <- grid$rho_idx[idx]
      distance_info <- calc_distance(s, dist_choice, rho_pos[j], psi_pos[i], gamma_param, uninformative_baf_threshold = uninformative_baf_threshold)
      return(list(i = i, j = j, val = distance_info$distance_value, minimise = distance_info$minimise))
    }, mc.cores = nthreads)

    for (res in results) {
      d[res$i, res$j] <- res$val
    }
    minimise <- results[[1]]$minimise
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

  run_grid_point <- function(idx) {
    psi <- grid$psi[idx]
    rho <- grid$rho[idx]

    res <- calc_distance_clonal(
      s, dist_choice, rho, psi, gamma_param, read_depth,
      siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR,
      uninformative_baf_threshold
    )
    return(res)
  }

  if (nthreads > 1 && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(seq_len(nrow(grid)), run_grid_point, mc.cores = nthreads)
  } else {
    results <- lapply(seq_len(nrow(grid)), run_grid_point)
  }

  # Extract values
  # Handle both list and atomic vector results from mclapply
  d_vals <- sapply(results, function(x) if (is.list(x)) x$distance_value else NA)
  ref_vals <- sapply(results, function(x) if (is.list(x)) x$max_clonal_segment else NA)
  maj_vals <- sapply(results, function(x) if (is.list(x)) x$ref_maj else NA)
  min_vals <- sapply(results, function(x) if (is.list(x)) x$ref_min else NA)

  dist_mat <- matrix(d_vals, nrow = length(psi_pos), ncol = length(rho_pos))
  ref_seg_mat <- matrix(ref_vals, nrow = length(psi_pos), ncol = length(rho_pos))
  ref_major_mat <- matrix(maj_vals, nrow = length(psi_pos), ncol = length(rho_pos))
  ref_minor_mat <- matrix(min_vals, nrow = length(psi_pos), ncol = length(rho_pos))

  rownames(dist_mat) <- psi_pos
  colnames(dist_mat) <- rho_pos

  return(list(
    distance_matrix = dist_mat,
    minimise = results[[1]]$minimise,
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
  uninformative_baf_threshold
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
  max_clonal_segment_size <- 0
  ref_maj <- NA
  ref_min <- NA

  # Filter informative segments
  informative_idx <- which(s[, "b"] > uninformative_baf_threshold)

  if (length(informative_idx) == 0) {
    # If no segments informative, we still need to return a structure
    return(list(distance_value = Inf, minimise = TRUE, max_clonal_segment = 0, ref_maj = NA, ref_min = NA))
  }

  for (i in informative_idx) {
    BAFreq <- s[i, "b"]
    LogR <- s[i, "r"]
    BAF_length <- s[i, "length"]
    BAF_size <- s[i, "size"]
    BAF_mean <- s[i, "mean"]
    BAF_sd <- s[i, "sd"]

    # Calculate Clonal Status
    segment_info <- is_segment_clonal(
      LogR = LogR, BAF_req = BAFreq, BAF_length = BAF_length,
      BAF_size = BAF_size, BAF_mean = BAF_mean, BAF_sd = BAF_sd,
      read_depth = read_depth, rho = rho, psi = psi, gamma_param = gamma_param,
      siglevel_BAF = siglevel_BAF, maxdist_BAF = maxdist_BAF,
      siglevel_LogR = siglevel_LogR, maxdist_LogR = maxdist_LogR
    )

    is_clonal <- segment_info$is_clonal
    nMaj <- segment_info$nMaj
    nMin <- segment_info$nMin
    is_balanced <- segment_info$balanced

    segment_size <- BAF_length
    genome_size <- genome_size + segment_size
    seg_count <- seg_count + 1

    if (is_clonal) {
      clonal_genome_size <- clonal_genome_size + segment_size
      if (max_clonal_segment_size < segment_size && !is_balanced) {
        max_clonal_segment <- i
        max_clonal_segment_size <- segment_size
        ref_maj <- nMaj
        ref_min <- nMin
      }
    }

    # Calculate Standard Error
    standard_error_info <- calc_standardised_error(
      LogR, BAFreq, BAF_length, BAF_size, BAF_mean, BAF_sd,
      rho, psi, gamma_param, maxdist_BAF
    )

    if (standard_error_info$included_segment > 0) {
      n_included_segments <- n_included_segments + 1
      sum1 <- sum1 + standard_error_info$tvar^2
    }

    # These sums follow the original iterative logic
    sum2 <- sum2 + (BAFreq - BAF_mean)^2
    sum3 <- sum3 + (segment_size * (BAFreq - BAF_mean)^2)

    # Calculate Log Likelihood Ratio
    ln_lratio <- calc_ln_likelihood_ratio(
      LogR, BAFreq, BAF_length, BAF_size, BAF_mean,
      read_depth, rho, psi, gamma_param, maxdist_BAF
    )
    sum_ln_lratio <- sum_ln_lratio + ln_lratio
  }

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
