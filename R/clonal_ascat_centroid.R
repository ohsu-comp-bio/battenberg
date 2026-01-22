####################################################################################################
#' This function is an alternative procedure for finding the optimum (psi, rho) pair.
#' This function first finds all the find all the global optima,
#' and then finds the centroid of this set of globla optima.
#' Then we find the global optimum which is nearest to the centroid.
#' (When the set of global optima is convex, we expect the selected optimum to be at the centroid.)
#' @param d A distance matrix
#' @param ref_seg_matrix The corresponding ref seg matrix that belongs to d
#' @param ref_major The corresponding major allele values with d
#' @param ref_minor The corresponding minor allele values with d
#' @param s A segmented BAF/LogR data.frame from \code{get_segment_info}
#' @param dist_choice Some distance metrics require adaptation of the data (i.e. log transform)
#' @param minimise Boolean whether we're minimising or maximising
#' @param new_bounds The rho/psi boundaries between we are searching for a solution. This is a named list with values psi_min, psi_max, rho_min, rho_max
#' @param distancepng String where the sunrise distance plot will be saved
#' @param gamma_param The platform gamma
#' @param siglevel_BAF The level at which BAF becomes significant TODO: this option is no longer used
#' @param maxdist_BAF TODO: this option is no longer used
#' @param siglevel_LogR The p-value at which logR becomes significant when establishing whether a segment should be subclonal
#' @param maxdist_LogR The maximum distance allowed as slack when establishing the significance. This allows for the case when a breakpoint is missed, the segment would then not automatically become subclonal
#' @param allow100percent Boolean whether to allow for a 100"\%" cellularity solution
#' @param uninformative_baf_threshold The threshold above which BAF becomes uninformative
#' @param read_depth TODO: this option is no longer used
#' @return A list with fields optima_info_without_ref and optima_info
#' @export
find_centroid_of_global_minima <- function(
  d, ref_seg_matrix, ref_major, ref_minor, s, dist_choice, minimise,
  new_bounds, distancepng, gamma_param, siglevel_BAF, maxdist_BAF,
  siglevel_LogR, maxdist_LogR, allow100percent, uninformative_baf_threshold,
  read_depth
) {
  if (!minimise) {
    d <- -d # This ensures that we "maximise" instead of "minimise"!
  }

  # Find height of global minima
  gmin <- min(d, na.rm = TRUE)

  # Find all global minima
  nropt <- 0
  optima <- list()

  # Pre-extract psi/rho values from grid
  psi_grid <- as.numeric(rownames(d))
  rho_grid <- as.numeric(colnames(d))

  for (i in 1:nrow(d)) {
    for (j in 1:ncol(d)) {
      if (!is.na(d[i, j]) && d[i, j] == gmin) {
        psi <- psi_grid[i]
        rho <- rho_grid[j]

        # Calculate ploidy
        term_base <- (rho - 1)
        term_psi <- ((1 - rho) * 2 + rho * psi)
        factor <- 2^(s[, "r"] / gamma_param)

        nA <- (term_base - (s[, "b"] - 1) * factor * term_psi) / rho
        nB <- (term_base + s[, "b"] * factor * term_psi) / rho

        ploidy <- sum((nA + nB) * s[, "length"], na.rm = TRUE) / sum(s[, "length"])

        # goodnessOfFit is the same as gmin in this implementation
        goodnessOfFit <- gmin

        nropt <- nropt + 1
        optima[[nropt]] <- list(gmin = gmin, i = i, j = j, ploidy = ploidy, gof = goodnessOfFit)
      }
    }
  }

  # Find a "centroid" of the set of global minima
  grid_x_vect <- sapply(optima, function(z) z$i)
  grid_y_vect <- sapply(optima, function(z) z$j)

  centre_x <- median(grid_x_vect)
  centre_y <- median(grid_y_vect)
  centre <- c(centre_x, centre_y)

  index <- 1
  sqrdist_min <- Inf
  for (i in 1:length(optima)) {
    grid_point <- c(optima[[i]]$i, optima[[i]]$j)
    sqrdist <- (grid_point[1] - centre[1])^2 + (grid_point[2] - centre[2])^2

    if (sqrdist <= sqrdist_min) {
      sqrdist_min <- sqrdist
      index <- i
    }
  }

  grid_x <- optima[[index]]$i
  grid_y <- optima[[index]]$j

  psi_opt1 <- psi_grid[grid_x]
  rho_opt1 <- min(rho_grid[grid_y], 1)
  ploidy_opt1 <- optima[[index]]$ploidy
  goodnessOfFit_opt1 <- optima[[index]]$gof

  ref_seg <- ref_seg_matrix[grid_x, grid_y]

  if (minimise) {
    dist_optima <- gmin
  } else {
    dist_optima <- -gmin
    goodnessOfFit_opt1 <- -goodnessOfFit_opt1
  }

  # First optima set (without reference segment override)
  optima_info_without_ref <- list(
    nropt = nropt, psi_opt1 = psi_opt1, rho_opt1 = rho_opt1,
    ploidy_opt1 = ploidy_opt1, ref_seg = ref_seg,
    goodnessOfFit_opt1 = goodnessOfFit_opt1
  )

  # Logic for reference segment override
  if (ref_seg == 0) {
    psi_opt1 <- 2
    rho_opt1 <- 1
    ploidy_opt1 <- 2
    goodnessOfFit_opt1 <- 1
  } else {
    ref_segment_info <- get_psi_rho_from_ref_seg(
      ref_seg, s, ref_major[grid_x, grid_y], ref_minor[grid_x, grid_y], gamma_param
    )

    psi_opt1 <- ref_segment_info$psi
    rho_opt1 <- ref_segment_info$rho
    ploidy_opt1 <- ref_segment_info$ploidy

    if (!is.na(rho_opt1)) {
      distance_info <- calc_distance_clonal(
        s, dist_choice, rho_opt1, psi_opt1, gamma_param, read_depth,
        siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR, uninformative_baf_threshold
      )
      goodnessOfFit_opt1 <- distance_info$distance_value
    } else {
      goodnessOfFit_opt1 <- Inf
    }
  }

  # Final optima set
  optima_info <- list(
    nropt = nropt, psi_opt1 = psi_opt1, rho_opt1 = rho_opt1,
    ploidy_opt1 = ploidy_opt1, ref_seg = ref_seg,
    goodnessOfFit_opt1 = goodnessOfFit_opt1
  )

  # Plotting
  if (!is.na(distancepng)) {
    rhos <- c(optima_info_without_ref$rho_opt1, rho_opt1)
    psis <- c(optima_info_without_ref$psi_opt1, psi_opt1)

    grDevices::png(filename = distancepng, width = 1000, height = 1000, res = 1000 / 7, type = "cairo")
    clonal_findcentroid_plot(minimise, dist_choice, -d, psis, rhos, new_bounds)
    grDevices::dev.off()
  }

  return(list(optima_info_without_ref = optima_info_without_ref, optima_info = optima_info))
}
