#' Prioritize candidate integer copy number states around a fractional state
#'
#' Returns candidate grid edges based on BAF and LogR position, following Battenberg's
#' original prioritization rules (LogR distance + simplicity).
#'
#' @param full logical; if TRUE return all 6 candidate edges (for subclonal search),
#'                   if FALSE return only the best edge (for clonal likelihood).
#' @return A list containing matrices `nMaj1`, `nMin1`, `nMaj2`, `nMin2` (NxM),
#'         and `nMaj`, `nMin` (Nx2) for the best edge corners.
#' @noRd
prioritizeCopyNumbers <- function(rho, psi, BAF_req, nMajor, nMinor, full = TRUE) {
  # Vectorized Inputs
  x <- floor(nMinor)
  y <- floor(nMajor)
  ntot <- nMajor + nMinor
  n <- length(BAF_req)

  # Pre-calculate BAF at key corners
  calc_baf <- function(nM, nm) {
    num <- 1 - rho + rho * nM
    den <- 2 - 2 * rho + rho * (nM + nm)
    lev <- num / den
    lev[nM == 0 & nm == 0] <- 0.5
    lev
  }

  lev3 <- calc_baf(y, x) # Corner C3
  lev2 <- calc_baf(y + 1, x + 1) # Corner C2

  case_1_2a <- BAF_req > lev3
  case_2c <- (!case_1_2a) & (BAF_req > lev2)
  logR_low <- ntot < (x + y + 1)

  # Initialize matrices for all 6 possible candidates
  # We use the offsets defined in original orderEdges logic
  m1 <- matrix(0, n, 6)
  n1 <- matrix(0, n, 6)
  m2 <- matrix(0, n, 6)
  n2 <- matrix(0, n, 6)

  # Helper to fill offsets for a logical mask
  fill_offsets <- function(mask, om1, on1, om2, on2) {
    # Guard against NAs in mask
    mask[is.na(mask)] <- FALSE
    if (any(mask)) {
      m1[mask, ] <<- sweep(matrix(om1, sum(mask), 6, byrow = TRUE), 1, y[mask], "+")
      n1[mask, ] <<- sweep(matrix(on1, sum(mask), 6, byrow = TRUE), 1, x[mask], "+")
      m2[mask, ] <<- sweep(matrix(om2, sum(mask), 6, byrow = TRUE), 1, y[mask], "+")
      n2[mask, ] <<- sweep(matrix(on2, sum(mask), 6, byrow = TRUE), 1, x[mask], "+")
    }
  }

  # Fill based on original Battenberg orderEdges logic
  fill_offsets(
    case_1_2a & logR_low,
    c(0, -1, 0, 1, 1, 1), c(0, 0, 0, 0, -1, 0),
    c(1, 1, 2, 1, 1, 1), c(0, 0, 0, 1, 1, 2)
  )
  fill_offsets(
    case_1_2a & (!logR_low),
    c(1, 1, 1, 0, -1, 0), c(0, -1, 0, 0, 0, 0),
    c(1, 1, 1, 1, 1, 2), c(1, 1, 2, 0, 0, 0)
  )
  fill_offsets(
    case_2c & logR_low,
    c(0, 0, 0, 1, 1, 1), c(0, -1, 0, 0, -1, 0),
    c(0, 0, 0, 1, 1, 1), c(1, 1, 2, 1, 1, 2)
  )
  fill_offsets(
    case_2c & (!logR_low),
    c(1, 1, 1, 0, 0, 0), c(0, -1, 0, 0, -1, 0),
    c(1, 1, 1, 0, 0, 0), c(1, 1, 2, 1, 1, 2)
  )
  fill_offsets(
    (!case_1_2a) & (!case_2c) & logR_low,
    c(0, 0, 0, 0, -1, 0), c(0, -1, 0, 1, 1, 1),
    c(0, 0, 0, 1, 1, 2), c(1, 1, 2, 1, 1, 1)
  )
  fill_offsets(
    (!case_1_2a) & (!case_2c) & (!logR_low),
    c(0, -1, 0, 0, 0, 0), c(1, 1, 1, 0, -1, 0),
    c(0, 0, 0, 0, 0, 0), c(1, 1, 2, 1, 1, 2)
  )

  # Validation: Avoid negative CNs
  invalid <- (m1 < 0 | n1 < 0 | m2 < 0 | n2 < 0)
  invalid[is.na(invalid)] <- TRUE # Treat NAs as invalid
  m1[invalid] <- NA
  n1[invalid] <- NA
  m2[invalid] <- NA
  n2[invalid] <- NA

  if (full) {
    return(list(
      nMaj1 = m1, nMin1 = n1, nMaj2 = m2, nMin2 = n2,
      nMaj = cbind(m1[, 1], m2[, 1]), nMin = cbind(n1[, 1], n2[, 1])
    ))
  } else {
    return(list(
      nMaj1 = m1[, 1, drop = FALSE], nMin1 = n1[, 1, drop = FALSE],
      nMaj2 = m2[, 1, drop = FALSE], nMin2 = n2[, 1, drop = FALSE],
      nMaj = cbind(m1[, 1], m2[, 1]), nMin = cbind(n1[, 1], n2[, 1])
    ))
  }
}
