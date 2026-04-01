#' Run code in parallel or serial based on debug status
#'
#' A helper function to abstract the pattern of switching between parallel
#' execution via foreach and serial execution via lapply.
#'
#' @param iterator A vector or list to iterate over (e.g., seq_along(x)).
#' @param func A function to apply to each element of the iterator.
#' @param libs Path to library paths for workers.
#'
#' @return A list of results from the applied function.
#' @keywords internal
run_with_error_handling <- function(iterator, func, libs, nthreads = 1) {
  if (length(iterator) == 0) {
    return(list())
  }

  # Set up foreach to use the registered backend
  # Use %dopar% if a backend is registered and nthreads > 1, else %do%
  `%op%` <- if (foreach::getDoParWorkers() > 1) foreach::`%dopar%` else foreach::`%do%`

  results <- foreach::foreach(i = iterator) %op% {
    # Set thread budget for this worker
    data.table::setDTthreads(nthreads)
    Sys.setenv(OMP_NUM_THREADS = nthreads, MKL_NUM_THREADS = nthreads, OPENBLAS_NUM_THREADS = nthreads)

    .libPaths(libs)

    # Execute the function and capture its result
    worker_result <- withCallingHandlers(
      {
        func(i)
      },
      error = function(e) {
        msg <- sprintf("!!! BATTENBERG ERROR IN PARALLEL WORKER NODE %s !!!\nMessage: %s\nStack Trace:", i, conditionMessage(e))
        calls <- sys.calls()
        for (j in rev(seq_along(calls))) {
          msg <- paste(msg, sprintf("%d: %s", j, deparse(calls[[j]])), sep = "\n")
        }
        msg <- paste(msg, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", sep = "\n")
        stop(msg, call. = FALSE)
      }
    )

    # Trigger garbage collection after each worker finishes its task to free up RAM
    gc()

    # The last expression in the loop body is what gets returned to the results list
    worker_result
  }
  return(results)
}

#' Safe wrapper for mclapply that prevents deadlocks
#'
#' This function disables data.table multi-threading before forking and restores it after.
#' This is critical to prevent hangs in Singularity/Linux environments.
#'
#' @param X A vector or list to iterate over.
#' @param FUN The function to be applied.
#' @param mc.cores The number of cores to use.
#' @param ... Additional arguments passed to mclapply.
#' @return A list of results.
bt_mclapply <- function(X, FUN, mc.cores = 1, ...) {
  # If we are already in a parallel worker (e.g., from sample-level parallelism),
  # or if 1 core is requested, we MUST run sequentially.
  is_nested <- FALSE
  if (requireNamespace("foreach", quietly = TRUE)) {
    is_nested <- foreach::getDoParWorkers() > 1
  }

  if (mc.cores <= 1 || is_nested) {
    return(lapply(X, FUN, ...))
  }

  # Ensure data.table multi-threading is off before forking to prevent deadlocks
  old_threads <- data.table::getDTthreads()
  data.table::setDTthreads(1)

  on.exit({
    data.table::setDTthreads(old_threads)
  })

  # mc.preschedule=FALSE is more stable in container environments with varying task sizes
  parallel::mclapply(X, FUN, mc.cores = mc.cores, mc.preschedule = FALSE, ...)
}

#' Helper for chromosome-aware smart downsampling for plot performance
#'
#' @param v The vector to downsample.
#' @param target Target number of points.
#' @return A vector of indices to keep.
bt_downsample_indices <- function(v, target) {
  n <- length(v)
  if (n <= target) {
    return(seq_along(v))
  }
  # We use a combined approach: uniform sampling + local extremes (min/max)
  # to preserve visual dips/peaks in LogR/BAF
  bin_size <- ceiling(n / (target / 2))

  # Use data.table for speed and concise grouping
  dt_ds <- data.table::data.table(val = as.numeric(v), id = seq_along(v))
  dt_ds[, bin := ceiling(id / bin_size)]
  keep <- dt_ds[, .(id_min = id[which.min(val)], id_max = id[which.max(val)]), by = bin]

  sort(unique(c(keep$id_min, keep$id_max)))
}
