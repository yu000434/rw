#' Predictive mean matching with donor row tracking
#'
#' This imputation method follows the same basic PMM donor-selection step as
#' `mice::mice.impute.pmm()`, but stores the observed donor row selected for
#' each imputed recipient. Use it inside `mice()` with `method = "pmmrw"` and
#' `tasks = "train"` when PMM donor-source correlation will be used by
#' `pool_rw()`.
#'
#' @inheritParams mice::mice.impute.pmm
#' @details
#' The current implementation is intentionally limited to numeric PMM with
#' donor row tracking. The recorded donor IDs let downstream variance
#' calculations group analysis-score contributions by observed donor source.
#' It requires `tasks = "train"` and currently does not support `exclude`,
#' `use.matcher = TRUE`, or `mlocal` values other than 1.
#' @return A vector of imputed values.
#' @export
mice.impute.pmmrw <- function(y, ry, x, wy = NULL, task = "impute", model = NULL,
                              exclude = NULL, ridge = 1e-05, matchtype = 1L,
                              donors = 5L, nbins = NULL, use.matcher = FALSE,
                              mlocal = 1L, ...) {
  if (!is.null(exclude)) {
    stop("`pmmrw` does not currently support `exclude`.", call. = FALSE)
  }
  if (is.factor(y)) {
    stop("`pmmrw` currently supports numeric variables only.", call. = FALSE)
  }
  if (isTRUE(use.matcher)) {
    stop("`pmmrw` currently supports `use.matcher = FALSE` only.", call. = FALSE)
  }
  if (!identical(as.integer(mlocal), 1L)) {
    stop("`pmmrw` currently supports `mlocal = 1` only.", call. = FALSE)
  }
  if (is.null(wy)) {
    wy <- !ry
  }
  if (task != "train") {
    stop("`pmmrw` must be used with `tasks = 'train'`.", call. = FALSE)
  }
  if (is.null(model) || !is.environment(model)) {
    stop("`model` must be an environment; use `tasks = 'train'`.", call. = FALSE)
  }

  row_id <- seq_along(y)

  x <- cbind(1, as.matrix(x))
  parm <- mice:::.norm.draw(y, ry, x, ridge = ridge, ...)
  beta_hat <- drop(parm$coef)
  beta_dot <- drop(parm$beta)

  if (matchtype == 0L) {
    beta_dot <- beta_hat
  } else if (matchtype == 2L) {
    beta_hat <- beta_dot
  }

  yhat_obs <- as.vector(x[ry, , drop = FALSE] %*% beta_hat)
  yhat_mis <- as.vector(x[wy, , drop = FALSE] %*% beta_dot)

  nbins <- initialize_nbins_pmmrw(nbins, length(yhat_obs), length(unique(yhat_obs)))
  donors <- initialize_donors_pmmrw(donors, length(yhat_obs))
  edges <- stats::quantile(
    yhat_obs,
    probs = seq(0, 1, length.out = nbins + 1L),
    type = 7L,
    na.rm = TRUE
  )
  lookup <- bin_yhat_donor_pmmrw(
    yhat = yhat_obs,
    y = y[ry],
    row_id = row_id[ry],
    k = donors,
    edges = edges
  )
  draws <- draw_neighbors_pmmrw(
    yhat = yhat_mis,
    edges = edges,
    lookup = lookup$value,
    donor_lookup = lookup$donor_id,
    mlocal = mlocal
  )

  donor_row <- rep(NA_integer_, length(y))
  donor_row[wy] <- draws$donor_id[, 1L]

  model$setup <- list(
    method = "pmmrw",
    n = length(yhat_obs),
    task = task,
    donors = donors,
    matchtype = matchtype,
    nbins = nbins,
    ridge = ridge
  )
  model$beta.hat <- beta_hat
  model$beta.dot <- beta_dot
  model$edges <- edges
  model$lookup <- lookup$value
  model$donor_lookup <- lookup$donor_id
  model$sigma.dot <- parm$sigma
  model$xnames <- colnames(x)
  model$donor_id <- donor_row

  draws$value[, 1L]
}

initialize_nbins_pmmrw <- function(nbins, n, nu) {
  if (is.null(nbins)) {
    nbins <- round(4 * log(n) + 1.5)
  }
  if (nbins > nu) {
    nbins <- nu
  }
  max(2L, nbins)
}

initialize_donors_pmmrw <- function(donors, n) {
  if (is.null(donors)) {
    donors <- round(n / 600 + 7)
  }
  max(1L, min(donors, n))
}

bin_yhat_donor_pmmrw <- function(yhat, y, row_id, k, edges) {
  stopifnot(length(yhat) == length(y), length(yhat) == length(row_id))

  sort_order <- order(yhat)
  yhat_sorted <- yhat[sort_order]
  y_sorted <- y[sort_order]
  row_sorted <- row_id[sort_order]

  bin <- findInterval(yhat_sorted, vec = edges, all.inside = TRUE)
  row_index <- seq_along(y_sorted)
  index_list <- split(row_index, bin)
  nbins <- length(edges) - 1L

  lookup_index <- t(vapply(seq_len(nbins), function(b) {
    index <- index_list[[as.character(b)]]
    if (length(index) == 0L) {
      sample(row_index, size = k, replace = TRUE)
    } else if (length(index) == 1L) {
      rep(index, k)
    } else {
      sample(index, size = k, replace = length(index) < k)
    }
  }, integer(k)))

  list(
    value = matrix(y_sorted[lookup_index], nrow = nbins),
    donor_id = matrix(row_sorted[lookup_index], nrow = nbins)
  )
}

draw_neighbors_pmmrw <- function(yhat, edges, lookup, donor_lookup, mlocal = 1L) {
  n <- length(yhat)
  nbins <- length(edges) - 1L

  bin <- findInterval(yhat, edges, rightmost.closed = TRUE, all.inside = TRUE)
  t0 <- edges[pmax(bin, 1L)]
  t1 <- edges[pmin(bin + 1L, nbins)]
  p_left <- ifelse(t1 > t0, (t1 - yhat) / (t1 - t0), 0.5)

  selected_bin <- ifelse(stats::runif(n) < p_left, bin, pmin(bin + 1L, nbins))
  indices <- matrix(sample(1L:ncol(lookup), n * mlocal, replace = TRUE), nrow = n)

  list(
    value = matrix(lookup[cbind(selected_bin, indices)], nrow = n, ncol = mlocal),
    donor_id = matrix(donor_lookup[cbind(selected_bin, indices)], nrow = n, ncol = mlocal)
  )
}

#' Extract PMM donor row IDs
#'
#' @param object A `mids` object created by `mice()` using `method = "pmmrw"`
#'   and `tasks = "train"`.
#' @param vars Optional PMM variable names. If omitted, all `pmmrw` variables
#'   are used.
#' @return For one PMM variable, an `n` by `m` integer matrix. Entries are `NA`
#'   for rows that were not imputed by PMM and the observed donor row number for
#'   PMM recipients. For multiple variables, a named list of such matrices.
#' @export
extract_donor_id <- function(object, vars = NULL) {
  if (!inherits(object, "mids")) {
    cli::cli_abort("{.arg object} must be a {.cls mids} object.")
  }
  if (is.null(object$models)) {
    cli::cli_abort("{.arg object} must contain models from {.code tasks = 'train'}.")
  }

  if (is.null(vars)) {
    vars <- names(object$method)[object$method == "pmmrw"]
  }
  if (length(vars) == 0L) {
    cli::cli_abort("No {.val pmmrw} variables were found.")
  }

  out <- lapply(vars, function(var) {
    if (!var %in% names(object$models)) {
      cli::cli_abort("No stored model found for {.field {var}}.")
    }
    donor_list <- lapply(seq_len(object$m), function(i) {
      donor <- object$models[[var]][[i]]$donor_id
      if (is.null(donor)) {
        cli::cli_abort(
          "No donor IDs were recorded for {.field {var}} imputation {i}."
        )
      }
      as.integer(donor)
    })
    mat <- do.call(cbind, donor_list)
    colnames(mat) <- paste0("imp", seq_len(object$m))
    mat
  })
  names(out) <- vars

  if (length(out) == 1L) {
    out[[1L]]
  } else {
    out
  }
}
