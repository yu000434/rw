smooth_pmm_working_scale <- function(eta_obs, eta_mis, donors) {
  delta <- outer(eta_obs, eta_mis, "-")
  dist2 <- delta * delta

  bw_scale <- min(stats::sd(eta_obs), stats::IQR(eta_obs) / 1.34,
                  na.rm = TRUE)
  if (!is.finite(bw_scale) || bw_scale <= 0) {
    bw_scale <- stats::sd(eta_obs)
  }
  h_floor <- max(1e-8 * stats::sd(eta_obs), 1e-12)
  h_rule <- max(0.9 * bw_scale * length(eta_obs)^(-1 / 5), h_floor)

  lambda0 <- rep(0.5 / h_rule^2, length(eta_mis))
  log_weight <- -sweep(dist2, 2L, lambda0, "*")
  log_weight <- sweep(log_weight, 2L, apply(log_weight, 2L, max), "-")
  w0 <- exp(log_weight)
  w0 <- sweep(w0, 2L, colSums(w0), "/")
  neff0 <- 1 / colSums(w0 * w0)

  # Base smooth-width rule used before matching to the MICE candidate width.
  h <- h_rule * pmax(donors / pmax(neff0, 1e-8), 1e-4)^0.75
  lambda <- 0.5 / pmax(h, h_floor)^2

  list(delta = delta, lambda = lambda)
}

bin_yhat_donor_pmmrw <- function(yhat, y, row_id, k, edges) {
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

draw_neighbors_pmmrw <- function(yhat, edges, lookup, donor_lookup) {
  n <- length(yhat)
  nbins <- length(edges) - 1L

  bin <- findInterval(yhat, edges, rightmost.closed = TRUE, all.inside = TRUE)
  t0 <- edges[pmax(bin, 1L)]
  t1 <- edges[pmin(bin + 1L, nbins)]
  p_left <- ifelse(t1 > t0, (t1 - yhat) / (t1 - t0), 0.5)
  p_left <- pmin(pmax(p_left, 0), 1)

  selected_bin <- ifelse(stats::runif(n) < p_left, bin, pmin(bin + 1L, nbins))
  indices <- sample(seq_len(ncol(lookup)), n, replace = TRUE)

  list(
    value = lookup[cbind(selected_bin, indices)],
    donor_id = donor_lookup[cbind(selected_bin, indices)],
    bin = bin,
    p_left = p_left,
    selected_bin = selected_bin
  )
}

#' Predictive mean matching for RW-PMM
#'
#' Uses the same PMM imputations as `mice::mice.impute.pmm()` while storing the
#' donor and working-score information needed for RW-PMM variance estimation
#' with `with_rw()` and `pool_rw()`.
#'
#' @inheritParams mice::mice.impute.pmm
#' @param donors The PMM donor-pool size.
#' @param matchtype,nbins,use.matcher PMM matching options.
#' @details The completed values follow standard MICE PMM. The additional
#'   working score is used only by `pool_rw()` for RW-PMM variance estimation.
#'   This method supports numeric variables and requires `tasks = "train"`.
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
    stop("`model` must be an environment; use `tasks = 'train'.", call. = FALSE)
  }

  n <- length(y)
  row_id <- seq_len(n)
  x <- cbind(`(Intercept)` = 1, as.matrix(x))
  norm_draw <- utils::getFromNamespace(".norm.draw", "mice")
  parm <- norm_draw(y, ry, x, ridge = ridge, ...)
  beta_est <- drop(parm$coef)
  beta_draw <- drop(parm$beta)
  beta_hat <- beta_est
  beta_dot <- beta_draw

  if (matchtype == 0L) {
    beta_dot <- beta_hat
  } else if (matchtype == 2L) {
    beta_hat <- beta_dot
  }

  x_obs <- x[ry, , drop = FALSE]
  x_mis <- x[wy, , drop = FALSE]
  yhat_obs <- as.vector(x_obs %*% beta_hat)
  yhat_mis <- as.vector(x_mis %*% beta_dot)

  if (is.null(nbins)) {
    nbins <- round(4 * log(length(yhat_obs)) + 1.5)
  }
  nbins <- max(2L, min(nbins, length(unique(yhat_obs))))
  if (is.null(donors)) {
    donors <- round(length(yhat_obs) / 600 + 7)
  }
  donors <- max(1L, min(donors, length(yhat_obs)))
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
    donor_lookup = lookup$donor_id
  )
  donor_pos <- match(draws$donor_id, row_id[ry])

  smooth <- smooth_pmm_working_scale(yhat_obs, yhat_mis, donors)

  picked <- cbind(donor_pos, seq_len(sum(wy)))
  # Anchor the smooth distance score to the candidate law actually used by
  # MICE: first choose an adjacent bin, then sample uniformly from its lookup.
  right_bin <- pmin(draws$bin + 1L, nrow(lookup$donor_id))
  left_id <- lookup$donor_id[draws$bin, , drop = FALSE]
  right_id <- lookup$donor_id[right_bin, , drop = FALSE]
  left_pos <- matrix(match(left_id, row_id[ry]), nrow = nrow(left_id))
  right_pos <- matrix(match(right_id, row_id[ry]), nrow = nrow(right_id))

  n_mis <- nrow(x_mis)
  left_delta <- matrix(
    yhat_obs[left_pos] - rep(yhat_mis, ncol(left_pos)), nrow = n_mis
  )
  right_delta <- matrix(
    yhat_obs[right_pos] - rep(yhat_mis, ncol(right_pos)), nrow = n_mis
  )
  mean_derivative <- matrix(0, n_mis, ncol(x_mis))
  for (j in seq_len(ncol(x_mis))) {
    left_x <- matrix(x_obs[, j][left_pos], nrow = n_mis)
    right_x <- matrix(x_obs[, j][right_pos], nrow = n_mis)
    left_mean <- rowMeans(2 * left_delta *
      (left_x - x_mis[, j]))
    right_mean <- rowMeans(2 * right_delta *
      (right_x - x_mis[, j]))
    mean_derivative[, j] <- draws$p_left * left_mean +
      (1 - draws$p_left) * right_mean
  }
  donor_width2 <- draws$p_left * rowMeans(left_delta^2) +
    (1 - draws$p_left) * rowMeans(right_delta^2)
  realized_derivative <- sweep(
    x_obs[donor_pos, , drop = FALSE] - x_mis,
    1L,
    2 * smooth$delta[picked],
    "*"
  )
  # The working width cannot be narrower than the realized MICE candidate-pool
  # width; the existing three-quarter rule remains the lower-width fallback.
  hard_lambda <- 0.5 / pmax(
    donor_width2,
    0.5 / smooth$lambda
  )
  score_mis <- sweep(realized_derivative - mean_derivative, 1L,
                     -hard_lambda, "*")

  # Add the score of MICE's smooth adjacent-bin choice. Quantile-edge
  # derivatives are computed with the same type-7 interpolation used above.
  probs <- seq(0, 1, length.out = nbins + 1L)
  hq <- (nrow(x_obs) - 1) * probs + 1
  lo <- floor(hq)
  hi <- ceiling(hq)
  gamma <- hq - lo
  x_ordered <- x_obs[order(yhat_obs), , drop = FALSE]
  edge_x <- (1 - gamma) * x_ordered[lo, , drop = FALSE] +
    gamma * x_ordered[hi, , drop = FALSE]
  t0_x <- edge_x[draws$bin, , drop = FALSE]
  t1_x <- edge_x[pmin(draws$bin + 1L, nbins), , drop = FALSE]
  t0 <- edges[pmax(draws$bin, 1L)]
  t1 <- edges[pmin(draws$bin + 1L, nbins)]
  interval_width <- t1 - t0
  same_bin <- draws$bin == pmin(draws$bin + 1L, nbins)
  interval_width[same_bin] <- 1
  dp <- sweep(
    sweep(t1_x, 1L, 1 - draws$p_left, "*") +
      sweep(t0_x, 1L, draws$p_left, "*") - x_mis,
    1L, interval_width, "/"
  )
  left_selected <- draws$selected_bin == draws$bin
  bin_score <- dp
  bin_score[left_selected, ] <- sweep(
    dp[left_selected, , drop = FALSE],
    1L, draws$p_left[left_selected], "/"
  )
  bin_score[!left_selected, ] <- sweep(
    -dp[!left_selected, , drop = FALSE],
    1L, 1 - draws$p_left[!left_selected], "/"
  )
  bin_score[same_bin, ] <- 0
  score_mis <- score_mis + bin_score

  sigma_hat <- as.numeric(parm$sigma)
  residual_obs <- y[ry] - as.vector(x_obs %*% beta_est)
  score_orig <- x_obs * (residual_obs / sigma_hat^2)
  information <- -crossprod(x_obs) / (sigma_hat^2 * n)
  d_obs <- t(-solve(information, t(score_orig)))

  pmm_score <- matrix(0, n, ncol(x), dimnames = list(NULL, colnames(x)))
  pmm_d <- matrix(0, n, ncol(x), dimnames = list(NULL, colnames(x)))
  pmm_score[wy, ] <- score_mis
  pmm_d[ry, ] <- d_obs

  donor_id <- rep(NA_integer_, n)
  donor_id[wy] <- draws$donor_id

  model$setup <- list(
    method = "pmmrw", n = sum(ry), task = task, donors = donors,
    matchtype = matchtype, nbins = nbins, ridge = ridge
  )
  model$beta.hat <- beta_hat
  model$beta.dot <- beta_dot
  model$edges <- edges
  model$lookup <- lookup$value
  model$donor_lookup <- lookup$donor_id
  model$sigma.dot <- sigma_hat
  model$xnames <- colnames(x)
  model$donor_id <- donor_id
  model$pmm_score <- pmm_score
  model$pmm_d <- pmm_d
  model$pmm_diagnostics <- c(
    mean_base_width = mean(sqrt(0.5 / smooth$lambda)),
    mean_donor_width = mean(sqrt(donor_width2)),
    mean_working_width = mean(sqrt(0.5 / hard_lambda))
  )

  draws$value
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

  if (length(out) == 1L) out[[1L]] else out
}
