softmax_columns_pmmrw <- function(log_weight) {
  z <- sweep(log_weight, 2L, apply(log_weight, 2L, max), "-")
  ez <- exp(z)
  sweep(ez, 2L, colSums(ez), "/")
}

smooth_pmm_weights <- function(eta_obs, eta_mis, donors) {
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
  w0 <- softmax_columns_pmmrw(-sweep(dist2, 2L, lambda0, "*"))
  neff0 <- 1 / colSums(w0 * w0)

  # Current smooth PMM bandwidth rule.
  h <- h_rule * pmax(donors / pmax(neff0, 1e-8), 1e-4)^0.75
  lambda <- 0.5 / pmax(h, h_floor)^2
  weight <- softmax_columns_pmmrw(-sweep(dist2, 2L, lambda, "*"))

  list(delta = delta, weight = weight, lambda = lambda)
}

draw_smooth_pmmrw <- function(weight) {
  cum_weight <- apply(weight, 2L, cumsum)
  donor <- colSums(sweep(cum_weight, 2L, stats::runif(ncol(weight)), "<")) + 1L
  pmin(donor, nrow(weight))
}

#' Smooth predictive mean matching for RW-PMM
#'
#' Selects observed donors using a smooth PMM donor-selection law and uses their
#' observed values as imputations. The selected donor rows and smooth PMM
#' working scores are stored for RW-PMM variance estimation with `with_rw()` and
#' `pool_rw()`.
#'
#' @inheritParams mice::mice.impute.pmm
#' @param donors Target donor width `k` used in the smooth bandwidth rule.
#' @param matchtype,nbins,use.matcher Retained for compatibility with the
#'   `mice` PMM interface; they do not replace the current smooth working law.
#' @details This method supports numeric variables and requires
#'   `tasks = "train"`.
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
  parm <- mice:::.norm.draw(y, ry, x, ridge = ridge, ...)
  beta_hat <- drop(parm$coef)
  beta_draw <- drop(parm$beta)

  x_obs <- x[ry, , drop = FALSE]
  x_mis <- x[wy, , drop = FALSE]
  eta_obs <- as.vector(x_obs %*% beta_draw)
  eta_mis <- as.vector(x_mis %*% beta_draw)
  smooth <- smooth_pmm_weights(eta_obs, eta_mis, donors)
  donor_pos <- draw_smooth_pmmrw(smooth$weight)

  picked <- cbind(donor_pos, seq_len(sum(wy)))
  two_delta <- 2 * smooth$delta
  weighted_derivative <- two_delta * smooth$weight
  mean_derivative <- crossprod(weighted_derivative, x_obs) -
    sweep(x_mis, 1L, colSums(weighted_derivative), "*")
  realized_derivative <- sweep(
    x_obs[donor_pos, , drop = FALSE] - x_mis,
    1L,
    2 * smooth$delta[picked],
    "*"
  )
  score_mis <- sweep(realized_derivative - mean_derivative, 1L,
                     -smooth$lambda, "*")

  sigma_hat <- as.numeric(parm$sigma)
  residual_obs <- y[ry] - as.vector(x_obs %*% beta_hat)
  score_orig <- x_obs * (residual_obs / sigma_hat^2)
  information <- -crossprod(x_obs) / (sigma_hat^2 * n)
  d_obs <- t(-solve(information, t(score_orig)))

  pmm_score <- matrix(0, n, ncol(x), dimnames = list(NULL, colnames(x)))
  pmm_d <- matrix(0, n, ncol(x), dimnames = list(NULL, colnames(x)))
  pmm_score[wy, ] <- score_mis
  pmm_d[ry, ] <- d_obs

  donor_id <- rep(NA_integer_, n)
  donor_id[wy] <- row_id[ry][donor_pos]

  model$setup <- list(
    method = "pmmrw", n = sum(ry), task = task, donors = donors,
    matchtype = matchtype, nbins = nbins, ridge = ridge
  )
  model$beta.hat <- beta_hat
  model$beta.dot <- beta_draw
  model$sigma.dot <- sigma_hat
  model$xnames <- colnames(x)
  model$donor_id <- donor_id
  model$pmm_score <- pmm_score
  model$pmm_d <- pmm_d

  y[ry][donor_pos]
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
