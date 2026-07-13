#' Pool results with Robins-Wang variance
#' 
#' @param object A with_rw object from with_rw()
#' @param donor_id Optional `n` by `m` matrix of PMM donor row IDs for
#'   donor-source analysis-score adjustment. If the imputation used
#'   `method = "pmmrw"` for one variable, these IDs are extracted automatically.
#'   Entries should be `NA` for rows that were not PMM recipients and the
#'   observed donor row number for rows imputed from a PMM donor.
#' @param ... Additional arguments (currently unused)
#' @return A pool_rw object with pooled estimates and Robins-Wang SEs
#' 
#' @examples
#' \dontrun{
#' library(mice)
#' imp <- mice(nhanes, method = "norm", m = 5, tasks = "train", print = FALSE)
#' fit <- with_rw(imp, lm(bmi ~ age + hyp))
#' pooled <- pool_rw(fit)
#' summary(pooled)
#' }
#' @export
pool_rw <- function(object, donor_id = NULL, ...) {
  if (!inherits(object, "with_rw")) {
    cli::cli_abort(
      "{.arg object} must be a {.cls with_rw} object, not {.cls {class(object)}}."
    )
  }
  
  m <- object$m
  n <- object$n
  results <- object$results
  methods <- character()
  if (!is.null(object$mids) && !is.null(object$mids$models)) {
    methods <- object$mids$method[names(object$mids$models)]
  }
  
  pmmrw_vars <- names(methods)[methods == "pmmrw"]
  if (is.null(donor_id) && length(pmmrw_vars) == 1L) {
    donor_id <- extract_donor_id(object$mids, pmmrw_vars)
  } else if (is.null(donor_id) && length(pmmrw_vars) > 1L) {
    cli::cli_abort(c(
      "Multiple {.val pmmrw} variables were found.",
      "i" = "Please supply the donor ID matrix to use with {.code pool_rw(object, donor_id = donor_id)}."
    ))
  }
  
  est <- pool_estimates(results)
  GAMMA <- compute_rw_variance(results, m, n, donor_id = donor_id)
  se <- sqrt(diag(GAMMA))
  family <- extract_model_family(results[[1]]$model)
  
  if (family == "gaussian") {
    df <- stats::df.residual(results[[1]]$model)
  } else {
    df <- NULL
  }
  
  out <- construct_pooled_output(est, se, df, family)
  
  structure(
    list(
      pooled = out,
      variance = GAMMA,
      m = m,
      n = n,
      donor_correlation = !is.null(donor_id),
      call = match.call()
    ),
    class = "pool_rw"
  )
}

#' @export
print.pool_rw <- function(x, ...) {
  cli::cli_h1("Robins-Wang Pooled Results")
  cli::cli_text("Number of imputations: {x$m}")
  cli::cli_text("Sample size: {x$n}")
  cli::cli_text("")
  print(x$pooled, digits = 4)
  invisible(x)
}

#' @export
summary.pool_rw <- function(object, ...) {
  print(object, ...)
}

#' @export
coef.pool_rw <- function(object, ...) {
  stats::setNames(object$pooled$estimate, object$pooled$term)
}

#' @export
vcov.pool_rw <- function(object, ...) {
  object$variance
}
