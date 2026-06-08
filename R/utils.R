
validate_mids_input <- function(data) {
  if (!inherits(data, "mids")) {
    cli::cli_abort(
      "{.arg data} must be a {.cls mids} object from {.pkg mice}, not {.cls {class(data)}}."
    )
  }
  
  if (is.null(data$models)) {
    cli::cli_abort(c(
      "{.arg data} must contain imputation models.",
      "i" = "Run {.fn mice} with {.code tasks = 'train'} to save models."
    ))
  }
  
  validate_imputation_methods(data)
}

validate_imputation_methods <- function(data) {
  supported_methods <- c("norm", "logreg", "pmmrw")
  
  imputed_vars <- names(data$models)
  methods_used <- data$method[imputed_vars]
  
  unsupported <- setdiff(methods_used, supported_methods)
  if (length(unsupported) > 0) {
    cli::cli_abort(c(
      "Unsupported imputation method{?s}: {.val {unsupported}}",
      "i" = "Robins-Wang variance currently only supports: {.val {supported_methods}}",
      "x" = "Other methods may not store the required variance information."
    ))
  }
}
