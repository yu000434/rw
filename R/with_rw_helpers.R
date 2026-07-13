extract_model_family <- function(model) {
  if (inherits(model, "lm") && !inherits(model, "glm")) {
    "gaussian"
  } else if (inherits(model, "glm")) {
    stats::family(model)$family
  } else {
    cli::cli_abort("Unsupported model class: {.cls {class(model)}}")
  }
}

extract_single_model <- function(data, var, p) {
  mod <- data$models[[var]][[p]]
  beta <- stats::setNames(as.numeric(mod$beta.dot), mod$xnames)
  imp_method <- data$method[var]
  
  if (imp_method == "logreg") {
    list(method = imp_method, family = "binomial", coefficients = beta)
  } else if (imp_method == "norm") {
    sigma2 <- extract_sigma2(mod, var)
    list(method = imp_method, family = "gaussian", coefficients = beta, sigma2 = sigma2)
  } else if (imp_method == "pmmrw") {
    list(
      method = imp_method,
      family = "gaussian",
      coefficients = beta,
      pmm_score = mod$pmm_score,
      pmm_d = mod$pmm_d
    )
  } else {
    cli::cli_abort(
      "Unsupported imputation method {.val {imp_method}} for variable {.field {var}}."
    )
  }
}

extract_imputation_models <- function(data, imputed_vars, m) {
  model_list <- vector("list", m)
  
  for (p in seq_len(m)) {
    model_list[[p]] <- lapply(imputed_vars, function(var) {
      extract_single_model(data, var, p)
    })
    names(model_list[[p]]) <- imputed_vars
  }
  
  model_list
}

extract_sigma2 <- function(mod, var) {
  if (!is.null(mod$sigma2)) {
    as.numeric(mod$sigma2)
  } else if (!is.null(mod$sigma.dot)) {
    as.numeric(mod$sigma.dot)^2
  } else if (!is.null(mod$sigma)) {
    as.numeric(mod$sigma)^2
  } else {
    cli::cli_abort(c(
      "Cannot extract {.field sigma2} from imputation model for {.field {var}}.",
      "i" = "Model object is missing {.field sigma2}, {.field sigma.dot}, and {.field sigma}."
    ))
  }
}

convert_logreg_factors <- function(data.i, imputed_vars) {
  for (var in imputed_vars) {
    if (is.factor(data.i[[var]])) {
      levels_char <- levels(data.i[[var]])
      if (all(levels_char %in% c("0", "1"))) {
        data.i[[var]] <- as.numeric(as.character(data.i[[var]]))
      } else {
        data.i[[var]] <- as.numeric(data.i[[var]]) - 1
      }
    }
  }
  data.i
}

fit_analysis_model <- function(expr, data.i, envir) {
  mod_analysis <- eval(expr = expr, envir = data.i, enclos = envir)
  if (is.expression(mod_analysis)) {
    mod_analysis <- eval(expr = mod_analysis, envir = data.i, enclos = envir)
  }
  
  if (!inherits(mod_analysis, c("lm", "glm"))) {
    cli::cli_abort(c(
      "Analysis model must be {.cls lm} or {.cls glm}, not {.cls {class(mod_analysis)}}.",
      "i" = "Use {.fn lm} or {.fn glm} for analysis models."
    ))
  }
  
  mod_analysis
}


compute_lm_components <- function(mod_analysis, X_analysis) {
  y_analysis <- stats::model.response(stats::model.frame(mod_analysis))
  pred <- stats::fitted(mod_analysis)
  resid <- y_analysis - pred
  sigma2 <- summary(mod_analysis)$sigma^2
  U_analysis <- X_analysis * resid / sigma2
  tau <- -crossprod(X_analysis) / sigma2
  
  list(U = U_analysis, tau = tau)
}

compute_glm_components <- function(mod_analysis, X_analysis) {
  y_analysis <- stats::model.response(stats::model.frame(mod_analysis))
  pred <- stats::fitted(mod_analysis)
  U_analysis <- X_analysis * (y_analysis - pred)
  
  mu <- stats::predict(mod_analysis, type = "response")
  w <- stats::family(mod_analysis)$mu.eta(
    stats::predict(mod_analysis, type = "link"))^2 /
    stats::family(mod_analysis)$variance(mu)
  tau <- -crossprod(X_analysis * sqrt(w))
  
  list(U = U_analysis, tau = tau)
}

extend_score_to_full_data <- function(U_analysis, mf, data.i, n, mod_analysis) {
  U_imp <- matrix(0, n, ncol(U_analysis))
  colnames(U_imp) <- names(stats::coef(mod_analysis))
  
  analysis_ids <- as.numeric(rownames(mf))
  if (is.null(analysis_ids) || any(is.na(analysis_ids))) {
    analysis_ids <- match(rownames(mf), rownames(data.i))
  }
  
  U_imp[analysis_ids, ] <- U_analysis
  U_imp
}

compute_analysis_components <- function(mod_analysis, data.i, n) {
  X_analysis <- stats::model.matrix(mod_analysis)
  mf <- stats::model.frame(mod_analysis)
  
  if (inherits(mod_analysis, "lm") && !inherits(mod_analysis, "glm")) {
    components <- compute_lm_components(mod_analysis, X_analysis)
  } else {
    components <- compute_glm_components(mod_analysis, X_analysis)
  }
  
  U_imp <- extend_score_to_full_data(components$U, mf, data.i, n, mod_analysis)
  
  list(
    U = U_imp,
    tau = components$tau,
    n_analysis = nrow(mf)
  )
}

prepare_imputed_data <- function(data, p, imputed_vars, n) {
  data.i <- mice::complete(data, action = p)
  data.i <- convert_logreg_factors(data.i, imputed_vars)
  
  # Add variable-specific imputation indicators
  imputed_masks <- compute_imputation_masks(data, imputed_vars, n)
  for (var in imputed_vars) {
    data.i[[paste0(".imputed_", var)]] <- imputed_masks[, var]
  }
  
  data.i
}

process_imputed_datasets <- function(data, expr, model_list, imputed_vars, m, n, envir) {
  results <- vector("list", m)
  
  for (p in seq_len(m)) {
    data.i <- prepare_imputed_data(data, p, imputed_vars, n)
    mod_analysis <- fit_analysis_model(expr, data.i, envir)
    
    imputation_components <- compute_imputation_components(
      data.i, model_list[[p]], imputed_vars
    )
    
    analysis_components <- compute_analysis_components(mod_analysis, data.i, n)
    
    results[[p]] <- c(
      list(model = mod_analysis),
      analysis_components,
      imputation_components
    )
  }
  
  results
}
