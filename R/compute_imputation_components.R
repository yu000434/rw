compute_imputation_masks <- function(data, imputed_vars, n) {
  # Return a matrix: n x num_vars, where each column indicates missingness for that variable
  masks <- sapply(imputed_vars, function(var) {
    as.integer(is.na(data$data[[var]]))
  })
  colnames(masks) <- imputed_vars
  masks
}

compute_score <- function(data, model_info, var_name) {
  # Extract and clean coefficient names
  beta <- model_info$coefficients
  nm <- names(beta)
  nm[nm == ""] <- "(Intercept)"
  names(beta) <- nm
  coef_names <- names(beta)
  
  # Build design matrix
  if ("(Intercept)" %in% coef_names) {
    nm_noint <- setdiff(coef_names, "(Intercept)")
    X <- cbind(
      `(Intercept)` = 1, 
      as.matrix(data[, nm_noint, drop = FALSE])
    )
    X <- X[, coef_names, drop = FALSE]
  } else {
    X <- as.matrix(data[, coef_names, drop = FALSE])
  }
  
  y <- data[[var_name]]
  family <- model_info$family[1]
  
  # Family-specific score calculation
  if (family == "gaussian") {
    sigma_est <- as.numeric(model_info$sigma2)
    pred <- as.vector(X %*% beta)
    resid <- y - pred
    S_beta <- sweep(X, 1, resid / sigma_est, `*`)
    S_sigma <- 0.5 * (-1 / sigma_est + resid^2 / sigma_est^2)
    S_u <- cbind(S_beta, S_sigma)
    return(S_u)
    
  } else if (family == "binomial") {
    lp <- as.vector(X %*% beta)
    p <- 1 / (1 + exp(-lp))
    S_u <- sweep(X, 1, y - p, `*`)
    return(S_u)
    
  } else {
    cli::cli_abort(
      "Unsupported family: {.val {family}}. Supported: {.val gaussian}, {.val binomial}."
    )
  }
}

compute_information <- function(data, model_info, var_name, imputed_flag) {
  # Extract and clean coefficient names
  beta <- model_info$coefficients
  nm <- names(beta)
  nm[nm == ""] <- "(Intercept)"
  names(beta) <- nm
  coef_names <- names(beta)
  
  # Build design matrix
  if ("(Intercept)" %in% coef_names) {
    nm_noint <- setdiff(coef_names, "(Intercept)")
    X <- cbind(
      `(Intercept)` = 1, 
      as.matrix(data[, nm_noint, drop = FALSE])
    )
    X <- X[, coef_names, drop = FALSE]
  } else {
    X <- as.matrix(data[, coef_names, drop = FALSE])
  }
  
  y <- data[[var_name]]
  obs_idx <- which(imputed_flag == 0)
  n <- nrow(data)
  family <- model_info$family[1]
  
  # Family-specific information calculation
  if (family == "gaussian") {
    sigma_est <- as.numeric(model_info$sigma2)
    pred <- as.vector(X %*% beta)
    X_obs <- X[obs_idx, , drop = FALSE]
    r_obs <- y[obs_idx] - pred[obs_idx]
    
    # Compute information matrix components
    S2_uu <- (-1 / sigma_est) * crossprod(X_obs)
    S2_u_sigma <- (-1 / sigma_est^2) * crossprod(X_obs, r_obs)
    S2_sigma_sigma <- sum(1 / (2 * sigma_est^2) - r_obs^2 / sigma_est^3)
    
    # Assemble full information matrix
    n_params <- ncol(X)
    S2_mod <- matrix(0, n_params + 1, n_params + 1)
    S2_mod[1:n_params, 1:n_params] <- S2_uu / n
    S2_mod[1:n_params, n_params + 1] <- as.vector(S2_u_sigma) / n
    S2_mod[n_params + 1, 1:n_params] <- as.vector(S2_u_sigma) / n
    S2_mod[n_params + 1, n_params + 1] <- S2_sigma_sigma / n
    
    return(S2_mod)
    
  } else if (family == "binomial") {
    lp <- as.vector(X %*% beta)
    p <- 1 / (1 + exp(-lp))
    w <- p * (1 - p)
    
    # Zero out weights for imputed observations
    w_obs <- w
    w_obs[imputed_flag == 1] <- 0
    
    # Compute information matrix 
    X_weighted <- sweep(X, 1, sqrt(w_obs), `*`)
    S2_mod <- -crossprod(X_weighted) / n
    
    return(S2_mod)
    
  } else {
    cli::cli_abort(
      "Unsupported family: {.val {family}}. Supported: {.val gaussian}, {.val binomial}."
    )
  }
}


compute_score_matrix <- function(data.i, model_list, imputed_vars) {
  S_u_list <- lapply(seq_along(imputed_vars), function(i) {
    compute_score(data.i, model_list[[i]], imputed_vars[i])
  })
  do.call(cbind, S_u_list)
}

compute_information_matrix <- function(data.i, model_list, imputed_vars) {
  S2_list <- lapply(seq_along(imputed_vars), function(i) {
    # Use variable-specific imputation flag
    var_imputed_flag <- data.i[[paste0(".imputed_", imputed_vars[i])]]
    compute_information(data.i, model_list[[i]], imputed_vars[i], var_imputed_flag)
  })
  
  build_block_diagonal(S2_list)
}

build_block_diagonal <- function(S2_list) {
  dims <- vapply(S2_list, ncol, integer(1))
  n_tot <- sum(dims)
  S2 <- matrix(0, n_tot, n_tot)
  idx <- 1
  
  for (i in seq_along(S2_list)) {
    rng <- idx:(idx + dims[i] - 1)
    S2[rng, rng] <- S2_list[[i]]
    idx <- idx + dims[i]
  }
  
  S2
}


compute_imputation_components <- function(data.i, model_list, imputed_vars) {
  # Compute score matrix for all variables
  S_u <- compute_score_matrix(data.i, model_list, imputed_vars)
  
  # Build S_mis_imp and S_orig using variable-specific masks
  n <- nrow(data.i)
  col_idx <- 1
  S_mis_imp_list <- list()
  S_orig_list <- list()
  
  for (i in seq_along(imputed_vars)) {
    var <- imputed_vars[i]
    model_info <- model_list[[i]]
    
    # Determine number of parameters for this variable
    n_params <- length(model_info$coefficients)
    if (model_info$family[1] == "gaussian") {
      n_params <- n_params + 1  # Add 1 for sigma parameter
    }
    
    # Get variable-specific imputation indicator
    var_imputed <- data.i[[paste0(".imputed_", var)]]
    
    # Extract columns for this variable
    var_cols <- col_idx:(col_idx + n_params - 1)
    
    # Create masks: 1 if imputed, 0 if observed
    ImputedMat_var <- matrix(var_imputed == 1, n, n_params, byrow = FALSE)

    if (identical(model_info$method, "pmmrw")) {
      # PMM copies observed donor values. The Gaussian model supplies the
      # matching index, but its normal residual score is not a missing-value
      # density score for the copied value.
      S_mis_imp_list[[i]] <- S_u[, var_cols, drop = FALSE] * 0
      S_orig_list[[i]] <- S_u[, var_cols, drop = FALSE] * 0
    } else {
      # S_mis_imp: scores for imputed observations only
      S_mis_imp_list[[i]] <- S_u[, var_cols, drop = FALSE] * ImputedMat_var

      # S_orig: scores for observed observations only
      S_orig_list[[i]] <- S_u[, var_cols, drop = FALSE] * (1 - ImputedMat_var)
    }
    
    col_idx <- col_idx + n_params
  }
  
  # Combine all variables
  S_mis_imp <- do.call(cbind, S_mis_imp_list)
  S_orig <- do.call(cbind, S_orig_list)
  
  # Compute information matrix and its inverse
  S2 <- compute_information_matrix(data.i, model_list, imputed_vars)
  Dmat <- solve(S2)
  
  d <- t(- Dmat %*% t(S_orig))
  
  list(S_mis_imp = S_mis_imp, d = d)
}
