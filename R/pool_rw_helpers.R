compute_omega <- function(u_bar, n) {
  crossprod(u_bar) / n
}

normalize_donor_id <- function(donor_id, m, n) {
  if (is.null(donor_id)) {
    return(NULL)
  }

  if (is.data.frame(donor_id)) {
    donor_id <- as.matrix(donor_id)
  }

  if (!is.matrix(donor_id) || nrow(donor_id) != n || ncol(donor_id) != m) {
    cli::cli_abort("{.arg donor_id} must be an {n} by {m} matrix.")
  }

  x_num <- suppressWarnings(matrix(as.numeric(donor_id), nrow = n, ncol = m))
  bad_na <- is.na(x_num) & !is.na(donor_id)
  bad_int <- !is.na(x_num) & x_num != floor(x_num)
  if (any(bad_na | bad_int)) {
    cli::cli_abort("{.arg donor_id} values must be integer row numbers or {.code NA}.")
  }

  bad_range <- !is.na(x_num) & (x_num < 1L | x_num > n)
  if (any(bad_range)) {
    cli::cli_abort("{.arg donor_id} values must be row numbers between 1 and {n}.")
  }

  out <- matrix(as.integer(x_num), nrow = n, ncol = m)
  dimnames(out) <- dimnames(donor_id)
  out
}

cluster_U_by_donor <- function(U, donor_id) {
  U_cluster <- as.matrix(U)
  recipients <- which(!is.na(donor_id))

  for (i in recipients) {
    j <- donor_id[[i]]
    U_cluster[j, ] <- U_cluster[j, ] + U_cluster[i, ]
    U_cluster[i, ] <- 0
  }

  U_cluster
}

compute_variance_components <- function(results, m, n) {
  kappa_sum <- 0
  alpha_sum <- 0
  d_bar_sum <- 0
  
  for (p in seq_len(m)) {
    U <- results[[p]]$U
    S_mis_imp <- results[[p]]$S_mis_imp
    d <- results[[p]]$d
    
    kappa_sum <- kappa_sum + crossprod(U, S_mis_imp)
    alpha_sum <- alpha_sum + crossprod(d)
    d_bar_sum <- d_bar_sum + d
  }
  
  #d <- results[[1]]$d
  
  list(
    kappa = kappa_sum / (n * m),
    alpha = alpha_sum / (n * m),
  #  alpha = crossprod(d) / n,
    d_bar = d_bar_sum / m
  #  d_bar = d
  )
}

compute_rw_variance <- function(results, m, n, donor_id = NULL) {
  donor_id <- normalize_donor_id(donor_id, m, n)

  U_sum <- Reduce(`+`, lapply(results, function(r) r$U))
  U_omega_sum <- Reduce(`+`, lapply(seq_len(m), function(p) {
    U <- results[[p]]$U
    if (is.null(donor_id)) U else cluster_U_by_donor(U, donor_id[, p])
  }))
  tau_sum <- Reduce(`+`, lapply(results, function(r) r$tau))
  
  u_bar <- U_sum / m
  u_bar_omega <- U_omega_sum / m
  omega <- compute_omega(u_bar_omega, n)
  
  components <- compute_variance_components(results, m, n)
  
  delta <- omega + 
    components$kappa %*% components$alpha %*% t(components$kappa) +
    (1 / n) * (components$kappa %*% t(components$d_bar) %*% u_bar_omega + 
                 t(components$kappa %*% t(components$d_bar) %*% u_bar_omega))
  
  tau <- tau_sum / (m * n)
  tau_inv <- solve(tau)
  (1 / n) * tau_inv %*% delta %*% t(tau_inv)
}

pool_estimates <- function(results) {
  template <- stats::coef(results[[1]]$model)
  coefs <- vapply(results, function(r) stats::coef(r$model),
                  numeric(length(template)))
  if (is.null(dim(coefs))) {
    coefs <- matrix(coefs, nrow = length(template))
    rownames(coefs) <- names(template)
  }
  rowMeans(coefs)
}

construct_pooled_output <- function(est, se, df, family) {
  test_stat <- est / se
  
  if (family == "gaussian") {
    crit_val <- stats::qt(0.975, df = df)
    p_value <- 2 * stats::pt(-abs(test_stat), df = df)
  } else {
    crit_val <- stats::qnorm(0.975)  
    p_value <- 2 * stats::pnorm(-abs(test_stat))
  }
  out <- data.frame(
    term = names(est),
    estimate = est,
    std.error = se,
    statistic = test_stat,
    p.value = p_value,
    conf.low = est - crit_val * se,
    conf.high = est + crit_val * se
  )
  out
}
