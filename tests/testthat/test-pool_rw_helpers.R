test_that("pool_estimates averages coefficients across imputations", {
  mod1 <- lm(mpg ~ wt, data = mtcars[1:10, ])
  mod2 <- lm(mpg ~ wt, data = mtcars[11:20, ])
  mod3 <- lm(mpg ~ wt, data = mtcars[21:30, ])
  
  results <- list(
    list(model = mod1),
    list(model = mod2),
    list(model = mod3)
  )
  
  result <- pool_estimates(results)
  
  expect_length(result, 2)
  expect_named(result, c("(Intercept)", "wt"))
  expect_equal(result[["(Intercept)"]], mean(c(coef(mod1)[1], coef(mod2)[1], coef(mod3)[1])))
  expect_equal(result[["wt"]], mean(c(coef(mod1)[2], coef(mod2)[2], coef(mod3)[2])))
})

test_that("pool_estimates handles single imputation", {
  mod <- lm(mpg ~ wt + hp, data = mtcars)
  results <- list(list(model = mod))
  
  result <- pool_estimates(results)
  
  expect_equal(result, coef(mod))
})

test_that("pool_estimates preserves coefficient names", {
  mod1 <- lm(mpg ~ wt + hp + cyl, data = mtcars[1:16, ])
  mod2 <- lm(mpg ~ wt + hp + cyl, data = mtcars[17:32, ])
  
  results <- list(
    list(model = mod1),
    list(model = mod2)
  )
  
  result <- pool_estimates(results)
  
  expect_named(result, c("(Intercept)", "wt", "hp", "cyl"))
})

test_that("compute_omega calculates crossprod correctly", {
  u_bar <- matrix(c(1, 2, 3, 4), nrow = 10, ncol = 2, byrow = TRUE)
  n <- 10
  
  result <- compute_omega(u_bar, n)
  
  expected <- crossprod(u_bar) / n
  expect_equal(result, expected)
})

test_that("compute_omega returns symmetric matrix", {
  u_bar <- matrix(rnorm(30), nrow = 10, ncol = 3)
  n <- 10
  
  result <- compute_omega(u_bar, n)
  
  expect_equal(result, t(result))
})

test_that("compute_omega handles single column", {
  u_bar <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
  n <- 5
  
  result <- compute_omega(u_bar, n)
  
  expect_equal(dim(result), c(1, 1))
  expect_equal(result[1, 1], sum(u_bar^2) / n)
})

test_that("cluster_U_by_donor aggregates recipient scores to donor rows", {
  U <- matrix(
    c(1, 0,
      2, 0,
      3, 0,
      4, 0),
    nrow = 4,
    byrow = TRUE
  )
  donor_id <- c(NA, 1, NA, 1)

  result <- cluster_U_by_donor(U, donor_id)

  expected <- matrix(
    c(7, 0,
      0, 0,
      3, 0,
      0, 0),
    nrow = 4,
    byrow = TRUE
  )
  expect_equal(result, expected)
})

test_that("normalize_donor_id accepts an n by m matrix", {
  donor_id <- matrix(c(NA, 1, NA, NA, 2, NA), nrow = 3)
  expect_equal(normalize_donor_id(donor_id, m = 2, n = 3), donor_id)
  expect_error(normalize_donor_id(c(NA, 1, NA), m = 1, n = 3),
               "must be an 3 by 1 matrix")
  expect_error(normalize_donor_id(matrix(c(NA, 1.5, NA), nrow = 3), m = 1, n = 3),
               "integer row numbers")
  expect_error(normalize_donor_id(matrix(c(NA, 4, NA), nrow = 3), m = 1, n = 3),
               "between 1 and 3")
})

test_that("compute_variance_components calculates all components", {
  results <- list(
    list(
      U = matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2),
      S_mis_imp = matrix(c(0.5, 0.6, 0.7, 0.8), nrow = 2, ncol = 2),
      d = matrix(c(0.1, 0.2), nrow = 1, ncol = 2)
    ),
    list(
      U = matrix(c(1.5, 2.5, 3.5, 4.5), nrow = 2, ncol = 2),
      S_mis_imp = matrix(c(0.6, 0.7, 0.8, 0.9), nrow = 2, ncol = 2),
      d = matrix(c(0.15, 0.25), nrow = 1, ncol = 2)
    )
  )
  m <- 2
  n <- 2
  
  result <- compute_variance_components(results, m, n)
  
  expect_named(result, c("kappa", "alpha", "d_bar"))
  expect_equal(dim(result$kappa), c(2, 2))
  expect_equal(dim(result$alpha), c(2, 2))
  expect_equal(dim(result$d_bar), c(1, 2))
})

test_that("compute_variance_components kappa is average of crossproducts", {
  U1 <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  S1 <- matrix(c(2, 0, 0, 2), nrow = 2, ncol = 2)
  U2 <- matrix(c(3, 0, 0, 3), nrow = 2, ncol = 2)
  S2 <- matrix(c(4, 0, 0, 4), nrow = 2, ncol = 2)
  
  results <- list(
    list(U = U1, S_mis_imp = S1, d = matrix(0, 1, 2)),
    list(U = U2, S_mis_imp = S2, d = matrix(0, 1, 2))
  )
  
  result <- compute_variance_components(results, 2, 2)
  
  expected_kappa <- (crossprod(U1, S1) + crossprod(U2, S2)) / (2 * 2)
  expect_equal(result$kappa, expected_kappa)
})

test_that("compute_variance_components alpha is average of d crossproducts", {
  d1 <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
  d2 <- matrix(c(2, 3, 4), nrow = 1, ncol = 3)
  
  results <- list(
    list(U = matrix(0, 1, 1), S_mis_imp = matrix(0, 1, 3), d = d1),
    list(U = matrix(0, 1, 1), S_mis_imp = matrix(0, 1, 3), d = d2)
  )
  
  result <- compute_variance_components(results, 2, 1)
  
  expected_alpha <- (crossprod(d1) + crossprod(d2)) / (1 * 2)
  expect_equal(result$alpha, expected_alpha)
})

test_that("compute_variance_components d_bar is average of d", {
  d1 <- matrix(c(1, 2), nrow = 1, ncol = 2)
  d2 <- matrix(c(3, 4), nrow = 1, ncol = 2)
  d3 <- matrix(c(5, 6), nrow = 1, ncol = 2)
  
  results <- list(
    list(U = matrix(0, 1, 1), S_mis_imp = matrix(0, 1, 2), d = d1),
    list(U = matrix(0, 1, 1), S_mis_imp = matrix(0, 1, 2), d = d2),
    list(U = matrix(0, 1, 1), S_mis_imp = matrix(0, 1, 2), d = d3)
  )
  
  result <- compute_variance_components(results, 3, 1)
  
  expect_equal(result$d_bar, matrix(c(3, 4), nrow = 1, ncol = 2))
})

test_that("construct_pooled_output creates correct data frame structure", {
  est <- c(intercept = 1.5, x = 2.0)
  se <- c(intercept = 0.5, x = 0.3)
  
  result <- construct_pooled_output(est, se, df = 100, family = "gaussian")
  
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"))
  expect_equal(nrow(result), 2)
})

test_that("construct_pooled_output calculates statistics correctly", {
  est <- c(x = 2.0)
  se <- c(x = 0.5)
  
  result <- construct_pooled_output(est, se, df = 100, family = "gaussian")
  
  expect_equal(result$estimate, 2.0)
  expect_equal(result$std.error, 0.5)
  expect_equal(result$statistic, 4.0)
})

test_that("construct_pooled_output calculates p-values correctly", {
  est <- c(x = 0.0)
  se <- c(x = 1.0)
  
  result <- construct_pooled_output(est, se, df = 100, family = "gaussian")
  
  expect_equal(result$p.value, 1.0)
})

test_that("construct_pooled_output calculates confidence intervals correctly", {
  est <- c(x = 10.0)
  se <- c(x = 2.0)
  
  result <- construct_pooled_output(est, se, df = NULL, family = "binomial")
  
  expect_equal(result$conf.low, 10.0 - stats::qnorm(0.975) * 2.0)
  expect_equal(result$conf.high, 10.0 + stats::qnorm(0.975) * 2.0)
})


test_that("construct_pooled_output preserves term names", {
  est <- c("(Intercept)" = 5.0, age = 0.5, bmi = 1.2)
  se <- c("(Intercept)" = 1.0, age = 0.1, bmi = 0.3)
  
  result <- construct_pooled_output(est, se, df = 100, family = "gaussian")
  
  expect_equal(result$term, c("(Intercept)", "age", "bmi"))
})

test_that("compute_rw_variance returns symmetric matrix", {
  U <- matrix(c(1, 2, 3, 4), nrow = 10, ncol = 2, byrow = TRUE)
  tau <- matrix(c(10, 2, 2, 8), nrow = 2, ncol = 2)
  S_mis_imp <- matrix(c(0.5, 0.6), nrow = 10, ncol = 2, byrow = TRUE)
  d <- matrix(c(0.1, 0.2), nrow = 10, ncol = 2, byrow = TRUE)
  
  results <- list(
    list(U = U, tau = tau, S_mis_imp = S_mis_imp, d = d),
    list(U = U * 1.1, tau = tau * 1.05, S_mis_imp = S_mis_imp * 0.9, d = d * 1.2)
  )
  
  result <- compute_rw_variance(results, 2, 10)
  
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("compute_rw_variance has correct dimensions", {
  p <- 2
  n <- 20
  m <- 3
  
  U <- matrix(rnorm(n * p), nrow = n, ncol = p)
  tau <- crossprod(matrix(rnorm(n * p), nrow = n, ncol = p))
  S_mis_imp <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  d <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  
  results <- lapply(1:m, function(i) {
    list(U = U, tau = tau, S_mis_imp = S_mis_imp, d = d)
  })
  
  result <- compute_rw_variance(results, m, n)
  
  expect_equal(dim(result), c(p, p))
})

test_that("compute_rw_variance handles single imputation", {
  U <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  tau <- matrix(c(10, 1, 1, 10), nrow = 2, ncol = 2)
  S_mis_imp <- matrix(c(0.5, 0.6, 0.7, 0.8), nrow = 2, ncol = 2)
  d <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, ncol = 2)
  
  results <- list(
    list(U = U, tau = tau, S_mis_imp = S_mis_imp, d = d)
  )
  
  result <- compute_rw_variance(results, 1, 2)
  
  expect_equal(dim(result), c(2, 2))
  expect_equal(result, t(result))
})

test_that("compute_rw_variance uses donor-clustered omega when donor_id is supplied", {
  n <- 3
  U <- matrix(
    c(1, 0,
      2, 0,
      3, 0),
    nrow = n,
    byrow = TRUE
  )
  U_cluster <- matrix(
    c(3, 0,
      0, 0,
      3, 0),
    nrow = n,
    byrow = TRUE
  )
  tau <- diag(2) * n
  S_mis_imp <- matrix(0, nrow = n, ncol = 1)
  d <- matrix(0, nrow = n, ncol = 1)
  results <- list(list(U = U, tau = tau, S_mis_imp = S_mis_imp, d = d))

  donor_id <- matrix(c(NA, 1, NA), ncol = 1)
  result <- compute_rw_variance(results, m = 1, n = n, donor_id = donor_id)

  expected_omega <- crossprod(U_cluster) / n
  expected <- expected_omega / n
  expect_equal(result, expected)
})

test_that("donor IDs adjust omega and pair the smooth cross term with donor-source U", {
  n <- 3
  U <- matrix(
    c(1, 2,
      3, 4,
      5, 6),
    nrow = n,
    byrow = TRUE
  )
  U_cluster <- matrix(
    c(4, 6,
      0, 0,
      5, 6),
    nrow = n,
    byrow = TRUE
  )
  tau <- diag(2) * n
  S_mis_imp <- matrix(c(0.2, 0.4, 0.6), ncol = 1)
  d <- matrix(c(0.1, 0.3, 0.5), ncol = 1)
  results <- list(list(U = U, tau = tau, S_mis_imp = S_mis_imp, d = d))
  donor_id <- matrix(c(NA, 1, NA), ncol = 1)

  result <- compute_rw_variance(results, m = 1, n = n, donor_id = donor_id)

  omega <- crossprod(U_cluster) / n
  kappa <- crossprod(U, S_mis_imp) / n
  alpha <- crossprod(d) / n
  cross <- (1 / n) * (
    kappa %*% t(d) %*% U_cluster +
      t(kappa %*% t(d) %*% U_cluster)
  )
  expected_delta <- omega + kappa %*% alpha %*% t(kappa) + cross
  expected <- expected_delta / n

  expect_equal(result, expected)
})

test_that("construct_pooled_output handles negative estimates", {
  est <- c(x = -2.5)
  se <- c(x = 0.5)
  
  result <- construct_pooled_output(est, se, df = NULL, family = "binomial")
  
  expect_equal(result$estimate, -2.5)
  expect_equal(result$statistic, -5.0)
  expect_true(result$p.value < 0.05)
  expect_true(result$conf.low < result$conf.high)
})

test_that("construct_pooled_output handles very small standard errors", {
  est <- c(x = 5.0)
  se <- c(x = 0.001)
  
  result <- construct_pooled_output(est, se, df = NULL, family = "binomial")
  
  expect_true(result$p.value < 0.001)
  expect_true(abs(result$conf.high - result$conf.low) < 0.01)
})

test_that("construct_pooled_output handles very large standard errors", {
  est <- c(x = 1.0)
  se <- c(x = 100.0)
  
  result <- construct_pooled_output(est, se, df = NULL, family = "binomial")
  
  expect_true(result$p.value > 0.9)
  expect_true(abs(result$conf.high - result$conf.low) > 300)
})

test_that("compute_variance_components handles single observation", {
  results <- list(
    list(
      U = matrix(c(1, 2), nrow = 1, ncol = 2),
      S_mis_imp = matrix(c(0.5, 0.6), nrow = 1, ncol = 2),
      d = matrix(c(0.1, 0.2), nrow = 1, ncol = 2)
    )
  )
  
  result <- compute_variance_components(results, 1, 1)
  
  expect_equal(dim(result$kappa), c(2, 2))
  expect_equal(dim(result$alpha), c(2, 2))
  expect_equal(dim(result$d_bar), c(1, 2))
})
