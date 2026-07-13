
test_that("compute_imputation_masks identifies imputed values", {
  data <- list(
    data = data.frame(
      x = c(1, NA, 3, NA),
      y = c(NA, 2, 3, 4)
    )
  )
  
  result <- compute_imputation_masks(data, c("x", "y"), 4)
  
  expect_equal(result[,1], c(0, 1, 0, 1))
  expect_equal(result[,2], c(1, 0, 0, 0))
})

test_that("compute_imputation_masks handles no missing values", {
  data <- list(
    data = data.frame(
      x = c(1, 2, 3),
      y = c(4, 5, 6)
    )
  )
  
  result <- compute_imputation_masks(data, c("x", "y"), 3)
  
  expect_equal(unname(result), matrix(rep(0, 6), nrow = 3))
})

test_that("compute_imputation_masks handles all missing values", {
  data <- list(
    data = data.frame(
      x = c(NA, NA, NA),
      y = c(NA, NA, NA)
    )
  )
  
  result <- compute_imputation_masks(data, c("x", "y"), 3)
  
  expect_equal(unname(result), matrix(rep(1, 6), nrow = 3))
})

test_that("compute_imputation_mask handles single variable", {
  data <- list(
    data = data.frame(
      x = c(1, NA, 3)
    )
  )
  
  result <- compute_imputation_masks(data, "x", 3)
  
  expect_equal(unname(result), matrix(c(0, 1, 0)))
})

test_that("compute_score handles gaussian family correctly", {
  data <- data.frame(
    x = c(1, 2, 3, 4),
    y = c(2.1, 4.2, 5.9, 8.1)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c("(Intercept)" = 0.1, x = 2.0),
    sigma2 = 1.0
  )
  
  result <- compute_score(data, model_info, "y")
  
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 3)
  expect_true(all(is.finite(result)))
})

test_that("compute_score handles binomial family correctly", {
  data <- data.frame(
    x = c(1, 2, 3, 4),
    y = c(0, 1, 1, 0)
  )
  
  model_info <- list(
    family = "binomial",
    coefficients = c("(Intercept)" = 0.5, x = -0.3)
  )
  
  result <- compute_score(data, model_info, "y")
  
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 2)
  expect_true(all(is.finite(result)))
})

test_that("compute_score handles intercept correctly", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(2, 4, 6)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c("(Intercept)" = 1.0, x = 1.5),
    sigma2 = 0.5
  )
  
  result <- compute_score(data, model_info, "y")
  
  expect_equal(ncol(result), 3)
})


test_that("compute_score rejects unsupported family", {
  data <- data.frame(x = 1:3, y = 1:3)
  
  model_info <- list(
    family = "poisson",
    coefficients = c("(Intercept)" = 1.0, x = 0.5)
  )
  
  expect_error(
    compute_score(data, model_info, "y"),
    "Unsupported family.*poisson"
  )
})

test_that("compute_score gaussian calculates residuals correctly", {
  data <- data.frame(
    x = c(0, 1, 2),
    y = c(1, 3, 5)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c("(Intercept)" = 1.0, x = 2.0),
    sigma2 = 1.0
  )
  
  result <- compute_score(data, model_info, "y")
  
  pred <- 1.0 + 2.0 * data$x
  resid <- data$y - pred
  
  expect_equal(result[, 1], resid / 1.0)
  expect_equal(result[, 2], data$x * resid / 1.0)
})

test_that("compute_score binomial calculates probabilities correctly", {
  data <- data.frame(
    x = c(0, 1),
    y = c(0, 1)
  )
  
  model_info <- list(
    family = "binomial",
    coefficients = c("(Intercept)" = 0.0, x = 0.0)
  )
  
  result <- compute_score(data, model_info, "y")
  
  expect_equal(unname(result[1, 1]), 0 - 0.5)
  expect_equal(unname(result[2, 1]), 1 - 0.5)
})

test_that("compute_score handles multiple predictors", {
  data <- data.frame(
    x1 = c(1, 2, 3),
    x2 = c(4, 5, 6),
    y = c(10, 15, 20)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c("(Intercept)" = 1.0, x1 = 2.0, x2 = 0.5),
    sigma2 = 1.0
  )
  
  result <- compute_score(data, model_info, "y")
  
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
})

test_that("compute_information handles gaussian family correctly", {
  data <- data.frame(
    x = c(1, 2, 3, 4),
    y = c(2, 4, 6, 8)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c("(Intercept)" = 0.0, x = 2.0),
    sigma2 = 1.0
  )
  
  imputed_flag <- c(0, 0, 1, 1)
  
  result <- compute_information(data, model_info, "y", imputed_flag)
  
  expect_equal(dim(result), c(3, 3))
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("compute_information handles binomial family correctly", {
  data <- data.frame(
    x = c(1, 2, 3, 4),
    y = c(0, 1, 1, 0)
  )
  
  model_info <- list(
    family = "binomial",
    coefficients = c("(Intercept)" = 0.0, x = 0.0)
  )
  
  imputed_flag <- c(0, 0, 1, 1)
  
  result <- compute_information(data, model_info, "y", imputed_flag)
  
  expect_equal(dim(result), c(2, 2))
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("compute_information returns symmetric matrix", {
  data <- data.frame(
    x = rnorm(20),
    y = rnorm(20)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c("(Intercept)" = 0.0, x = 1.0),
    sigma2 = 1.0
  )
  
  imputed_flag <- rep(0, 20)
  
  result <- compute_information(data, model_info, "y", imputed_flag)
  
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("compute_information zeros out imputed observations", {
  data <- data.frame(
    x = c(1, 2, 3, 4),
    y = c(0, 1, 1, 0)
  )
  
  model_info <- list(
    family = "binomial",
    coefficients = c("(Intercept)" = 0.0, x = 0.0)
  )
  
  imputed_all <- c(0, 0, 0, 0)
  imputed_half <- c(0, 0, 1, 1)
  
  result_all <- compute_information(data, model_info, "y", imputed_all)
  result_half <- compute_information(data, model_info, "y", imputed_half)
  
  expect_true(all(abs(result_all) >= abs(result_half)))
})

test_that("compute_information rejects unsupported family", {
  data <- data.frame(x = 1:3, y = 1:3)
  
  model_info <- list(
    family = "poisson",
    coefficients = c("(Intercept)" = 1.0, x = 0.5)
  )
  
  expect_error(
    compute_information(data, model_info, "y", rep(0, 3)),
    "Unsupported family.*poisson"
  )
})

test_that("compute_information gaussian includes sigma parameter", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(2, 4, 6)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c("(Intercept)" = 0.0, x = 2.0),
    sigma2 = 1.0
  )
  
  imputed_flag <- rep(0, 3)
  
  result <- compute_information(data, model_info, "y", imputed_flag)
  
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
})

test_that("build_block_diagonal creates correct structure", {
  mat1 <- matrix(c(1, 0, 0, 1), 2, 2)
  mat2 <- matrix(c(2, 0, 0, 2), 2, 2)
  
  result <- build_block_diagonal(list(mat1, mat2))
  
  expect_equal(dim(result), c(4, 4))
  expect_equal(result[1:2, 1:2], mat1)
  expect_equal(result[3:4, 3:4], mat2)
  expect_equal(result[1:2, 3:4], matrix(0, 2, 2))
  expect_equal(result[3:4, 1:2], matrix(0, 2, 2))
})

test_that("build_block_diagonal handles single matrix", {
  mat <- matrix(1:9, 3, 3)
  
  result <- build_block_diagonal(list(mat))
  
  expect_equal(result, mat)
})

test_that("build_block_diagonal handles different sized matrices", {
  mat1 <- matrix(1, 2, 2)
  mat2 <- matrix(2, 3, 3)
  mat3 <- matrix(3, 1, 1)
  
  result <- build_block_diagonal(list(mat1, mat2, mat3))
  
  expect_equal(dim(result), c(6, 6))
  expect_equal(result[1:2, 1:2], mat1)
  expect_equal(result[3:5, 3:5], mat2)
  expect_equal(result[6, 6], 3)
})

test_that("compute_score_matrix combines scores correctly", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(2, 4, 6),
    z = c(0, 1, 1)
  )
  
  model_list <- list(
    list(
      family = "gaussian",
      coefficients = c("(Intercept)" = 0.0, x = 2.0),
      sigma2 = 1.0
    ),
    list(
      family = "binomial",
      coefficients = c("(Intercept)" = 0.0, y = 0.5)
    )
  )
  
  result <- compute_score_matrix(data, model_list, c("y", "z"))
  
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 5)
})

test_that("compute_score_matrix handles single variable", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(2, 4, 6)
  )
  
  model_list <- list(
    list(
      family = "gaussian",
      coefficients = c("(Intercept)" = 0.0, x = 2.0),
      sigma2 = 1.0
    )
  )
  
  result <- compute_score_matrix(data, model_list, "y")
  
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
})

test_that("compute_information_matrix combines information correctly", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(2, 4, 6),
    z = c(0, 1, 1),
    .imputed = c(0, 0, 1)
  )
  
  model_list <- list(
    list(
      family = "gaussian",
      coefficients = c("(Intercept)" = 0.0, x = 2.0),
      sigma2 = 1.0
    ),
    list(
      family = "binomial",
      coefficients = c("(Intercept)" = 0.0, y = 0.5)
    )
  )
  
  result <- compute_information_matrix(data, model_list, c("y", "z"))
  
  expect_equal(dim(result), c(5, 5))
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("compute_imputation_components returns correct structure", {
  data.i <- data.frame(
    x = c(1, 2, 3, 4),
    y = c(2, 4, 7, 3),
    z = c(0, 1, 1, 0),
    .imputed_z = c(0, 1, 1, 0),
    .imputed_y = c(0, 1, 0, 1)
  )
  
  model_list <- list(
    list(
      family = "gaussian",
      coefficients = c("(Intercept)" = 0.0, x = 2.0),
      sigma2 = 1.0
    ),
    list(
      family = "binomial",
      coefficients = c("(Intercept)" = 0.0, y = 0.5)
    )
  )
  
  result <- compute_imputation_components(data.i, model_list, c("y", "z"))
  
  expect_named(result, c("S_mis_imp", "d"))
  expect_equal(nrow(result$S_mis_imp), 4)
  expect_equal(nrow(result$d), 4)
})

test_that("compute_imputation_components masks non-imputed observations", {
  data.i <- data.frame(
    x = c(1, 2, 4),
    y = c(2, 4, 6),
    .imputed_y = c(0, 0, 0) 
  )
  
  model_list <- list(
    list(
      family = "gaussian",
      coefficients = c("(Intercept)" = 0.0, x = 2.0),
      sigma2 = 1.0
    )
  )
  
  result <- compute_imputation_components(data.i, model_list, "y")
  
  expect_true(all(result$S_mis_imp == 0))
})

test_that("pmmrw contributes its stored smooth score and d", {
  data.i <- data.frame(
    x = c(1, 2, 3, 4),
    y = c(2, 4, 7, 3),
    .imputed_y = c(0, 1, 0, 1)
  )

  pmm_score <- matrix(1:8, nrow = 4)
  pmm_d <- matrix(9:16, nrow = 4)

  model_list <- list(
    list(
      method = "pmmrw",
      family = "gaussian",
      coefficients = c("(Intercept)" = 0.0, x = 2.0),
      pmm_score = pmm_score,
      pmm_d = pmm_d
    )
  )

  result <- compute_imputation_components(data.i, model_list, "y")

  expect_equal(result$S_mis_imp, pmm_score)
  expect_equal(result$d, pmm_d)
})


test_that("compute_score handles no intercept model", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(2, 4, 6)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c(x = 2.0),
    sigma2 = 1.0
  )
  
  result <- compute_score(data, model_info, "y")
  
  expect_equal(ncol(result), 2)
})

test_that("compute_information handles no intercept model", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(2, 4, 6)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c(x = 2.0),
    sigma2 = 1.0
  )
  
  imputed_flag <- rep(0, 3)
  
  result <- compute_information(data, model_info, "y", imputed_flag)
  
  expect_equal(dim(result), c(2, 2))
})

test_that("compute_score gaussian sigma column is correct", {
  data <- data.frame(
    x = c(1, 2),
    y = c(3, 5)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c("(Intercept)" = 1.0, x = 2.0),
    sigma2 = 4.0
  )
  
  result <- compute_score(data, model_info, "y")
  
  pred <- 1.0 + 2.0 * data$x
  resid <- data$y - pred
  expected_sigma_score <- 0.5 * (-1 / 4.0 + resid^2 / 16.0)
  
  expect_equal(result[, 3], expected_sigma_score)
})

test_that("build_block_diagonal preserves matrix values", {
  set.seed(123)
  mat1 <- matrix(rnorm(9), 3, 3)
  mat2 <- matrix(rnorm(4), 2, 2)
  
  result <- build_block_diagonal(list(mat1, mat2))
  
  expect_equal(result[1:3, 1:3], mat1)
  expect_equal(result[4:5, 4:5], mat2)
})

test_that("compute_score_matrix preserves row count", {
  n <- 10
  data <- data.frame(
    x = rnorm(n),
    y = rnorm(n),
    z = sample(0:1, n, replace = TRUE)
  )
  
  model_list <- list(
    list(family = "gaussian", coefficients = c("(Intercept)" = 0, x = 1), sigma2 = 1),
    list(family = "binomial", coefficients = c("(Intercept)" = 0, y = 0.5))
  )
  
  result <- compute_score_matrix(data, model_list, c("y", "z"))
  
  expect_equal(nrow(result), n)
})

test_that("compute_information_matrix is always symmetric", {
  data <- data.frame(
    x = rnorm(15),
    y = rnorm(15),
    z = sample(0:1, 15, replace = TRUE),
    .imputed = sample(0:1, 15, replace = TRUE)
  )
  
  model_list <- list(
    list(family = "gaussian", coefficients = c("(Intercept)" = 0, x = 1), sigma2 = 1),
    list(family = "binomial", coefficients = c("(Intercept)" = 0, y = 0.5))
  )
  
  result <- compute_information_matrix(data, model_list, c("y", "z"))
  
  expect_equal(result, t(result), tolerance = 1e-10)
})
