test_that("pmmrw reproduces ordinary MICE PMM and records its donors", {
  skip_if_not_installed("mice")

  data <- mice::nhanes
  method <- mice::make.method(data)
  method[] <- ""
  method["bmi"] <- "pmmrw"
  method_standard <- method
  method_standard["bmi"] <- "pmm"

  pred <- mice::make.predictorMatrix(data)
  pred[,] <- 0
  pred["bmi", "age"] <- 1

  imp <- mice::mice(
    data,
    method = method,
    predictorMatrix = pred,
    m = 5,
    tasks = "train",
    print = FALSE,
    seed = 1
  )
  imp_standard <- mice::mice(
    data,
    method = method_standard,
    predictorMatrix = pred,
    m = 5,
    print = FALSE,
    seed = 1
  )

  donor_id <- extract_donor_id(imp, "bmi")
  missing <- which(is.na(data$bmi))
  observed <- which(!is.na(data$bmi))
  model <- imp$models$bmi[[1]]

  expect_identical(imp$imp$bmi, imp_standard$imp$bmi)
  expect_equal(dim(donor_id), c(nrow(data), 5))
  expect_true(all(is.na(donor_id[observed, ])))
  expect_true(all(donor_id[missing, ] %in% observed))
  for (p in seq_len(imp$m)) {
    expect_equal(
      as.numeric(imp$imp$bmi[, p]),
      data$bmi[donor_id[missing, p]]
    )
  }
  expect_equal(dim(model$pmm_score), c(nrow(data), 2))
  expect_equal(dim(model$pmm_d), c(nrow(data), 2))
  expect_true(all(is.finite(model$pmm_score)))
  expect_true(any(model$pmm_score != 0))
  expect_true(any(model$pmm_d != 0))
  expect_named(
    model$pmm_diagnostics,
    c("mean_bandwidth", "mean_effective_donors")
  )
  expect_true(all(model$pmm_diagnostics > 0))
})

test_that("stored PMM score equals the all-donor working-score derivative", {
  skip_if_not_installed("mice")

  data <- mice::nhanes[c("age", "bmi")]
  method <- mice::make.method(data)
  method[] <- ""
  method["bmi"] <- "pmmrw"
  pred <- mice::make.predictorMatrix(data)
  pred[,] <- 0
  pred["bmi", "age"] <- 1

  imp <- mice::mice(
    data,
    method = method,
    predictorMatrix = pred,
    m = 1,
    tasks = "train",
    print = FALSE,
    seed = 12
  )

  model <- imp$models$bmi[[1]]
  observed <- which(!is.na(data$bmi))
  missing <- which(is.na(data$bmi))
  x <- cbind(`(Intercept)` = 1, age = data$age)
  x_obs <- x[observed, , drop = FALSE]
  x_mis <- x[missing, , drop = FALSE]
  eta_obs <- as.vector(x_obs %*% model$beta.work)
  eta_mis <- as.vector(x_mis %*% model$beta.work)
  working <- smooth_pmm_working_law(
    eta_obs,
    eta_mis,
    model$setup$donors
  )
  donor_pos <- match(model$donor_id[missing], observed)
  picked <- cbind(donor_pos, seq_along(missing))
  two_delta_weight <- 2 * working$delta * working$weight
  mean_derivative <- crossprod(two_delta_weight, x_obs) -
    sweep(x_mis, 1, colSums(two_delta_weight), "*")
  realized_derivative <- sweep(
    x_obs[donor_pos, , drop = FALSE] - x_mis,
    1,
    2 * working$delta[picked],
    "*"
  )
  expected <- sweep(
    realized_derivative - mean_derivative,
    1,
    -working$lambda,
    "*"
  )

  expect_equal(colSums(working$weight), rep(1, length(missing)))
  expect_equal(model$pmm_score[missing, ], expected)
})

test_that("pool_rw automatically uses pmmrw donor IDs", {
  skip_if_not_installed("mice")

  data <- mice::nhanes
  method <- mice::make.method(data)
  method[] <- ""
  method["bmi"] <- "pmmrw"

  pred <- mice::make.predictorMatrix(data)
  pred[,] <- 0
  pred["bmi", "age"] <- 1

  imp <- mice::mice(
    data,
    method = method,
    predictorMatrix = pred,
    m = 2,
    tasks = "train",
    print = FALSE,
    seed = 1
  )

  fit <- with_rw(imp, lm(bmi ~ age))
  pooled <- pool_rw(fit)

  expect_true(pooled$donor_correlation)
  expect_s3_class(pooled, "pool_rw")
  expect_true(any(fit$results[[1]]$S_mis_imp != 0))
  expect_true(any(fit$results[[1]]$d != 0))
})

test_that("smooth working law uses the three-quarter baseline width", {
  eta_obs <- c(-1, 0, 2, 3)
  eta_mis <- c(0.25, 1.5)
  donors <- 4

  result <- smooth_pmm_working_law(eta_obs, eta_mis, donors)

  delta <- outer(eta_obs, eta_mis, "-")
  dist2 <- delta^2
  scale <- min(sd(eta_obs), IQR(eta_obs) / 1.34)
  h_rule <- 0.9 * scale * length(eta_obs)^(-1 / 5)
  lambda0 <- rep(0.5 / h_rule^2, length(eta_mis))
  log_weight0 <- -sweep(dist2, 2, lambda0, "*")
  log_weight0 <- sweep(log_weight0, 2, apply(log_weight0, 2, max), "-")
  weight0 <- exp(log_weight0)
  weight0 <- sweep(weight0, 2, colSums(weight0), "/")
  effective0 <- 1 / colSums(weight0^2)
  bandwidth <- h_rule * pmax(donors / effective0, 1e-4)^0.75

  expect_equal(result$bandwidth, bandwidth)
  expect_equal(result$lambda, 0.5 / bandwidth^2)
  expect_equal(colSums(result$weight), rep(1, length(eta_mis)))
})
