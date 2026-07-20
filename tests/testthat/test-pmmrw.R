test_that("pmmrw records donor IDs during mice imputation", {
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
    m = 2,
    tasks = "train",
    print = FALSE,
    seed = 1
  )
  imp_standard <- mice::mice(
    data,
    method = method_standard,
    predictorMatrix = pred,
    m = 2,
    tasks = "train",
    print = FALSE,
    seed = 1
  )

  donor_id <- extract_donor_id(imp, "bmi")
  model <- imp$models$bmi[[1]]

  expect_equal(dim(donor_id), c(nrow(data), 2))
  expect_true(all(is.na(donor_id[!is.na(data$bmi), ])))
  expect_true(all(!is.na(donor_id[is.na(data$bmi), ])))
  expect_true(all(donor_id[is.na(data$bmi), ] %in% which(!is.na(data$bmi))))
  expect_identical(imp$imp$bmi, imp_standard$imp$bmi)
  expect_equal(dim(model$pmm_score), c(nrow(data), 2))
  expect_equal(dim(model$pmm_d), c(nrow(data), 2))
  expect_true(all(is.finite(model$pmm_score)))
  expect_true(any(model$pmm_score != 0))
  expect_true(any(model$pmm_d != 0))
  expect_named(
    model$pmm_diagnostics,
    c("mean_base_width", "mean_donor_width", "mean_working_width")
  )
  expect_true(all(model$pmm_diagnostics > 0))
  expect_gte(
    model$pmm_diagnostics[["mean_working_width"]],
    model$pmm_diagnostics[["mean_base_width"]]
  )
})

test_that("PMM neighbor draws retain their MICE selection-law state", {
  yhat <- c(0.25, 1.25)
  edges <- c(0, 1, 2, 3)
  lookup <- matrix(1:9, nrow = 3, byrow = TRUE)
  donor_lookup <- matrix(11:19, nrow = 3, byrow = TRUE)

  set.seed(1)
  draws <- draw_neighbors_pmmrw(yhat, edges, lookup, donor_lookup)

  expect_equal(draws$bin, c(1L, 2L))
  expect_equal(draws$p_left, c(0.75, 0.75))
  expect_true(all(draws$selected_bin %in% c(1L, 2L, 3L)))
  expect_true(all(draws$donor_id %in% donor_lookup))
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

test_that("base smooth scale uses the three-quarter width rule", {
  eta_obs <- c(-1, 0, 2, 3)
  eta_mis <- c(0.25, 1.5)
  donors <- 5

  result <- smooth_pmm_working_scale(eta_obs, eta_mis, donors)

  delta <- outer(eta_obs, eta_mis, "-")
  dist2 <- delta^2
  scale <- min(sd(eta_obs), IQR(eta_obs) / 1.34)
  h_rule <- 0.9 * scale * length(eta_obs)^(-1 / 5)
  lambda0 <- rep(0.5 / h_rule^2, length(eta_mis))
  log_w0 <- -sweep(dist2, 2, lambda0, "*")
  log_w0 <- sweep(log_w0, 2, apply(log_w0, 2, max), "-")
  w0 <- sweep(exp(log_w0), 2, colSums(exp(log_w0)), "/")
  neff0 <- 1 / colSums(w0^2)
  h <- h_rule * pmax(donors / neff0, 1e-4)^0.75

  expect_equal(result$lambda, 0.5 / h^2)
})
