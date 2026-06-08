test_that("pmmrw records donor IDs during mice imputation", {
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

  donor_id <- extract_donor_id(imp, "bmi")

  expect_equal(dim(donor_id), c(nrow(data), 2))
  expect_true(all(is.na(donor_id[!is.na(data$bmi), ])))
  expect_true(all(!is.na(donor_id[is.na(data$bmi), ])))
  expect_true(all(donor_id[is.na(data$bmi), ] %in% which(!is.na(data$bmi))))
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
})
