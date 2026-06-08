test_that("extract_sigma2 handles different sigma formats", {
  mod1 <- list(sigma2 = 4.5)
  expect_equal(extract_sigma2(mod1, "bmi"), 4.5)
  
  mod2 <- list(sigma.dot = 2.0)
  expect_equal(extract_sigma2(mod2, "bmi"), 4.0)
  
  mod3 <- list(sigma = 3.0)
  expect_equal(extract_sigma2(mod3, "bmi"), 9.0)
  
  mod4 <- list(beta.dot = c(1, 2))
  expect_error(
    extract_sigma2(mod4, "bmi"),
    "Cannot extract.*sigma2"
  )
})

test_that("extract_single_model handles logreg correctly", {
  data <- list(
    models = list(
      hyp = list(list(
        beta.dot = c(0.5, -1.2),
        xnames = c("(Intercept)", "age")
      ))
    ),
    method = c(hyp = "logreg")
  )
  
  result <- extract_single_model(data, "hyp", 1)
  
  expect_equal(result$family, "binomial")
  expect_named(result$coefficients, c("(Intercept)", "age"))
  expect_equal(result$coefficients[[1]], 0.5)
  expect_equal(result$coefficients[[2]], -1.2)
})

test_that("extract_single_model handles norm correctly", {
  data <- list(
    models = list(
      bmi = list(list(
        beta.dot = c(25, 2.5),
        xnames = c("(Intercept)", "age"),
        sigma2 = 16.0
      ))
    ),
    method = c(bmi = "norm")
  )
  
  result <- extract_single_model(data, "bmi", 1)
  
  expect_equal(result$family, "gaussian")
  expect_equal(result$sigma2, 16.0)
  expect_named(result$coefficients, c("(Intercept)", "age"))
})

test_that("extract_single_model rejects unsupported methods", {
  data <- list(
    models = list(
      var = list(list(beta.dot = 1, xnames = "x"))
    ),
    method = c(var = "pmm")
  )
  
  expect_error(
    extract_single_model(data, "var", 1),
    "Unsupported imputation method.*pmm"
  )
})

test_that("extract_imputation_models processes multiple imputations", {
  data <- list(
    models = list(
      bmi = list(
        list(beta.dot = c(25, 2), xnames = c("int", "age"), sigma2 = 16),
        list(beta.dot = c(26, 2.1), xnames = c("int", "age"), sigma2 = 17)
      )
    ),
    method = c(bmi = "norm")
  )
  
  result <- extract_imputation_models(data, "bmi", 2)
  
  expect_length(result, 2)
  expect_named(result[[1]], "bmi")
  expect_named(result[[2]], "bmi")
  expect_equal(result[[1]]$bmi$sigma2, 16)
  expect_equal(result[[2]]$bmi$sigma2, 17)
})

test_that("extract_imputation_models handles multiple variables", {
  data <- list(
    models = list(
      bmi = list(
        list(beta.dot = c(25, 2), xnames = c("int", "age"), sigma2 = 16)
      ),
      hyp = list(
        list(beta.dot = c(0.5, -1.2), xnames = c("int", "bmi"))
      )
    ),
    method = c(bmi = "norm", hyp = "logreg")
  )
  
  result <- extract_imputation_models(data, c("bmi", "hyp"), 1)
  
  expect_length(result, 1)
  expect_named(result[[1]], c("bmi", "hyp"))
  expect_equal(result[[1]]$bmi$family, "gaussian")
  expect_equal(result[[1]]$hyp$family, "binomial")
})

test_that("convert_logreg_factors handles 0/1 factor levels", {
  data <- data.frame(
    x = factor(c("0", "1", "0", "1"), levels = c("0", "1")),
    y = 1:4
  )
  
  result <- convert_logreg_factors(data, "x")
  
  expect_type(result$x, "double")
  expect_equal(result$x, c(0, 1, 0, 1))
})

test_that("convert_logreg_factors handles generic binary factors", {
  data <- data.frame(
    x = factor(c("no", "yes", "no", "yes"), levels = c("no", "yes")),
    y = 1:4
  )
  
  result <- convert_logreg_factors(data, "x")
  
  expect_type(result$x, "double")
  expect_equal(result$x, c(0, 1, 0, 1))
})

test_that("convert_logreg_factors leaves non-factors unchanged", {
  data <- data.frame(
    x = c(1.5, 2.5, 3.5),
    y = 1:3
  )
  
  result <- convert_logreg_factors(data, "x")
  
  expect_equal(result$x, c(1.5, 2.5, 3.5))
})

test_that("convert_logreg_factors handles multiple variables", {
  data <- data.frame(
    x = factor(c("0", "1", "0"), levels = c("0", "1")),
    y = factor(c("a", "b", "a"), levels = c("a", "b")),
    z = c(1.5, 2.5, 3.5)
  )
  
  result <- convert_logreg_factors(data, c("x", "y"))
  
  expect_type(result$x, "double")
  expect_type(result$y, "double")
  expect_equal(result$x, c(0, 1, 0))
  expect_equal(result$y, c(0, 1, 0))
  expect_equal(result$z, c(1.5, 2.5, 3.5))
})

test_that("fit_analysis_model fits lm correctly", {
  data <- data.frame(
    y = c(2, 4, 6, 8),
    x = c(1, 2, 3, 4)
  )
  
  expr <- quote(lm(y ~ x))
  result <- fit_analysis_model(expr, data, parent.frame())
  
  expect_s3_class(result, "lm")
  expect_equal(coef(result)[[2]], 2, tolerance = 1e-10)
})

test_that("fit_analysis_model fits glm correctly", {
  data <- data.frame(
    y = c(0, 1, 1, 0, 1),
    x = c(1, 2, 3, 1, 0)
  )
  
  expr <- quote(glm(y ~ x, family = binomial()))
  result <- fit_analysis_model(expr, data, parent.frame())
  
  expect_s3_class(result, "glm")
  expect_s3_class(result, "lm")
})

test_that("fit_analysis_model rejects non-lm/glm models", {
  data <- data.frame(x = 1:10, y = 1:10)
  
  expr <- quote(mean(y))
  
  expect_error(
    fit_analysis_model(expr, data, parent.frame()),
    "Analysis model must be.*lm.*or.*glm"
  )
})

test_that("compute_lm_components calculates U and tau correctly", {
  data <- data.frame(y = c(2.1, 3.8, 6.2, 7.7), x = c(1, 2, 3, 4))
  mod <- lm(y ~ x, data = data)
  X <- model.matrix(mod)
  
  result <- compute_lm_components(mod, X)
  sigma2 <- summary(mod)$sigma^2
  
  expect_named(result, c("U", "tau"))
  expect_equal(nrow(result$U), 4)
  expect_equal(ncol(result$U), 2)
  expect_equal(dim(result$tau), c(2, 2))
  expect_equal(result$tau, -crossprod(X) / sigma2)
})

test_that("compute_lm_components U has correct residual structure", {
  data <- data.frame(y = c(1.2, 2.7, 5.1, 6.8), x = c(1, 2, 3, 4))
  mod <- lm(y ~ x, data = data)
  X <- model.matrix(mod)
  
  result <- compute_lm_components(mod, X)
  
  resid <- residuals(mod)
  sigma2 <- summary(mod)$sigma^2
  expected_U <- X * resid / sigma2
  
  expect_equal(result$U, expected_U)
})

test_that("compute_glm_components calculates U and tau correctly", {
  data <- data.frame(
    y = c(0, 1, 1, 0, 1),
    x = c(1, 2, 3, 1, 0)
  )
  mod <- glm(y ~ x, data = data, family = binomial())
  X <- model.matrix(mod)
  
  result <- compute_glm_components(mod, X)
  
  expect_named(result, c("U", "tau"))
  expect_equal(nrow(result$U), 5)
  expect_equal(ncol(result$U), 2)
  expect_equal(dim(result$tau), c(2, 2))
  expect_equal(result$tau, t(result$tau))
})

test_that("compute_glm_components U reflects residuals", {
  data <- data.frame(
    y = c(0, 1, 1, 0, 1),
    x = c(1, 2, 3, 1, 0)
  )
  mod <- glm(y ~ x, data = data, family = binomial())
  X <- model.matrix(mod)
  
  result <- compute_glm_components(mod, X)
  
  y <- model.response(model.frame(mod))
  pred <- fitted(mod)
  expected_U <- X * (y - pred)
  
  expect_equal(result$U, expected_U)
})

test_that("extend_score_to_full_data pads with zeros correctly", {
  data <- data.frame(
    y = c(2, 4, 6, 8, 10),
    x = c(1, 2, 3, 4, 5)
  )
  mod <- lm(y ~ x, data = data[2:4, ])
  mf <- model.frame(mod)
  U_analysis <- model.matrix(mod) * residuals(mod)
  
  result <- extend_score_to_full_data(U_analysis, mf, data, 5, mod)
  
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 2)
  expect_equal(unname(result[1, ]), c(0, 0))
  expect_equal(unname(result[5, ]), c(0, 0))
  expect_true(any(unname(result[2:4, ]) != 0))
})

test_that("extend_score_to_full_data preserves column names", {
  data <- data.frame(y = 1:5, x = 1:5)
  mod <- lm(y ~ x, data = data)
  mf <- model.frame(mod)
  U_analysis <- model.matrix(mod) * residuals(mod)
  
  result <- extend_score_to_full_data(U_analysis, mf, data, 5, mod)
  
  expect_equal(colnames(result), names(coef(mod)))
})

test_that("extend_score_to_full_data handles missing rownames", {
  data <- data.frame(y = 1:5, x = 1:5)
  mod <- lm(y ~ x, data = data)
  mf <- model.frame(mod)
  rownames(mf) <- NULL
  
  U_analysis <- model.matrix(mod) * residuals(mod)
  
  result <- extend_score_to_full_data(U_analysis, mf, data, 5, mod)
  
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 2)
})

test_that("compute_analysis_components delegates to lm function", {
  data_lm <- data.frame(y = c(2, 4, 6), x = c(1, 2, 3))
  mod_lm <- lm(y ~ x, data = data_lm)
  
  result_lm <- compute_analysis_components(mod_lm, data_lm, 3)
  
  expect_named(result_lm, c("U", "tau", "n_analysis"))
  expect_equal(result_lm$n_analysis, 3)
  expect_equal(nrow(result_lm$U), 3)
  expect_equal(ncol(result_lm$U), 2)
})

test_that("compute_analysis_components delegates to glm function", {
  data_glm <- data.frame(
    y = c(0, 1, 1, 0, 1),
    x = c(1, 2, 3, 1, 0)
  )  
  mod_glm <- glm(y ~ x, data = data_glm, family = binomial())
  
  result_glm <- compute_analysis_components(mod_glm, data_glm, 5)
  
  expect_named(result_glm, c("U", "tau", "n_analysis"))
  expect_equal(result_glm$n_analysis, 5)
  expect_equal(nrow(result_glm$U), 5)
  expect_equal(ncol(result_glm$U), 2)
})

test_that("compute_analysis_components handles models with subset", {
  data <- data.frame(y = 1:10, x = 1:10)
  mod <- lm(y ~ x, data = data, subset = 2:8)
  
  result <- compute_analysis_components(mod, data, 10)
  
  expect_equal(result$n_analysis, 7)
  expect_equal(nrow(result$U), 10)
  expect_equal(sum(result$U[1, ]), 0)
  expect_equal(sum(result$U[9, ]), 0)
  expect_equal(sum(result$U[10, ]), 0)
})
