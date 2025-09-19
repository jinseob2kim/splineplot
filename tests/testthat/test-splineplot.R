test_that("splineplot works with GAM Cox models", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("survival")

  # Generate test data
  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  lp <- -0.06*(x - 35) + 0.0009*(x - 35)^3/(6^2)
  time <- rexp(n, rate = exp(lp))
  status <- rbinom(n, 1, 0.8)
  dat <- data.frame(x, time, status)

  # Fit GAM Cox model using weights method
  fit <- mgcv::gam(time ~ s(x),
                    family = mgcv::cox.ph(), weights = status, data = dat)

  # Test basic plot creation
  p <- splineplot(fit, dat)
  expect_s3_class(p, "ggplot")

  # Test with custom parameters
  p2 <- splineplot(fit, dat, refx = 35, ylim = c(0.5, 2.0))
  expect_s3_class(p2, "ggplot")

  # Test with weights method (preferred for GAM Cox)
  fit2 <- mgcv::gam(time ~ s(x),
                     family = mgcv::cox.ph(),
                     weights = status,
                     data = dat)
  p3 <- splineplot(fit2, dat)
  expect_s3_class(p3, "ggplot")
})

test_that("splineplot works with GAM logistic models", {
  skip_if_not_installed("mgcv")

  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  lp <- -0.06*(x - 35) + 0.0009*(x - 35)^3/(6^2)
  y <- rbinom(n, 1, plogis(lp))
  dat <- data.frame(x, y)

  # Fit GAM logistic
  fit <- mgcv::gam(y ~ s(x), family = binomial(), data = dat)
  p <- splineplot(fit, dat)
  expect_s3_class(p, "ggplot")

  # Test log scale
  p2 <- splineplot(fit, dat, log_scale = TRUE)
  expect_s3_class(p2, "ggplot")
})

test_that("splineplot works with GAM poisson models", {
  skip_if_not_installed("mgcv")

  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  lp <- -0.06*(x - 35) + 0.0009*(x - 35)^3/(6^2)
  y <- rpois(n, lambda = exp(lp/2))
  dat <- data.frame(x, y)

  fit <- mgcv::gam(y ~ s(x), family = poisson(), data = dat)
  p <- splineplot(fit, dat)
  expect_s3_class(p, "ggplot")
})

test_that("splineplot works with GAM gaussian models", {
  skip_if_not_installed("mgcv")

  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  y <- 2 + 0.1*x + rnorm(n, 0, 0.5)
  dat <- data.frame(x, y)

  fit <- mgcv::gam(y ~ s(x), data = dat)
  p <- splineplot(fit, dat)
  expect_s3_class(p, "ggplot")
})

test_that("splineplot works with GLM models using ns()", {
  skip_if_not_installed("splines")

  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  lp <- -0.06*(x - 35) + 0.0009*(x - 35)^3/(6^2)
  y <- rbinom(n, 1, plogis(lp))
  dat <- data.frame(x, y)

  # GLM with ns()
  fit <- glm(y ~ splines::ns(x, df = 4), family = binomial(), data = dat)
  p <- splineplot(fit, dat)
  expect_s3_class(p, "ggplot")

  # Check reference point
  df <- extract_spline_data(fit, dat, "x", refx = 35,
                            list(type = "glm", family = "binomial", ylabel = "OR"),
                            log_scale = FALSE, ci_level = 0.95)
  ref_idx <- which.min(abs(df$x - 35))
  expect_true(abs(df$y[ref_idx] - 1) < 0.1)  # OR should be close to 1
})

test_that("splineplot works with GLM models using bs()", {
  skip_if_not_installed("splines")

  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  lp <- -0.06*(x - 35)
  y <- rpois(n, lambda = exp(lp))
  dat <- data.frame(x, y)

  fit <- glm(y ~ splines::bs(x, df = 4), family = poisson(), data = dat)
  p <- splineplot(fit, dat)
  expect_s3_class(p, "ggplot")
})

test_that("splineplot works with linear models", {
  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  y <- 2 + 0.1*x + rnorm(n, 0, 0.5)
  dat <- data.frame(x, y)

  fit <- lm(y ~ splines::ns(x, df = 3), data = dat)
  p <- splineplot(fit, dat, xvar = "x")  # Need to specify xvar for lm
  expect_s3_class(p, "ggplot")
})

test_that("splineplot works with Cox models", {
  skip_if_not_installed("survival")
  skip_if_not_installed("splines")

  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  time <- rexp(n, rate = 0.1)
  status <- rbinom(n, 1, 0.8)
  dat <- data.frame(time, status, x)

  # Cox with ns()
  fit <- survival::coxph(survival::Surv(time, status) ~ splines::ns(x, df = 4),
                         data = dat)
  p <- splineplot(fit, dat)
  expect_s3_class(p, "ggplot")

  # Cox with bs()
  fit2 <- survival::coxph(survival::Surv(time, status) ~ splines::bs(x, df = 4),
                          data = dat)
  p2 <- splineplot(fit2, dat)
  expect_s3_class(p2, "ggplot")
})

test_that("splineplot handles interaction terms", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("survival")

  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  group <- factor(sample(c("A", "B"), n, replace = TRUE))
  time <- rexp(n, rate = 0.1)
  status <- rbinom(n, 1, 0.8)
  dat <- data.frame(x, group, time, status)

  # GAM with interaction
  fit <- mgcv::gam(time ~ s(x, by = group),
                    family = mgcv::cox.ph(), weights = status, data = dat)
  p <- splineplot(fit, dat)
  expect_s3_class(p, "ggplot")
})

test_that("splineplot options work correctly", {
  skip_if_not_installed("mgcv")

  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  y <- rbinom(n, 1, 0.5)
  dat <- data.frame(x, y)

  fit <- mgcv::gam(y ~ s(x), family = binomial(), data = dat)

  # Test ribbon_ci option
  p1 <- splineplot(fit, dat, ribbon_ci = FALSE)
  expect_s3_class(p1, "ggplot")

  p2 <- splineplot(fit, dat, ribbon_ci = TRUE)
  expect_s3_class(p2, "ggplot")

  # Test show_hist option
  p3 <- splineplot(fit, dat, show_hist = FALSE)
  expect_s3_class(p3, "ggplot")

  # Test show_ref_point option
  p4 <- splineplot(fit, dat, show_ref_point = FALSE)
  expect_s3_class(p4, "ggplot")

  # Test custom labels
  p5 <- splineplot(fit, dat,
                   xlab = "Age",
                   ylab = "Odds Ratio (95% CI)",
                   ylab_right = "Distribution (%)")
  expect_s3_class(p5, "ggplot")

  # Test xlim and ylim
  p6 <- splineplot(fit, dat,
                   xlim = c(20, 50),
                   ylim = c(0.5, 2.0))
  expect_s3_class(p6, "ggplot")
})

test_that("extract_spline_data returns correct structure", {
  skip_if_not_installed("mgcv")

  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  y <- rbinom(n, 1, 0.5)
  dat <- data.frame(x, y)

  fit <- mgcv::gam(y ~ s(x), family = binomial(), data = dat)

  df <- extract_spline_data(fit, dat, "x", refx = 35,
                            list(type = "gam", family = "binomial", ylabel = "OR"),
                            log_scale = FALSE, ci_level = 0.95)

  expect_true(is.data.frame(df))
  expect_true(all(c("x", "y", "lcl", "ucl", "ylabel") %in% names(df)))
  expect_equal(nrow(df), 200)  # Default grid size

  # Check reference point
  ref_idx <- which.min(abs(df$x - 35))
  expect_true(abs(df$y[ref_idx] - 1) < 0.1)  # Should be close to 1 for OR
})

test_that("detect functions work correctly", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("splines")

  # Test with GAM
  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  y <- rbinom(n, 1, 0.5)
  dat <- data.frame(x, y)

  fit_gam <- mgcv::gam(y ~ s(x), family = binomial(), data = dat)

  # These are internal functions, so we need to use :::
  model_info <- splineplot:::detect_model_info(fit_gam)
  expect_equal(model_info$type, "gam")
  expect_equal(model_info$family, "binomial")
  expect_equal(model_info$ylabel, "Odds Ratio")

  spline_info <- splineplot:::detect_spline_terms(fit_gam, model_info)
  expect_equal(spline_info$xvar, "x")
  expect_null(spline_info$by_var)

  # Test with GLM
  fit_glm <- glm(y ~ splines::ns(x, df = 4), family = binomial(), data = dat)

  model_info2 <- splineplot:::detect_model_info(fit_glm)
  expect_equal(model_info2$type, "glm")
  expect_equal(model_info2$family, "binomial")

  spline_info2 <- splineplot:::detect_spline_terms(fit_glm, model_info2)
  expect_equal(spline_info2$xvar, "x")
})

test_that("error handling works correctly", {
  # Test with unsupported model
  dat <- data.frame(x = 1:10, y = 1:10)
  fit <- list(class = "unknown")

  expect_error(splineplot(fit, dat),
              "Unsupported model type")

  # Test with missing spline terms
  fit2 <- lm(y ~ x, data = dat)
  expect_error(splineplot(fit2, dat),
              "No spline terms found")
})

test_that("splineplot works with models containing covariates", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("survival")
  skip_if_not_installed("splines")

  set.seed(123)
  n <- 200
  x <- rnorm(n, mean = 50, sd = 10)
  sex <- factor(sample(c("Male", "Female"), n, replace = TRUE))
  age_group <- factor(sample(c("Young", "Middle", "Old"), n, replace = TRUE))
  bmi <- rnorm(n, mean = 25, sd = 3)

  # Create outcomes
  lp <- -0.05*(x - 50) + 0.001*(x - 50)^2 +
        ifelse(sex == "Male", 0.3, 0) +
        0.02*bmi

  time <- rexp(n, rate = exp(lp/2))
  status <- rbinom(n, 1, 0.8)
  binary <- rbinom(n, 1, plogis(lp))
  count <- rpois(n, lambda = exp(lp/3))

  dat <- data.frame(x, sex, age_group, bmi, time, status, binary, count)

  # GAM Cox with covariates
  fit_gam_cox <- mgcv::gam(time ~ s(x) + sex + bmi,
                           family = mgcv::cox.ph(),
                           weights = status,
                           data = dat)
  p1 <- splineplot(fit_gam_cox, dat)
  expect_s3_class(p1, "ggplot")

  # GAM logistic with covariates
  fit_gam_logit <- mgcv::gam(binary ~ s(x) + sex + age_group + bmi,
                             family = binomial(),
                             data = dat)
  p2 <- splineplot(fit_gam_logit, dat)
  expect_s3_class(p2, "ggplot")

  # GLM with ns() and covariates
  fit_glm_ns <- glm(binary ~ ns(x, df = 4) + sex + bmi,
                    family = binomial(),
                    data = dat)
  p3 <- splineplot(fit_glm_ns, dat)
  expect_s3_class(p3, "ggplot")

  # GLM with bs() and covariates
  fit_glm_bs <- glm(count ~ bs(x, df = 4) + sex + age_group,
                    family = poisson(),
                    data = dat)
  p4 <- splineplot(fit_glm_bs, dat)
  expect_s3_class(p4, "ggplot")

  # Cox with ns() and covariates
  fit_cox_ns <- survival::coxph(Surv(time, status) ~ ns(x, df = 4) + sex + bmi,
                                data = dat)
  p5 <- splineplot(fit_cox_ns, dat)
  expect_s3_class(p5, "ggplot")

  # Cox with bs() and covariates
  fit_cox_bs <- survival::coxph(Surv(time, status) ~ bs(x, df = 3) + sex + age_group + bmi,
                                data = dat)
  p6 <- splineplot(fit_cox_bs, dat)
  expect_s3_class(p6, "ggplot")

  # Linear model with splines and covariates
  fit_lm <- lm(bmi ~ ns(x, df = 3) + sex,
               data = dat)
  p7 <- splineplot(fit_lm, dat, xvar = "x")
  expect_s3_class(p7, "ggplot")

  # GAM with interaction and covariates
  fit_gam_interact <- mgcv::gam(time ~ s(x, by = sex) + age_group + bmi,
                                family = mgcv::cox.ph(),
                                weights = status,
                                data = dat)
  p8 <- splineplot(fit_gam_interact, dat)
  expect_s3_class(p8, "ggplot")
})

test_that("splineplot handles multiple spline terms correctly", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("splines")
  skip_if_not_installed("survival")

  set.seed(123)
  n <- 200
  age <- rnorm(n, mean = 50, sd = 10)
  bmi <- rnorm(n, mean = 25, sd = 3)
  sbp <- rnorm(n, mean = 120, sd = 15)

  # Create outcomes
  lp <- -0.05*(age - 50) + 0.001*(age - 50)^2 +
        -0.02*(bmi - 25) + 0.002*(bmi - 25)^2 +
        0.01*sbp

  binary <- rbinom(n, 1, plogis(lp))
  time <- rexp(n, rate = exp(lp/3))
  status <- rbinom(n, 1, 0.8)

  dat <- data.frame(age, bmi, sbp, binary, time, status)

  # GAM with multiple smooth terms - should use first one (age) by default
  fit_gam_multi <- mgcv::gam(binary ~ s(age) + s(bmi) + sbp,
                             family = binomial(),
                             data = dat)
  p1 <- splineplot(fit_gam_multi, dat)
  expect_s3_class(p1, "ggplot")

  # Should be able to specify which spline to plot
  p2 <- splineplot(fit_gam_multi, dat, xvar = "bmi")
  expect_s3_class(p2, "ggplot")

  # GLM with multiple spline terms - should use first one by default
  fit_glm_multi <- glm(binary ~ ns(age, df = 4) + ns(bmi, df = 3) + sbp,
                      family = binomial(),
                      data = dat)
  p3 <- splineplot(fit_glm_multi, dat)
  expect_s3_class(p3, "ggplot")

  # Should be able to specify which spline to plot
  p4 <- splineplot(fit_glm_multi, dat, xvar = "bmi")
  expect_s3_class(p4, "ggplot")

  # Cox with multiple spline terms
  fit_cox_multi <- survival::coxph(Surv(time, status) ~ ns(age, df = 4) + bs(bmi, df = 3) + sbp,
                                   data = dat)
  p5 <- splineplot(fit_cox_multi, dat)
  expect_s3_class(p5, "ggplot")

  # Should be able to specify which spline to plot
  p6 <- splineplot(fit_cox_multi, dat, xvar = "bmi")
  expect_s3_class(p6, "ggplot")

  # Mixed spline and interaction
  dat$sex <- factor(sample(c("M", "F"), n, replace = TRUE))
  fit_gam_mixed <- mgcv::gam(binary ~ s(age, by = sex) + s(bmi) + sbp,
                             family = binomial(),
                             data = dat)
  # Should detect interaction with first spline
  p7 <- splineplot(fit_gam_mixed, dat)
  expect_s3_class(p7, "ggplot")

  # Should be able to plot the non-interaction spline
  p8 <- splineplot(fit_gam_mixed, dat, xvar = "bmi")
  expect_s3_class(p8, "ggplot")
})

test_that("Cox pspline has limited but functional support", {
  skip_if_not_installed("survival")

  set.seed(123)
  n <- 100
  x <- rnorm(n, mean = 35, sd = 6)
  time <- rexp(n, rate = 0.1)
  status <- rbinom(n, 1, 0.8)
  dat <- data.frame(time, status, x)

  # Cox with pspline - this should work but with limitations
  fit <- survival::coxph(survival::Surv(time, status) ~
                         survival::pspline(x, df = 4), data = dat)

  # Should not error, but SE at reference won't be exactly 0
  expect_no_error({
    p <- splineplot(fit, dat, xvar = "x")
  })
  expect_s3_class(p, "ggplot")
})