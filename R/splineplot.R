# Global variables used in ggplot2 aes() calls
utils::globalVariables(c("x", "y", "group", "lcl", "ucl", "xmin", "xmax", "ymin", "ymax"))

#' Spline Plot for GAM and GLM Models
#'
#' Create ggplot2 visualizations of smooth or spline effects from GAM and GLM models.
#' Supports Linear, Logistic, Poisson, and Cox models with interaction terms.
#' Handles GAM smooth terms (s(), te(), ti()), GLM splines (ns(), bs()), and Cox pspline().
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_rect geom_hline
#' @importFrom ggplot2 annotate coord_cartesian scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 labs theme theme_bw element_blank element_line element_text
#' @importFrom ggplot2 margin sec_axis scale_color_manual scale_fill_manual
#' @importFrom stats coef vcov predict median qnorm setNames terms model.matrix delete.response
#' @importFrom graphics hist
#' @importFrom utils head tail
#'
#' @param fit A fitted model object (gam, glm, lm, coxph)
#' @param data The data frame used to fit the model
#' @param xvar Character string specifying the variable name for x-axis (default: first spline term)
#' @param by_var Character string specifying the interaction variable (default: auto-detect from model)
#' @param refx Reference value for the x variable (default: median)
#' @param term_index For GAM with multiple smooth terms, which term to plot (default: 1)
#' @param bins Number of bins for histogram (default: 12)
#' @param xlim X-axis limits (default: range of x variable)
#' @param ylim Y-axis limits (default: auto-determined, e.g., c(0.25, 2.0) for HR/OR/RR)
#' @param show_hist Logical, whether to show histogram (default: TRUE)
#' @param log_scale Logical, whether to use log scale for OR/RR/HR (default: FALSE)
#' @param ci_level Confidence interval level (default: 0.95)
#' @param show_ref_point Logical, whether to show reference point marker (default: TRUE)
#' @param colors Named vector of colors for by_var levels
#' @param ribbon_ci Logical, whether to use ribbon style for CI (default: FALSE, uses dotted lines)
#' @param xlab Custom x-axis label (default: xvar name)
#' @param ylab Custom y-axis label (default: auto-determined based on model type)
#' @param ylab_right Custom right y-axis label for histogram (default: "Percent of Population")
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' n <- 200
#' x <- rnorm(n, mean = 50, sd = 10)
#' lp <- -0.05*(x - 50) + 0.001*(x - 50)^2
#' y <- rbinom(n, 1, plogis(lp))
#' dat <- data.frame(x = x, y = y)
#'
#' # GLM with natural splines
#' library(splines)
#' fit_glm <- glm(y ~ ns(x, df = 4), family = binomial(), data = dat)
#' p <- splineplot(fit_glm, dat)
#'
#' \donttest{
#' # GAM example (requires mgcv)
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   fit_gam <- mgcv::gam(y ~ s(x), family = binomial(), data = dat)
#'   p2 <- splineplot(fit_gam, dat)
#' }
#'
#' # Cox model example (requires survival)
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   time <- rexp(n, rate = exp(lp/2))
#'   status <- rbinom(n, 1, 0.8)
#'   dat$time <- time
#'   dat$status <- status
#'   fit_cox <- survival::coxph(survival::Surv(time, status) ~ ns(x, df = 4),
#'                               data = dat)
#'   p3 <- splineplot(fit_cox, dat)
#' }
#' }
splineplot <- function(fit, data, xvar = NULL, by_var = NULL, refx = NULL,
                      term_index = 1, bins = 12,
                      xlim = NULL, ylim = NULL, show_hist = NULL,
                      log_scale = FALSE, ci_level = 0.95,
                      show_ref_point = TRUE, colors = NULL, ribbon_ci = FALSE,
                      xlab = NULL, ylab = NULL, ylab_right = "Percent of Population") {

  # Determine model type and family first
  model_info <- detect_model_info(fit)

  # Auto-detect spline terms if not provided
  spline_info <- detect_spline_terms(fit, model_info, xvar)

  if (is.null(xvar)) {
    xvar <- spline_info$xvar
    if (is.null(xvar)) {
      stop("No spline terms found in the model. Please specify xvar explicitly.")
    }
    message("Using '", xvar, "' as x variable")
  }

  if (is.null(by_var) && !is.null(spline_info$by_var)) {
    by_var <- spline_info$by_var
    message("Detected interaction with '", by_var, "'")
  }

  # Input validation
  stopifnot(xvar %in% names(data))
  if (!is.null(by_var)) {
    stopifnot(by_var %in% names(data))
  }

  # Set default refx after xvar is determined
  if (is.null(refx)) {
    refx <- stats::median(data[[xvar]], na.rm = TRUE)
    message("Using refx = ", round(refx, 2), " (median of ", xvar, ")")
  }

  # Auto-determine show_hist if not specified
  if (is.null(show_hist)) {
    show_hist <- TRUE  # Always TRUE by default for both single and interaction
  }

  # Extract spline data based on model type and interaction
  if (!is.null(by_var)) {
    # Handle interaction terms
    spline_data <- extract_spline_interaction(fit, data, xvar, by_var, refx,
                                             model_info, term_index, log_scale, ci_level)
    p <- plot_spline_interaction(spline_data, xvar, by_var, xlim, ylim, colors,
                                show_hist, data, bins, model_info, log_scale, refx, show_ref_point, ribbon_ci,
                                xlab, ylab, ylab_right)
  } else {
    # Single smooth/spline term
    spline_data <- extract_spline_data(fit, data, xvar, refx,
                                      model_info, term_index, log_scale, ci_level)
    p <- plot_spline_single(spline_data, xvar, xlim, ylim, show_hist, data, bins,
                           model_info, log_scale, refx, show_ref_point, ribbon_ci,
                           xlab, ylab, ylab_right)
  }

  return(p)
}

#' Detect Model Information
#'
#' @noRd
detect_model_info <- function(fit) {
  model_type <- NULL
  model_family <- NULL
  ylabel <- NULL

  if (inherits(fit, "gam")) {
    model_type <- "gam"
    if (!is.null(fit$family)) {
      family_name <- fit$family$family
      if (grepl("Cox", family_name, ignore.case = TRUE)) {
        model_family <- "cox"
        ylabel <- "Hazard Ratio"
      } else if (family_name == "binomial") {
        model_family <- "binomial"
        ylabel <- "Odds Ratio"
      } else if (family_name == "poisson") {
        model_family <- "poisson"
        ylabel <- "Rate Ratio"
      } else if (family_name == "gaussian") {
        model_family <- "gaussian"
        ylabel <- "Effect"
      } else if (family_name == "quasipoisson") {
        model_family <- "poisson"
        ylabel <- "Rate Ratio"
      } else {
        model_family <- "gaussian"
        ylabel <- "Effect"
      }
    }
  } else if (inherits(fit, "coxph")) {
    model_type <- "coxph"
    model_family <- "cox"
    ylabel <- "Hazard Ratio"
  } else if (inherits(fit, "glm")) {
    model_type <- "glm"
    family_name <- fit$family$family
    if (family_name == "binomial") {
      model_family <- "binomial"
      ylabel <- "Odds Ratio"
    } else if (family_name == "poisson") {
      model_family <- "poisson"
      ylabel <- "Rate Ratio"
    } else if (family_name == "quasipoisson") {
      model_family <- "poisson"
      ylabel <- "Rate Ratio"
    } else if (family_name == "gaussian") {
      model_family <- "gaussian"
      ylabel <- "Effect"
    } else {
      model_family <- "gaussian"
      ylabel <- "Effect"
    }
  } else if (inherits(fit, "lm")) {
    model_type <- "lm"
    model_family <- "gaussian"
    ylabel <- "Effect"
  } else {
    stop("Unsupported model type. Use gam, glm, lm, or coxph models.")
  }

  return(list(
    type = model_type,
    family = model_family,
    ylabel = ylabel
  ))
}

#' Detect Spline Terms and Interactions
#'
#' @noRd
detect_spline_terms <- function(fit, model_info = NULL, xvar = NULL, by_var = NULL) {
  if (is.null(model_info)) {
    model_info <- detect_model_info(fit)
  }

  detected_xvar <- xvar
  detected_by_var <- by_var

  if (model_info$type == "gam") {
    # For GAM, look at smooth terms
    if (length(fit$smooth) > 0) {
      for (i in seq_along(fit$smooth)) {
        smooth <- fit$smooth[[i]]

        # Check if this smooth matches our xvar (if specified)
        if (!is.null(xvar)) {
          if (xvar %in% smooth$term) {
            # Check for by variable in this smooth
            if (!is.null(smooth$by) && smooth$by != "NA" && smooth$by != "NULL") {
              detected_by_var <- smooth$by
            }
            break
          }
        } else {
          # No xvar specified, use first smooth term
          if (is.null(detected_xvar)) {
            detected_xvar <- smooth$term[1]
            # Check for by variable
            if (!is.null(smooth$by) && smooth$by != "NA" && smooth$by != "NULL") {
              detected_by_var <- smooth$by
            }
          }
        }
      }
    }
  } else if (model_info$type == "coxph") {
    # For Cox models, parse formula for spline terms
    formula_str <- deparse(fit$formula)

    if (is.null(detected_xvar)) {
      # Look for pspline(var), ns(var), bs(var)
      patterns <- c("pspline\\(([^,\\)]+)", "ns\\(([^,\\)]+)", "bs\\(([^,\\)]+)")
      for (pattern in patterns) {
        matches <- regmatches(formula_str, gregexpr(pattern, formula_str))
        if (length(matches[[1]]) > 0) {
          var_match <- regmatches(matches[[1]][1], regexec(pattern, matches[[1]][1]))
          if (length(var_match[[1]]) > 1) {
            detected_xvar <- trimws(var_match[[1]][2])
            break
          }
        }
      }
    }

    # Check for interaction with detected xvar
    if (!is.null(detected_xvar) && is.null(detected_by_var)) {
      # Look for patterns like xvar:factor or factor:xvar
      interaction_pattern <- paste0("(", detected_xvar, ":[^\\s\\+\\)]+|[^\\s\\+\\(]+:", detected_xvar, ")")
      interaction_matches <- regmatches(formula_str, gregexpr(interaction_pattern, formula_str))
      if (length(interaction_matches[[1]]) > 0) {
        # Extract the other variable in the interaction
        parts <- strsplit(interaction_matches[[1]][1], ":")[[1]]
        detected_by_var <- parts[parts != detected_xvar][1]
      }
    }
  } else if (model_info$type %in% c("glm", "lm")) {
    # For GLM/LM, parse formula for spline terms
    formula_str <- deparse(fit$formula)

    if (is.null(detected_xvar)) {
      # Look for ns(var), bs(var)
      patterns <- c("ns\\(([^,\\)]+)", "bs\\(([^,\\)]+)")
      for (pattern in patterns) {
        matches <- regmatches(formula_str, gregexpr(pattern, formula_str))
        if (length(matches[[1]]) > 0) {
          var_match <- regmatches(matches[[1]][1], regexec(pattern, matches[[1]][1]))
          if (length(var_match[[1]]) > 1) {
            detected_xvar <- trimws(var_match[[1]][2])
            break
          }
        }
      }
    }

    # Check for interaction with spline term
    if (!is.null(detected_xvar) && is.null(detected_by_var)) {
      # Look for patterns with interaction
      interaction_pattern <- paste0("[nb]s\\(", detected_xvar, "[^\\)]*\\)\\s*[\\*:]\\s*([^\\s\\+]+)")
      interaction_matches <- regmatches(formula_str, regexec(interaction_pattern, formula_str))
      if (length(interaction_matches[[1]]) > 1) {
        detected_by_var <- trimws(interaction_matches[[1]][2])
      }
    }
  }

  return(list(
    xvar = detected_xvar,
    by_var = detected_by_var
  ))
}

#' Extract Spline Data
#'
#' Extract predictions and confidence intervals from fitted models
#'
#' @param fit Fitted model object
#' @param data Data frame
#' @param xvar Variable name
#' @param refx Reference value
#' @param model_info Model information list
#' @param term_index Which smooth term to use (for multiple s() terms)
#' @param log_scale Whether to use log scale
#' @param ci_level Confidence level
#'
#' @return Data frame with predictions
#' @export
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' n <- 100
#' x <- rnorm(n, mean = 50, sd = 10)
#' y <- rbinom(n, 1, plogis(-0.05*(x - 50)))
#' dat <- data.frame(x = x, y = y)
#'
#' # Fit GLM with splines
#' library(splines)
#' fit <- glm(y ~ ns(x, df = 4), family = binomial(), data = dat)
#'
#' # Extract spline data
#' model_info <- list(type = "glm", family = "binomial", ylabel = "Odds Ratio")
#' df <- extract_spline_data(fit, dat, "x", refx = 50, model_info,
#'                           log_scale = FALSE, ci_level = 0.95)
#' head(df)
#'
extract_spline_data <- function(fit, data, xvar, refx,
                               model_info, term_index = 1,
                               log_scale = FALSE, ci_level = 0.95) {

  xv <- data[[xvar]]
  # refx should already be set by splineplot() function

  # Check for multiple smooth terms in GAM
  if (model_info$type == "gam") {
    smooth_terms <- fit$smooth
    xvar_terms <- which(sapply(smooth_terms, function(s) xvar %in% s$term))
    if (length(xvar_terms) > 1) {
      if (term_index > length(xvar_terms)) {
        warning(paste("term_index", term_index, "exceeds number of smooth terms. Using first term."))
        term_index <- 1
      }
      # Use the specified term
      selected_term <- xvar_terms[term_index]
    }
  }

  # Create prediction grid
  gx <- seq(min(xv, na.rm = TRUE), max(xv, na.rm = TRUE), length.out = 200)
  nd <- setNames(data.frame(gx), xvar)

  # Add other variables at their reference values
  for (var in names(data)) {
    if (var != xvar && !(var %in% names(nd))) {
      if (is.numeric(data[[var]])) {
        nd[[var]] <- median(data[[var]], na.rm = TRUE)
      } else if (is.factor(data[[var]])) {
        nd[[var]] <- levels(data[[var]])[1]
      } else {
        nd[[var]] <- data[[var]][1]
      }
    }
  }

  nd_ref <- nd[1, , drop = FALSE]
  nd_ref[[xvar]] <- refx

  z_score <- qnorm((1 + ci_level) / 2)

  if (model_info$type == "gam") {
    # GAM predictions
    Xp <- predict(fit, newdata = nd, type = "lpmatrix")
    Xr <- predict(fit, newdata = nd_ref, type = "lpmatrix")
    b <- coef(fit)
    V <- vcov(fit)

    C <- sweep(Xp, 2, Xr[1, ], FUN = "-")
    diff <- as.numeric(C %*% b)
    se_diff <- sqrt(rowSums((C %*% V) * C))

  } else if (model_info$type %in% c("glm", "lm")) {
    # GLM/LM predictions using model matrix approach
    # Use terms object to get predictor-only model matrix
    tt <- terms(fit)
    Terms <- delete.response(tt)

    # Create model matrices for new data
    mm_nd <- model.matrix(Terms, data = nd)
    mm_ref <- model.matrix(Terms, data = nd_ref)

    b <- coef(fit)
    V <- vcov(fit)

    # Calculate differences from reference
    C <- sweep(mm_nd, 2, mm_ref[1, ], FUN = "-")
    diff <- as.numeric(C %*% b)
    se_diff <- sqrt(rowSums((C %*% V) * C))

  } else if (model_info$type == "coxph") {
    # Check if model contains pspline terms
    has_pspline <- any(grepl("ps\\(", names(fit$coefficients)))

    if (has_pspline) {
      # For pspline, we use a different approach since the basis is not directly accessible
      # Use predict with type="terms" which gives us the smooth function values
      pred <- predict(fit, newdata = nd, type = "terms", se.fit = TRUE)
      pred_ref <- predict(fit, newdata = nd_ref, type = "terms", se.fit = TRUE)

      # Find which term contains our variable
      term_names <- colnames(pred$fit)
      xvar_term <- grep(paste0("pspline\\(", xvar), term_names, value = TRUE)

      if (length(xvar_term) > 0) {
        # Calculate difference from reference
        diff <- pred$fit[, xvar_term] - as.numeric(pred_ref$fit[, xvar_term])

        # For SE, we approximate by assuming independence (not perfect but acceptable)
        # The SE at reference will be naturally small but not exactly 0
        se_at_points <- if (is.matrix(pred$se.fit)) pred$se.fit[, xvar_term] else pred$se.fit
        se_at_ref <- if (is.matrix(pred_ref$se.fit)) pred_ref$se.fit[, xvar_term] else pred_ref$se.fit

        # Approximate SE of difference (assumes independence)
        se_diff <- sqrt(se_at_points^2 + se_at_ref^2)

        # Find the point closest to refx and minimize SE there
        ref_idx <- which.min(abs(gx - refx))
        if (length(ref_idx) > 0 && ref_idx <= length(se_diff)) {
          # At the reference point, the SE should be minimal
          se_diff[ref_idx] <- min(se_diff) * 0.01  # Nearly 0 but not exactly
        }
      } else {
        # Fallback: use linear predictor for the whole model
        pred_lp <- predict(fit, newdata = nd, type = "lp", se.fit = TRUE)
        pred_ref_lp <- predict(fit, newdata = nd_ref, type = "lp", se.fit = TRUE)

        diff <- pred_lp$fit - pred_ref_lp$fit[1]
        se_diff <- sqrt(pred_lp$se.fit^2 + pred_ref_lp$se.fit^2)

        # Minimize SE at reference
        ref_idx <- which.min(abs(gx - refx))
        if (length(ref_idx) > 0) {
          se_diff[ref_idx] <- min(se_diff) * 0.01
        }
      }
    } else {
      # Regular Cox model without pspline (e.g., ns, bs)
      tt <- terms(fit)
      Terms <- delete.response(tt)

      mm_nd <- model.matrix(Terms, data = nd)
      mm_ref <- model.matrix(Terms, data = nd_ref)

      # Remove intercept if present (Cox models don't have intercept)
      if (ncol(mm_nd) > 0 && colnames(mm_nd)[1] == "(Intercept)") {
        mm_nd <- mm_nd[, -1, drop = FALSE]
        mm_ref <- mm_ref[, -1, drop = FALSE]
      }

      b <- coef(fit)
      V <- vcov(fit)

      # Calculate differences from reference
      C <- sweep(mm_nd, 2, mm_ref[1, ], FUN = "-")
      diff <- as.numeric(C %*% b)
      se_diff <- sqrt(rowSums((C %*% V) * C))
    }
  }

  # Calculate y values based on model family and log_scale
  # ALWAYS use log scale internally for ratio metrics to avoid scale issues
  if (model_info$family %in% c("cox", "binomial", "poisson")) {
    # Always keep as log scale for plotting
    y <- diff
    lcl <- diff - z_score * se_diff
    ucl <- diff + z_score * se_diff
    # Label changes based on log_scale
    if (log_scale) {
      ylabel <- paste("Log", model_info$ylabel)
    } else {
      ylabel <- model_info$ylabel
    }
  } else {
    # For gaussian/linear models
    y <- diff
    lcl <- diff - z_score * se_diff
    ucl <- diff + z_score * se_diff
    ylabel <- model_info$ylabel
  }

  df_curve <- data.frame(
    x = gx,
    y = y,
    lcl = lcl,
    ucl = ucl,
    ylabel = ylabel
  )

  return(df_curve)
}

#' Extract Spline Data with Interaction
#'
#' @noRd
extract_spline_interaction <- function(fit, data, xvar, by_var, refx,
                                      model_info, term_index = 1,
                                      log_scale = FALSE, ci_level = 0.95) {

  # refx should already be set by splineplot() function

  by_levels <- if (is.factor(data[[by_var]])) {
    levels(data[[by_var]])
  } else {
    unique(data[[by_var]])
  }

  z_score <- qnorm((1 + ci_level) / 2)

  result_list <- list()

  for (level in by_levels) {
    # Create prediction data for this level
    xv <- data[[xvar]]
    gx <- seq(min(xv, na.rm = TRUE), max(xv, na.rm = TRUE), length.out = 200)
    nd <- data.frame(x = gx)
    names(nd) <- xvar
    # Ensure by_var is factor with correct levels
    if (is.factor(data[[by_var]])) {
      nd[[by_var]] <- factor(rep(level, length(gx)), levels = levels(data[[by_var]]))
    } else {
      nd[[by_var]] <- level
    }

    # Add other variables
    for (var in names(data)) {
      if (!(var %in% c(xvar, by_var)) && !(var %in% names(nd))) {
        if (is.numeric(data[[var]])) {
          nd[[var]] <- median(data[[var]], na.rm = TRUE)
        } else if (is.factor(data[[var]])) {
          nd[[var]] <- factor(levels(data[[var]])[1], levels = levels(data[[var]]))
        } else {
          nd[[var]] <- data[[var]][1]
        }
      }
    }

    nd_ref <- nd[1, , drop = FALSE]
    nd_ref[[xvar]] <- refx

    if (model_info$type == "gam") {
      Xp <- predict(fit, newdata = nd, type = "lpmatrix")
      Xr <- predict(fit, newdata = nd_ref, type = "lpmatrix")
      b <- coef(fit)
      V <- vcov(fit)

      C <- sweep(Xp, 2, Xr[1, ], FUN = "-")
      diff <- as.numeric(C %*% b)
      se_diff <- sqrt(rowSums((C %*% V) * C))

    } else if (model_info$type == "coxph") {
      # Check for pspline terms
      has_pspline <- any(grepl("pspline", names(fit$coefficients)))

      if (has_pspline && grepl("pspline", xvar)) {
        # Handle pspline
        pred <- predict(fit, newdata = nd, type = "terms", se.fit = TRUE)
        pred_ref <- predict(fit, newdata = nd_ref, type = "terms", se.fit = TRUE)

        term_names <- colnames(pred$fit)
        xvar_term <- grep(paste0("pspline\\(", xvar), term_names, value = TRUE)

        if (length(xvar_term) > 0) {
          diff <- pred$fit[, xvar_term] - as.numeric(pred_ref$fit[, xvar_term])
          se_diff <- pred$se.fit[, xvar_term]
          ref_idx <- which.min(abs(gx - refx))
          if (length(ref_idx) > 0) se_diff[ref_idx] <- 0
        } else {
          # Fallback
          tt <- terms(fit)
          Terms <- delete.response(tt)
          mm_nd <- model.matrix(Terms, data = nd)
          mm_ref <- model.matrix(Terms, data = nd_ref)
          b <- coef(fit)
          V <- vcov(fit)
          C <- sweep(mm_nd, 2, mm_ref[1, ], FUN = "-")
          diff <- as.numeric(C %*% b)
          se_diff <- sqrt(rowSums((C %*% V) * C))
        }
      } else {
        # Regular Cox model
        tt <- terms(fit)
        Terms <- delete.response(tt)
        mm_nd <- model.matrix(Terms, data = nd)
        mm_ref <- model.matrix(Terms, data = nd_ref)

        # Remove intercept if present (Cox models don't have intercept)
        if (ncol(mm_nd) > 0 && colnames(mm_nd)[1] == "(Intercept)") {
          mm_nd <- mm_nd[, -1, drop = FALSE]
          mm_ref <- mm_ref[, -1, drop = FALSE]
        }

        b <- coef(fit)
        V <- vcov(fit)
        C <- sweep(mm_nd, 2, mm_ref[1, ], FUN = "-")
        diff <- as.numeric(C %*% b)
        se_diff <- sqrt(rowSums((C %*% V) * C))
      }

    } else {
      # GLM/LM using model matrix approach
      # Use terms object to get predictor-only model matrix
      tt <- terms(fit)
      Terms <- delete.response(tt)

      mm_nd <- model.matrix(Terms, data = nd)
      mm_ref <- model.matrix(Terms, data = nd_ref)

      b <- coef(fit)
      V <- vcov(fit)

      C <- sweep(mm_nd, 2, mm_ref[1, ], FUN = "-")
      diff <- as.numeric(C %*% b)
      se_diff <- sqrt(rowSums((C %*% V) * C))
    }

    # Calculate y values based on model family and log_scale
    # ALWAYS use log scale internally for ratio metrics to avoid scale issues
    if (model_info$family %in% c("cox", "binomial", "poisson")) {
      # Always keep as log scale for plotting
      y <- diff
      lcl <- diff - z_score * se_diff
      ucl <- diff + z_score * se_diff
    } else {
      # For gaussian/linear models
      y <- diff
      lcl <- diff - z_score * se_diff
      ucl <- diff + z_score * se_diff
    }

    result_list[[as.character(level)]] <- data.frame(
      x = gx,
      y = y,
      lcl = lcl,
      ucl = ucl,
      group = as.character(level)
    )
  }

  df_curve <- do.call(rbind, result_list)

  # Add ylabel column based on log_scale
  if (model_info$family %in% c("cox", "binomial", "poisson")) {
    if (log_scale) {
      df_curve$ylabel <- paste("Log", model_info$ylabel)
    } else {
      df_curve$ylabel <- model_info$ylabel
    }
  } else {
    df_curve$ylabel <- model_info$ylabel
  }

  return(df_curve)
}

#' Plot Single Spline
#'
#' @noRd
plot_spline_single <- function(df_curve, xvar, xlim, ylim, show_hist, data, bins,
                              model_info, log_scale, refx, show_ref_point, ribbon_ci = FALSE,
                              xlab = NULL, ylab = NULL, ylab_right = "Percent of Population") {

  xv <- data[[xvar]]

  if (is.null(ylim)) {
    # Auto-determine ylim - data is always in log scale for ratio metrics
    ylim <- range(c(df_curve$lcl, df_curve$ucl), na.rm = TRUE)
    ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
  } else if (model_info$family %in% c("cox", "binomial", "poisson") && !log_scale) {
    # User provided ylim in ratio scale (e.g., c(0.25, 4)), convert to log
    ylim <- log(ylim)
  }

  # Always define axis limits for consistency
  left_min <- ylim[1]
  left_max <- ylim[2]

  # Scale mapping for dual axes (only needed for histogram)
  to_left <- function(pct) left_min + pct/100 * (left_max - left_min)
  from_left <- function(y) 100 * (y - left_min)/(left_max - left_min)

  # Always use floating X-axis for consistency
  gap_size <- (left_max - left_min) * 0.03
  x_axis_y <- left_min - gap_size

  # Histogram data
  if (show_hist) {
    h <- hist(xv, breaks = bins, plot = FALSE, include.lowest = TRUE, right = FALSE)
    pct <- 100 * h$counts / sum(h$counts)

    # Histogram always starts at left_min (0% of secondary axis)
    df_hist <- data.frame(
      xmin = head(h$breaks, -1),
      xmax = tail(h$breaks, -1),
      ymin = left_min,
      ymax = to_left(pct)
    )
  } else {
    df_hist <- NULL
  }

  # Set x-axis range
  if (!is.null(xlim)) {
    x_min <- xlim[1]
    x_max <- xlim[2]
  } else {
    x_min <- min(xv, na.rm = TRUE)
    x_max <- max(xv, na.rm = TRUE)
  }

  axis_lwd <- 0.5

  # Define Y-axis breaks early for manual ticks
  # Use custom ylab if provided, otherwise use default from model
  if (!is.null(ylab)) {
    ylabel <- ylab
  } else {
    ylabel <- df_curve$ylabel[1]
  }

  if (model_info$family %in% c("cox", "binomial", "poisson") && !log_scale) {
    # For ratio scales - create breaks in original scale, then convert to log
    exp_ylim <- exp(ylim)
    if (exp_ylim[2] <= 2) {
      ratio_breaks <- c(0.25, 0.5, 1, 2)
    } else if (exp_ylim[2] <= 4) {
      ratio_breaks <- c(0.25, 0.5, 1, 2, 4)
    } else if (exp_ylim[2] <= 10) {
      ratio_breaks <- c(0.1, 0.25, 0.5, 1, 2, 5, 10)
    } else {
      ratio_breaks <- c(0.1, 0.5, 1, 2, 5, 10, 20)
    }
    # Filter to range and convert to log scale
    ratio_breaks <- ratio_breaks[ratio_breaks >= exp_ylim[1] & ratio_breaks <= exp_ylim[2]]
    y_breaks <- log(ratio_breaks)
    y_labels <- format(ratio_breaks, drop0trailing = TRUE)
  } else {
    # For effect/log scales
    y_breaks <- pretty(ylim, n = 5)
    y_labels <- format(y_breaks, drop0trailing = TRUE)
  }

  # Build plot
  p <- ggplot()

  # Add histogram if requested
  if (show_hist && !is.null(df_hist)) {
    # Filter histogram bars to xlim if specified
    if (!is.null(xlim)) {
      df_hist <- df_hist[df_hist$xmax >= xlim[1] & df_hist$xmin <= xlim[2], ]
      df_hist$xmin <- pmax(df_hist$xmin, xlim[1])
      df_hist$xmax <- pmin(df_hist$xmax, xlim[2])
    }
    p <- p + geom_rect(data = df_hist,
                      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                      fill = "grey80", color = "grey60", linewidth = 0.3, alpha = 0.85)
  }

  # Filter curve data to xlim if specified
  df_curve_plot <- df_curve
  if (!is.null(xlim)) {
    df_curve_plot <- df_curve_plot[df_curve_plot$x >= xlim[1] & df_curve_plot$x <= xlim[2], ]
  }
  # Set values outside ylim to NA so lines are not drawn outside range
  # This breaks the line where it goes out of bounds
  df_curve_plot$y[df_curve_plot$y < ylim[1] | df_curve_plot$y > ylim[2]] <- NA
  df_curve_plot$lcl[df_curve_plot$lcl < ylim[1] | df_curve_plot$lcl > ylim[2]] <- NA
  df_curve_plot$ucl[df_curve_plot$ucl < ylim[1] | df_curve_plot$ucl > ylim[2]] <- NA

  # Add confidence intervals based on ribbon_ci option
  if (ribbon_ci) {
    p <- p +
      # Ribbon style confidence intervals
      geom_ribbon(data = df_curve_plot, aes(x = x, ymin = lcl, ymax = ucl),
                  alpha = 0.2, fill = "grey50") +
      # Main curve
      geom_line(data = df_curve_plot, aes(x = x, y = y), linewidth = 1.2)
  } else {
    p <- p +
      # Dotted line style confidence intervals (default)
      geom_line(data = df_curve_plot, aes(x = x, y = lcl), linetype = "dotted") +
      geom_line(data = df_curve_plot, aes(x = x, y = ucl), linetype = "dotted") +
      # Main curve
      geom_line(data = df_curve_plot, aes(x = x, y = y), linewidth = 1.2)
  }

  # Add reference line - always at 0 since data is in log scale
  ref_y <- 0  # log(1) = 0 for ratio metrics, 0 for linear metrics
  p <- p +
    geom_hline(yintercept = ref_y, linetype = "dashed", linewidth = 0.35)

  # Add reference point (◆)
  if (show_ref_point && !is.null(refx)) {
    if (is.null(xlim) || (refx >= xlim[1] && refx <= xlim[2])) {
      p <- p + annotate("point", x = refx, y = 0, shape = 18, size = 3, colour = "black")
    }
  }

  # Need extra space below x-axis for ticks
  tick_len <- (left_max - left_min) * 0.02

  p <- p +
    # Axis settings
    coord_cartesian(xlim = c(x_min, x_max),
                   ylim = c(x_axis_y - tick_len - 0.02, ylim[2]),
                   clip = "off") +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    labs(x = if (!is.null(xlab)) xlab else xvar)

  # Y-axis settings (y_breaks already defined earlier)

  if (show_hist) {
    if (model_info$family %in% c("cox", "binomial", "poisson") && !log_scale) {
      # Use custom labels for ratio scales
      p <- p + scale_y_continuous(
        limits = c(x_axis_y - tick_len - 0.02, ylim[2]),
        breaks = y_breaks,
        labels = y_labels,
        name = ylabel,
        sec.axis = sec_axis(~ from_left(.), name = ylab_right)
      )
    } else {
      p <- p + scale_y_continuous(
        limits = c(x_axis_y - tick_len - 0.02, ylim[2]),
        breaks = y_breaks,
        labels = y_labels,
        name = ylabel,
        sec.axis = sec_axis(~ from_left(.), name = ylab_right)
      )
    }
  } else {
    if (model_info$family %in% c("cox", "binomial", "poisson") && !log_scale) {
      p <- p + scale_y_continuous(
        limits = c(x_axis_y - tick_len - 0.02, ylim[2]),
        breaks = y_breaks,
        labels = y_labels,
        name = ylabel
      )
    } else {
      p <- p + scale_y_continuous(
        limits = c(x_axis_y - tick_len - 0.02, ylim[2]),
        breaks = y_breaks,
        labels = y_labels,
        name = ylabel
      )
    }
  }

  p <- p +
    # Theme settings
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),  # X축 틱 숨김
      axis.ticks.y = element_line(colour = "black", linewidth = axis_lwd),  # Y축 틱 자동
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 15)),
      plot.margin = margin(t = 10, r = 10, b = 20, l = 10)
    )

  # Always draw floating X-axis for consistency
  p <- p + annotate("segment",
                   x = x_min, xend = x_max,
                   y = x_axis_y, yend = x_axis_y,
                   colour = "black", linewidth = axis_lwd)

  # Y-axis left
  p <- p + annotate("segment",
                   x = x_min, xend = x_min,
                   y = left_min, yend = left_max,
                   colour = "black", linewidth = axis_lwd)

  # Y-axis right
  p <- p + annotate("segment",
                   x = x_max, xend = x_max,
                   y = left_min, yend = left_max,
                   colour = "black", linewidth = axis_lwd)

  # Manual X-axis ticks only (Y축 틱은 theme에서 자동으로 처리)
  x_breaks <- pretty(c(x_min, x_max), n = 5)
  x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  tick_len <- (left_max - left_min) * 0.02

  for(xb in x_breaks) {
    p <- p + annotate("segment",
                     x = xb, xend = xb,
                     y = x_axis_y, yend = x_axis_y - tick_len,  # Draw ticks downward
                     colour = "black", linewidth = axis_lwd)
  }

  return(p)
}

#' Plot Spline with Interaction
#'
#' @noRd
plot_spline_interaction <- function(df_curve, xvar, by_var, xlim, ylim, colors,
                                   show_hist, data, bins, model_info, log_scale,
                                   refx, show_ref_point, ribbon_ci = FALSE,
                                   xlab = NULL, ylab = NULL, ylab_right = "Percent of Population") {

  if (is.null(ylim)) {
    # Auto-determine ylim - data is always in log scale for ratio metrics
    ylim <- range(c(df_curve$lcl, df_curve$ucl), na.rm = TRUE)
    ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
  } else if (model_info$family %in% c("cox", "binomial", "poisson") && !log_scale) {
    # User provided ylim in ratio scale (e.g., c(0.25, 4)), convert to log
    ylim <- log(ylim)
  }

  # Default colors if not provided
  if (is.null(colors)) {
    n_groups <- length(unique(df_curve$group))
    colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")[1:n_groups]
    names(colors) <- unique(df_curve$group)
  }

  # Filter data to xlim if specified
  df_curve_plot <- df_curve
  if (!is.null(xlim)) {
    df_curve_plot <- df_curve_plot[df_curve_plot$x >= xlim[1] & df_curve_plot$x <= xlim[2], ]
  }
  # Set values outside ylim to NA so lines are not drawn outside range
  # This breaks the line where it goes out of bounds
  df_curve_plot$y[df_curve_plot$y < ylim[1] | df_curve_plot$y > ylim[2]] <- NA
  df_curve_plot$lcl[df_curve_plot$lcl < ylim[1] | df_curve_plot$lcl > ylim[2]] <- NA
  df_curve_plot$ucl[df_curve_plot$ucl < ylim[1] | df_curve_plot$ucl > ylim[2]] <- NA

  # Always define axis limits for consistency
  left_min <- ylim[1]
  left_max <- ylim[2]

  # Scale mapping for dual axes (only needed for histogram)
  to_left <- function(pct) left_min + pct/100 * (left_max - left_min)
  from_left <- function(y) 100 * (y - left_min)/(left_max - left_min)

  # X-axis position - always use gap for clean separation
  gap_size <- (left_max - left_min) * 0.03
  if (show_hist) {
    # Use floating X-axis when histogram is shown
    x_axis_y <- left_min - gap_size
  } else {
    # Small gap even without histogram to separate axes
    x_axis_y <- left_min - gap_size
  }

  # Build plot
  p <- ggplot()

  # Add histogram if requested for interaction plots
  if (show_hist) {
    xv <- data[[xvar]]
    h <- hist(xv, breaks = bins, plot = FALSE, include.lowest = TRUE, right = FALSE)
    pct <- 100 * h$counts / sum(h$counts)

    # Histogram always starts at left_min (0% of secondary axis)
    df_hist <- data.frame(
      xmin = head(h$breaks, -1),
      xmax = tail(h$breaks, -1),
      ymin = left_min,
      ymax = to_left(pct)
    )

    # Filter histogram bars to xlim if specified
    if (!is.null(xlim)) {
      df_hist <- df_hist[df_hist$xmax >= xlim[1] & df_hist$xmin <= xlim[2], ]
      df_hist$xmin <- pmax(df_hist$xmin, xlim[1])
      df_hist$xmax <- pmin(df_hist$xmax, xlim[2])
    }

    p <- p + geom_rect(data = df_hist,
                      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                      fill = "grey80", color = "grey60", linewidth = 0.3,
                      alpha = 0.85)
  }

  # Add confidence intervals based on ribbon_ci option
  if (ribbon_ci) {
    p <- p +
      # Ribbon style confidence intervals
      geom_ribbon(data = df_curve_plot, aes(x = x, ymin = lcl, ymax = ucl, fill = group, group = group),
                 alpha = 0.2, color = NA) +
      # Main curves
      geom_line(data = df_curve_plot, aes(x = x, y = y, group = group, color = group), linewidth = 1.2)
  } else {
    p <- p +
      # Dotted line style confidence intervals (default)
      geom_line(data = df_curve_plot, aes(x = x, y = lcl, group = group, color = group), linetype = "dotted") +
      geom_line(data = df_curve_plot, aes(x = x, y = ucl, group = group, color = group), linetype = "dotted") +
      # Main curves
      geom_line(data = df_curve_plot, aes(x = x, y = y, group = group, color = group), linewidth = 1.2)
  }

  # Add reference line - always at 0 since data is in log scale
  ref_y <- 0  # log(1) = 0 for ratio metrics, 0 for linear metrics
  p <- p + geom_hline(yintercept = ref_y, linetype = "dashed", linewidth = 0.35)

  # Add reference point (◆) if requested (and within xlim)
  if (show_ref_point && !is.null(refx)) {
    if (is.null(xlim) || (refx >= xlim[1] && refx <= xlim[2])) {
      p <- p + annotate("point", x = refx, y = ref_y, shape = 18, size = 3, colour = "black")
    }
  }

  # Apply colors
  p <- p +
    scale_color_manual(values = colors, name = by_var) +
    scale_fill_manual(values = colors, name = by_var)

  # Y-axis label
  if (!is.null(ylab)) {
    ylabel <- ylab
  } else {
    # Use default label from df_curve (already set correctly in extract function)
    ylabel <- df_curve$ylabel[1]
  }

  # Define x_min and x_max for coord_cartesian
  if (!is.null(xlim)) {
    x_min <- xlim[1]
    x_max <- xlim[2]
  } else {
    x_min <- min(df_curve$x, na.rm = TRUE)
    x_max <- max(df_curve$x, na.rm = TRUE)
  }

  axis_lwd <- 0.5

  # Define y_breaks for manual ticks
  if (model_info$family %in% c("cox", "binomial", "poisson") && !log_scale) {
    # For ratio scales - create breaks in original scale, then convert to log
    exp_ylim <- exp(ylim)
    if (exp_ylim[2] <= 2) {
      ratio_breaks <- c(0.25, 0.5, 1, 2)
    } else if (exp_ylim[2] <= 4) {
      ratio_breaks <- c(0.25, 0.5, 1, 2, 4)
    } else if (exp_ylim[2] <= 10) {
      ratio_breaks <- c(0.1, 0.25, 0.5, 1, 2, 5, 10)
    } else {
      ratio_breaks <- c(0.1, 0.5, 1, 2, 5, 10, 20)
    }
    # Filter to range and convert to log scale
    ratio_breaks <- ratio_breaks[ratio_breaks >= exp_ylim[1] & ratio_breaks <= exp_ylim[2]]
    y_breaks <- log(ratio_breaks)
    y_labels <- format(ratio_breaks, drop0trailing = TRUE)
  } else {
    # For effect/log scales
    y_breaks <- pretty(ylim, n = 5)
    y_labels <- format(y_breaks, drop0trailing = TRUE)
  }

  # Need extra space below x-axis for ticks
  tick_len <- (left_max - left_min) * 0.02
  # Adjust ylim to include floating X-axis and ticks
  plot_ylim <- c(x_axis_y - tick_len - 0.02, left_max)

  p <- p +
    coord_cartesian(xlim = c(x_min, x_max),
                   ylim = plot_ylim,
                   clip = "off") +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    labs(x = if (!is.null(xlab)) xlab else xvar)

  # Y-axis settings with secondary axis for histogram
  if (show_hist) {
    p <- p + scale_y_continuous(
      limits = plot_ylim,
      breaks = y_breaks,
      labels = y_labels,
      name = ylabel,
      expand = c(0, 0),
      sec.axis = sec_axis(~ from_left(.), name = ylab_right)
    )
  } else {
    p <- p + scale_y_continuous(
      limits = plot_ylim,
      breaks = y_breaks,
      labels = y_labels,
      name = ylabel,
      expand = c(0, 0)
    )
  }

  p <- p +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),  # X축 틱 숨김
      axis.ticks.y = element_line(colour = "black", linewidth = axis_lwd),  # Y축 틱 자동
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 15)),
      legend.position = "right",
      plot.margin = margin(t = 10, r = 10, b = 20, l = 10)
    )

  # Manual axis drawing for interaction plots

  p <- p + annotate("segment",
                   x = x_min, xend = x_max,
                   y = x_axis_y, yend = x_axis_y,
                   colour = "black", linewidth = axis_lwd)

  # Y-axis left
  p <- p + annotate("segment",
                   x = x_min, xend = x_min,
                   y = left_min, yend = left_max,
                   colour = "black", linewidth = axis_lwd)

  # Y-axis right
  p <- p + annotate("segment",
                   x = x_max, xend = x_max,
                   y = left_min, yend = left_max,
                   colour = "black", linewidth = axis_lwd)

  # Manual X-axis ticks only (Y축 틱은 theme에서 자동으로 처리)
  x_breaks <- pretty(c(x_min, x_max), n = 5)
  x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  tick_len <- (left_max - left_min) * 0.02

  for(xb in x_breaks) {
    p <- p + annotate("segment",
                     x = xb, xend = xb,
                     y = x_axis_y, yend = x_axis_y - tick_len,  # Draw ticks downward
                     colour = "black", linewidth = axis_lwd)
  }

  return(p)
}

