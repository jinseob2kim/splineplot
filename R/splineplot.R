#' Spline Plot for GAM and GLM Models
#'
#' Create ggplot2 visualizations of smooth or spline effects from GAM and GLM models.
#' Supports Linear, Logistic, Poisson, and Cox models with interaction terms.
#' Handles GAM smooth terms (s(), te(), ti()), GLM splines (ns(), bs()), and Cox pspline().
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
#' \dontrun{
#' library(mgcv)
#' library(survival)
#'
#' # GAM example with Cox
#' fit <- gam(Surv(time, status) ~ s(age), family = cox.ph(), data = mydata)
#' splineplot(fit, mydata, xvar = "age")
#'
#' # Logistic regression with spline
#' fit2 <- glm(y ~ ns(age, 4), family = binomial, data = mydata)
#' splineplot(fit2, mydata, xvar = "age", log_scale = TRUE)
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
  if (model_info$family %in% c("cox", "binomial", "poisson")) {
    if (log_scale) {
      y <- diff
      lcl <- diff - z_score * se_diff
      ucl <- diff + z_score * se_diff
      ylabel <- paste("Log", model_info$ylabel)
    } else {
      y <- exp(diff)
      lcl <- exp(diff - z_score * se_diff)
      ucl <- exp(diff + z_score * se_diff)
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
    nd[[by_var]] <- level

    # Add other variables
    for (var in names(data)) {
      if (!(var %in% c(xvar, by_var)) && !(var %in% names(nd))) {
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
    if (model_info$family %in% c("cox", "binomial", "poisson")) {
      if (log_scale) {
        y <- diff
        lcl <- diff - z_score * se_diff
        ucl <- diff + z_score * se_diff
      } else {
        y <- exp(diff)
        lcl <- exp(diff - z_score * se_diff)
        ucl <- exp(diff + z_score * se_diff)
      }
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
    # Auto-determine ylim based on model type and scale
    if (model_info$family %in% c("cox", "binomial", "poisson") && !log_scale) {
      # For ratio scales (HR, OR, RR)
      data_range <- range(c(df_curve$lcl, df_curve$ucl), na.rm = TRUE)
      # Round to nice values
      ylim <- c(
        max(0.1, floor(data_range[1] * 10) / 10),
        ceiling(data_range[2] * 2) / 2
      )
      # Common ranges
      if (ylim[1] >= 0.2 && ylim[2] <= 2.5) {
        ylim <- c(0.25, 2.0)
      } else if (ylim[1] >= 0.1 && ylim[2] <= 5) {
        ylim <- c(0.1, 5.0)
      }
    } else {
      # For effect/log scales
      ylim <- range(c(df_curve$lcl, df_curve$ucl), na.rm = TRUE)
      ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
    }
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

    df_hist <- data.frame(
      xmin = head(h$breaks, -1),
      xmax = tail(h$breaks, -1),
      y = to_left(pct)
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
    # For ratio scales - create nice breaks
    if (ylim[2] <= 2) {
      y_breaks <- seq(ylim[1], ylim[2], by = 0.25)
    } else if (ylim[2] <= 5) {
      y_breaks <- seq(ylim[1], ylim[2], by = 0.5)
    } else if (ylim[2] <= 10) {
      y_breaks <- seq(ceiling(ylim[1]), floor(ylim[2]), by = 1)
    } else {
      y_breaks <- pretty(ylim, n = 5)
    }
    # Remove breaks outside ylim
    y_breaks <- y_breaks[y_breaks >= ylim[1] & y_breaks <= ylim[2]]
  } else {
    # For effect/log scales
    y_breaks <- pretty(ylim, n = 5)
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
                      aes(xmin = xmin, xmax = xmax, ymin = x_axis_y, ymax = y),
                      fill = "grey80", color = "grey60", linewidth = 0.3, alpha = 0.85)
  }

  # Filter curve data to xlim if specified
  df_curve_plot <- df_curve
  if (!is.null(xlim)) {
    df_curve_plot <- df_curve[df_curve$x >= xlim[1] & df_curve$x <= xlim[2], ]
  }

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

  # Add reference line based on model family and scale
  if (model_info$family %in% c("cox", "binomial", "poisson") && !log_scale) {
    ref_y <- 1
    p <- p +
      geom_hline(yintercept = ref_y, linetype = "dashed", linewidth = 0.35)

    # Add reference point (◆) for HR/OR/RR
    if (show_ref_point && !is.null(refx)) {
      if (is.null(xlim) || (refx >= xlim[1] && refx <= xlim[2])) {
        p <- p + annotate("point", x = refx, y = 1, shape = 18, size = 3, colour = "black")
      }
    }
  } else {
    ref_y <- 0
    p <- p +
      geom_hline(yintercept = ref_y, linetype = "dashed", linewidth = 0.35)

    # Add reference point for effect scale
    if (show_ref_point && !is.null(refx)) {
      if (is.null(xlim) || (refx >= xlim[1] && refx <= xlim[2])) {
        p <- p + annotate("point", x = refx, y = 0, shape = 18, size = 3, colour = "black")
      }
    }
  }

  p <- p +
    # Axis settings
    coord_cartesian(xlim = c(x_min, x_max),
                   ylim = c(x_axis_y - 0.05, ylim[2]),
                   clip = "off") +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    labs(x = if (!is.null(xlab)) xlab else xvar)

  # Y-axis settings (y_breaks already defined earlier)

  if (show_hist) {
    p <- p + scale_y_continuous(
      limits = c(x_axis_y - 0.05, ylim[2]),
      breaks = y_breaks,
      labels = function(x) ifelse(x >= ylim[1], format(x, nsmall = 2), ""),
      name = ylabel,
      sec.axis = sec_axis(~ from_left(.), name = ylab_right)
    )
  } else {
    p <- p + scale_y_continuous(
      limits = c(x_axis_y - 0.05, ylim[2]),  # Same as coord_cartesian
      breaks = y_breaks,
      labels = function(x) ifelse(x >= ylim[1], format(x, nsmall = 2), ""),
      name = ylabel
    )
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
                     y = x_axis_y, yend = x_axis_y - tick_len,
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
    # Auto-determine ylim based on model type and scale
    if (model_info$family %in% c("cox", "binomial", "poisson") && !log_scale) {
      # For ratio scales (HR, OR, RR)
      data_range <- range(c(df_curve$lcl, df_curve$ucl), na.rm = TRUE)
      # Round to nice values
      ylim <- c(
        max(0.1, floor(data_range[1] * 10) / 10),
        ceiling(data_range[2] * 2) / 2
      )
      # Common ranges
      if (ylim[1] >= 0.2 && ylim[2] <= 2.5) {
        ylim <- c(0.25, 2.0)
      } else if (ylim[1] >= 0.1 && ylim[2] <= 5) {
        ylim <- c(0.1, 5.0)
      }
    } else {
      # For effect/log scales
      ylim <- range(c(df_curve$lcl, df_curve$ucl), na.rm = TRUE)
      ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
    }
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
    df_curve_plot <- df_curve[df_curve$x >= xlim[1] & df_curve$x <= xlim[2], ]
  }

  # Build plot
  p <- ggplot(df_curve_plot, aes(x = x, group = group, color = group))

  # Add histogram if requested for interaction plots
  if (show_hist) {
    xv <- data[[xvar]]
    h <- hist(xv, breaks = bins, plot = FALSE, include.lowest = TRUE, right = FALSE)

    # Scale histogram to fit in bottom 20% of plot
    hist_max <- max(h$counts)
    hist_scale <- diff(ylim) * 0.2 / hist_max
    hist_base <- ylim[1]

    df_hist <- data.frame(
      xmin = head(h$breaks, -1),
      xmax = tail(h$breaks, -1),
      ymax = hist_base + h$counts * hist_scale,
      ymin = rep(hist_base, length(h$counts))
    )

    p <- p + geom_rect(data = df_hist,
                      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                      fill = "grey80", color = "grey60", linewidth = 0.3,
                      alpha = 0.85, inherit.aes = FALSE)
  }

  # Add confidence intervals based on ribbon_ci option
  if (ribbon_ci) {
    p <- p +
      # Ribbon style confidence intervals
      geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = group),
                 alpha = 0.2, color = NA) +
      # Main curves
      geom_line(aes(y = y), linewidth = 1.2)
  } else {
    p <- p +
      # Dotted line style confidence intervals (default)
      geom_line(aes(x = x, y = lcl, group = group, color = group), linetype = "dotted") +
      geom_line(aes(x = x, y = ucl, group = group, color = group), linetype = "dotted") +
      # Main curves
      geom_line(aes(y = y), linewidth = 1.2)
  }

  # Add reference line and points
  if (model_info$family %in% c("cox", "binomial", "poisson") && !log_scale) {
    ref_y <- 1
    p <- p + geom_hline(yintercept = ref_y, linetype = "dashed", linewidth = 0.35)
  } else {
    ref_y <- 0
    p <- p + geom_hline(yintercept = ref_y, linetype = "dashed", linewidth = 0.35)
  }

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
    ylabel <- if (log_scale && model_info$family %in% c("cox", "binomial", "poisson")) {
      paste("Log", model_info$ylabel)
    } else {
      model_info$ylabel
    }
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

  p <- p +
    coord_cartesian(xlim = c(x_min, x_max), ylim = ylim, clip = "off") +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, name = ylabel) +
    labs(x = if (!is.null(xlab)) xlab else xvar) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = axis_lwd),  # Show auto ticks
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )

  # Manual axis drawing for interaction plots

  # X-axis at bottom
  p <- p + annotate("segment",
                   x = x_min, xend = x_max,
                   y = ylim[1], yend = ylim[1],
                   colour = "black", linewidth = axis_lwd)

  # Y-axis left
  p <- p + annotate("segment",
                   x = x_min, xend = x_min,
                   y = ylim[1], yend = ylim[2],
                   colour = "black", linewidth = axis_lwd)

  # Y-axis right
  p <- p + annotate("segment",
                   x = x_max, xend = x_max,
                   y = ylim[1], yend = ylim[2],
                   colour = "black", linewidth = axis_lwd)

  # No manual ticks for interaction plots - let ggplot handle them

  return(p)
}

