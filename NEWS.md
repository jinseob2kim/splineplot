# splineplot 0.2.0

## Major Improvements
* **Y-axis scaling**: Implemented proper ratio scale display for HR/OR/RR when `log_scale=FALSE`
  - Shows actual ratio values (e.g., 0.5, 1, 2) instead of log values
  - Automatically selects appropriate breaks based on data range
  - Removes trailing zeros from labels (1.50 → 1.5)
* **X-axis improvements**:
  - Fixed X-axis tick visibility and direction (now properly pointing downward)
  - Adjusted plot limits to ensure ticks are always visible
* **Histogram alignment**: Base of histogram now correctly aligns with secondary Y-axis 0%
* **Internal refactoring**: Simplified data handling by always using log scale internally for ratio metrics
  - `log_scale` parameter now only affects Y-axis label display
  - More consistent and predictable behavior

## Bug Fixes
* Fixed missing Y-axis labels in GAM interaction plots
* Corrected reference line position (always at y=0 for log scale, which represents ratio=1)
* Fixed ylabel column missing in `extract_spline_interaction()` output

# splineplot 0.1.1

## Bug Fixes
* Fixed Y-axis tick marks display issue in interaction plots
* Fixed tick marks protruding from axes when histogram is shown
* Corrected secondary Y-axis scale for "Percent of Population" in interaction plots
* Fixed X-axis positioning with floating axis for histogram display
* Improved axis tick alignment for both single and interaction plots

## Documentation
* Updated GAM survival examples to use recommended `weights` parameter format
* Removed unnecessary logo reference from README
* Clarified that GAM Cox models should use `time ~ s(predictor), weights = status` format

# splineplot 0.1.0

## Initial Release

### Major Features
* Unified interface for visualizing spline effects from GAM and GLM models
* Support for multiple model types:
  - GAM models from `mgcv` package with `s()`, `te()`, `ti()` smooth terms
  - GLM/LM models with `ns()` and `bs()` splines from `splines` package
  - Cox proportional hazards models from `survival` package
* Automatic detection of:
  - Model type and family
  - Spline terms
  - Interaction variables
* Support for various outcome types:
  - Hazard Ratios (HR) for Cox models
  - Odds Ratios (OR) for logistic models
  - Rate Ratios (RR) for Poisson models
  - Effects for linear/Gaussian models

### Visualization Features
* Publication-ready ggplot2 output
* Customizable confidence intervals:
  - Dotted lines (default)
  - Ribbon/shaded style
* Built-in histogram showing data distribution
* Reference point marking with automatic SE = 0
* Support for interaction terms with by-variable
* Log scale option for ratio outcomes
* Customizable axis labels and limits

### Technical Features
* Automatic handling of different spline basis functions
* Proper reference value centering with SE = 0
* Support for both `Surv()` and weights methods in GAM Cox models
* Limited support for `pspline()` in Cox models

### Known Limitations
* `pspline()` in Cox models has limited support due to internal structure
* Recommend using `ns()` or `bs()` with Cox models for optimal results