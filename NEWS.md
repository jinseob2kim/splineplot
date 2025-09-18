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