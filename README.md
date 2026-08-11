<!-- badges: start -->
[![arXiv](https://img.shields.io/badge/arXiv-2603.15884-b31b1b.svg)](https://arxiv.org/abs/2603.15884)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

# DoseOptDesign

A utility score framework for dose optimization studies with binary
efficacy-safety endpoints. Implements sample size determination and
selection-induced bias characterization as described in:


> Gu, X., Xu, C., Xu, L., & Yuan, Y. (2026). A Utility Score Framework for Dose
> Optimization Studies with Binary Efficacy-Safety Endpoints: Sample Size
> Determination and Bias Characterization. *arXiv preprint arXiv:2603.15884*. https://arxiv.org/abs/2603.15884

## Citation

If you use this package in your research, please cite:
```bibtex
@misc{gu2025utility,
  title={A Utility Score Framework for Dose Optimization Studies with Binary 
         Efficacy-Safety Endpoints: Sample Size Determination and Bias Characterization},
  author={Gu, Xuemin and Xu, Cong and Xu, Lei and Yuan, Ying},
  year={2025},
  eprint={2603.15884},
  archivePrefix={arXiv},
  primaryClass={stat.ME},
  url={https://arxiv.org/abs/2603.15884}
}
```
## Documentation

- **Getting started vignette**: `vignette("getting-started", package = "DoseOptDesign")`
- **Function help**: `?calc_sample_size_utility_approx`
- **Package overview**: `?DoseOptDesign`
- **arXiv manuscript**: https://arxiv.org/abs/2603.15884

## Installation

```r
# Install from local source
install.packages("path/to/DoseOptDesign", repos = NULL, type = "source")

# Or via devtools/remotes from GitHub
# devtools::install_github("xuemingu/DoseOptDesign")
```

## Quick Start

### 1. Sample Size Calculation

```r
library(DoseOptDesign)

# Utility-based design (approximate)
design <- calc_sample_size_utility_approx(
  pL = 0.3, qL = 0.5, delta = 0.10, d = 0.15, phi = 0,
  alpha_L = 0.8, alpha_H = 0.8
)
print(design)

# Utility-based design (exact multinomial)
design_exact <- calc_sample_size_utility_exact(
  pL = 0.3, qL = 0.5, delta = 0.10, d = 0.15, phi = 0,
  alpha_L = 0.8, alpha_H = 0.8
)
print(design_exact)

# ROSE design (efficacy-only, for comparison)
rose <- calc_sample_size_rose_approx(pL = 0.3, delta = 0.10)
print(rose)
```

### 2. Analytical Bias and Type I Error

```r
# Compute selection-induced bias under null
bias <- calc_analytical_bias(
  p = 0.3, q = 0.8, phi = 0,
  u = c(1, 0.8, 0.2, 0),
  n1 = 60, n2 = 140
)
bias$bias_utility_combined   # Cov(X,U)-based plugin estimate
bias$bias_response_combined  # Response-only max bound

# Compute Type I error inflation
t1e <- calc_analytical_type1_error(
  p0 = 0.3,
  bias_combined = bias$bias_utility_combined,
  n_total = 200, alpha = 0.025
)
t1e$type1_z    # Z-test inflated Type I error
t1e$type1_bin  # Binomial test inflated Type I error
```

### 3. Monte Carlo Simulation

```r
# Full simulation with TTE endpoints (requires copula, survival packages)
sim <- compare_bias_methods_fast_enhanced_v2(
  pL = 0.3, pH = 0.3, qL = 0.8, qH = 0.8, phi = 0,
  N1 = 60, N2 = 140,
  u = c(1, 0.8, 0.2, 0),
  perform_tte_analysis = TRUE,
  corr_efficacy_tte = 0.3,
  nSim = 10000,
  Alpha = 0.025, Alpha_tte = 0.025,
  return_raw = FALSE
)
```

### 4. Shiny App

```r
DoseOptDesign::run_app()
```

## Key Functions

| Function | Description |
|---|---|
| `calc_sample_size_utility_approx()` | Normal approximation sample size |
| `calc_sample_size_utility_exact()` | Exact multinomial sample size |
| `calc_sample_size_rose_approx()` | ROSE design (approximate) |
| `calc_sample_size_rose_exact()` | ROSE design (exact) |
| `calc_analytical_bias()` | Selection-induced bias formulas |
| `calc_analytical_type1_error()` | Type I error inflation estimates |
| `compare_bias_methods_fast_enhanced_v2()` | Monte Carlo simulation |
| `calc_pi()` | Joint probabilities from marginals |
| `calc_utility()` | Utility scores from trade-off ratio |

## Dependencies

**Required:** R >= 4.0.0, stats

**Optional (for simulation/TTE):** copula, survival, future, future.apply, parallelly

**Optional (for Shiny app):** shiny, shinythemes, DT, shinyBS
