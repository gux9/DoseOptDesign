# dosopt: Dose Optimization Design Using Utility Score Framework

<!-- badges: start -->
[![R-CMD-check](https://github.com/xuemin-gu/dosopt/workflows/R-CMD-check/badge.svg)](https://github.com/xuemin-gu/dosopt/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/dosopt)](https://CRAN.R-project.org/package=dosopt)
<!-- badges: end -->

## Overview

**dosopt** implements a frequentist framework for one-stage randomized dose optimization studies with binary efficacy and safety endpoints under FDA's Project Optimus initiative. The package provides:

- **Utility-based sample size calculation** (approximate normal and exact multinomial methods)
- **ROSE design** (efficacy-only special case; Wang et al. 2025)
- **Multi-scenario direct approach** for robust sample size determination
- **Selection-induced bias** characterization for binary and time-to-event endpoints
- **Type I error inflation** formulas for downstream confirmatory trials
- **Monte Carlo simulation** for design validation
- **Interactive Shiny application** for design exploration

## Installation

```r
# Install from CRAN
install.packages("dosopt")

# Or install the development version from GitHub
# install.packages("devtools")
devtools::install_github("xuemin-gu/dosopt")
```

## Quick Start

### Utility-Based Sample Size (Approximate)

```r
library(dosopt)

# Compute sample size per arm using margin-based utility
# Scenario: pL = 0.3, qL = 0.7, delta = 0.10, d = 0.15
# Both arms need PCS >= 80%

res <- calc_sample_size_utility_approx(
  pL    = 0.3,   # efficacy rate at lower dose
  qL    = 0.7,   # toxicity rate at lower dose
  delta = 0.10,  # minimum clinically meaningful efficacy margin
  d     = 0.15,  # maximum acceptable toxicity margin
  phi   = 0,     # efficacy-safety correlation (0 = independence)
  alpha_L = 0.8, # target PCS for lower dose
  alpha_H = 0.8  # target PCS for higher dose
)
print(res)
```

### Exact Multinomial Sample Size

```r
res_exact <- calc_sample_size_utility_exact(
  pL    = 0.3, qL = 0.7,
  delta = 0.10, d = 0.15,
  phi   = 0,
  alpha_L = 0.8, alpha_H = 0.8,
  max_n = 300
)
print(res_exact)
```

### ROSE Design (Efficacy-Only)

```r
# Replicates Wang et al. (2025) Table 1
res_rose <- calc_sample_size_rose_approx(
  pL    = 0.4,
  delta = 0.15,
  alpha_L = 0.8, alpha_H = 0.8
)
print(res_rose)
```

### Bias Characterization

```r
# Analytical bias after two-stage selection
bias_res <- calc_bias(
  pL    = 0.3, qL = 0.7,
  delta = 0.10, d = 0.15,
  phi   = 0, n1 = 60, n2 = 120
)
print(bias_res)
```

### Type I Error Inflation (Z-test)

```r
t1e_res <- calc_type1_error(
  bias_stage1 = bias_res$bias_stage1,
  n1 = 60, n2 = 120,
  p0 = 0.3, alpha = 0.05,
  test = "z"
)
print(t1e_res)
```

### Interactive Shiny Application

```r
run_dosopt_app()
```

## Methodology

The utility score for a patient with outcome $(X, Y)$ (efficacy, no toxicity) is:

$$U = u_1 X(1-Y) + u_2 X Y + u_3 (1-X)(1-Y) + u_4 (1-X)Y$$

Under the independence constraint $u_1 - u_2 - u_3 + u_4 = 0$, with canonical weights $u_2 = 1/(1+r)$, $u_3 = r/(1+r)$ where $r = \delta/d$.

The jointly optimal sample size satisfies:

$$n = \left[\frac{z_{\alpha_L}\sqrt{v_L} - z_{1-\alpha_H}\sqrt{v_H}}{\Delta\mu_H - \Delta\mu_L}\right]^2$$

Selection-induced bias in the confirmatory stage:

$$\text{Bias} \approx \frac{\text{Cov}(X, U)}{\sigma_U \sqrt{n_1}} \cdot \frac{1}{\sqrt{\pi}} \exp\!\left(-\frac{\lambda_u^2 n_1}{4\sigma_U^2}\right)$$

See the package vignettes for full mathematical details and worked examples.

## Vignettes

| Vignette | Description |
|----------|-------------|
| `vignette("introduction", package = "dosopt")` | Quick start and basic examples |
| `vignette("sample-size-methods", package = "dosopt")` | Comparing approximate, exact, and direct approaches |
| `vignette("bias-and-type1", package = "dosopt")` | Bias propagation and Type I error inflation |
| `vignette("tte-endpoints", package = "dosopt")` | Time-to-event endpoints and surrogacy framework |

## Citation

If you use **dosopt** in your research, please cite:

> Gu X, Xu C, Xu L, Yuan Y (2025). "A Utility Score-Based Dose Optimization Design Under FDA's Project Optimus Initiative." *Manuscript in preparation.*

```bibtex
@Article{dosopt2025,
  author  = {Xuemin Gu and Cong Xu and Lei Xu and Ying Yuan},
  title   = {A Utility Score-Based Dose Optimization Design Under {FDA}'s {Project Optimus} Initiative},
  journal = {Manuscript in preparation},
  year    = {2025},
}
```

## Related Work

- Wang et al. (2025). ROSE: Randomized dose Optimization design with a Selection rule based on Efficacy.
- FDA Project Optimus: <https://www.fda.gov/patients/drug-development-process/step-3-clinical-research#Optimus>
- BOIN design: Yuan Y, Hess KR, Hilsenbeck SG, Gilbert MR (2016). *J Clin Oncol*.

## License

GPL-3 © Xuemin Gu
