# dosopt 0.1.0

## Initial Release

### New Features

* `calc_sample_size_utility_approx()`: Normal approximation sample size for utility-based dose optimization (Modes 1–3).
* `calc_sample_size_utility_exact()`: Exact multinomial PMF-based sample size via dynamic programming.
* `calc_sample_size_rose_approx()` and `calc_sample_size_rose_exact()`: Sample size for the ROSE design (efficacy-only; Wang et al. 2025).
* `sample_size_direct_approx()` and `sample_size_direct_exact()`: Multi-scenario direct approach (Mode 4) with optional `lambda_u` optimization.
* `calc_bias()`: Analytical selection-induced bias for binary endpoints after two-stage dose optimization.
* `calc_bias_TTE()`: Bias bounds for time-to-event endpoints using surrogacy correlation.
* `calc_type1_error()`: Type I error inflation for Z-test and exact binomial test.
* `calc_type1_error_landmark()`: Type I error for landmark survival analysis.
* `calc_type1_error_exp()`: Type I error for one-arm exponential log-scale test.
* `calc_type1_error_cox()`: Type I error for two-arm Cox regression via delta-method chain.
* `calc_expected_events()`: Expected events under exponential survival with uniform accrual.
* `simulate_dose_opt()`: Full two-stage Monte Carlo simulation with optional parallel processing.
* `calc_pi()`, `phi_bounds()`, `calc_utility()`, `calc_utility_moments()`: Core computational utilities.
* `run_dosopt_app()`: Launch interactive Shiny application.
* S3 `print` methods for `dose_design`, `rose_design`, and `ss_direct` objects.
