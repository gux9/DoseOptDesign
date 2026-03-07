test_that("calc_pi returns valid probabilities", {
  pi <- calc_pi(0.3, 0.8, 0)
  expect_equal(length(pi), 4)
  expect_true(all(pi >= 0))
  expect_equal(sum(pi), 1, tolerance = 1e-10)
  
  # Check marginals
  expect_equal(pi[1] + pi[2], 0.3, tolerance = 1e-10)  # p

  expect_equal(pi[1] + pi[3], 0.8, tolerance = 1e-10)  # q
})

test_that("calc_pi handles correlation", {
  pi_pos <- calc_pi(0.4, 0.6, 0.3)
  pi_neg <- calc_pi(0.4, 0.6, -0.3)
  pi_zero <- calc_pi(0.4, 0.6, 0)
  
  # Positive correlation increases pi1 (joint success)
  expect_true(pi_pos[1] > pi_zero[1])
  expect_true(pi_neg[1] < pi_zero[1])
})

test_that("calc_utility produces correct scores", {
  u <- calc_utility(1)  # r = 1: equal weight
  expect_equal(u, c(1, 0.5, 0.5, 0))
  
  u2 <- calc_utility(2/3)  # r = delta/d = 0.10/0.15
  expect_equal(u2, c(1, 0.6, 0.4, 0), tolerance = 1e-10)
})

test_that("calc_utility_moments computes correctly", {
  pi <- c(0.24, 0.06, 0.56, 0.14)  # p=0.3, q=0.8, phi=0
  u <- c(1, 0.6, 0.4, 0)
  m <- calc_utility_moments(pi, u)
  
  expected_mu <- sum(u * pi)
  expected_sigma2 <- sum(u^2 * pi) - expected_mu^2
  
  expect_equal(m$mu, expected_mu, tolerance = 1e-10)
  expect_equal(m$sigma2, expected_sigma2, tolerance = 1e-10)
})

test_that("ROSE approximate matches published Table 1", {
  # From Wang et al. (2025), Table 1: pL=0.4, delta=0.15, alpha=0.8 -> n=58
  rose <- calc_sample_size_rose_approx(pL = 0.4, delta = 0.15,
                                        alpha_L = 0.8, alpha_H = 0.8)
  expect_equal(rose$n, 58)
})

test_that("Utility approx is smaller than ROSE for same parameters", {
  util <- calc_sample_size_utility_approx(
    pL = 0.3, qL = 0.5, delta = 0.10, d = 0.15, phi = 0,
    alpha_L = 0.8, alpha_H = 0.8
  )
  rose <- calc_sample_size_rose_approx(pL = 0.3, delta = 0.10,
                                        alpha_L = 0.8, alpha_H = 0.8)
  expect_true(util$n < rose$n)
})

test_that("calc_analytical_bias returns correct structure", {
  bias <- calc_analytical_bias(
    p = 0.3, q = 0.8, phi = 0,
    u = c(1, 0.8, 0.2, 0),
    n1 = 60, n2 = 140
  )
  
  expect_true(is.list(bias))
  expect_true(all(c("bias_utility_combined", "bias_response_combined",
                     "bias_U_score_combined", "cov_XU", "sigma_U") %in% names(bias)))
  
  # Utility bias should be less than response max bound
  expect_true(bias$bias_utility_combined <= bias$bias_response_combined)
  
  # Both biases should be positive
  expect_true(bias$bias_utility_combined > 0)
  expect_true(bias$bias_response_combined > 0)
})

test_that("calc_analytical_type1_error inflates above nominal", {
  bias <- calc_analytical_bias(
    p = 0.3, q = 0.8, phi = 0,
    u = c(1, 0.8, 0.2, 0),
    n1 = 60, n2 = 140
  )
  
  t1e <- calc_analytical_type1_error(
    p0 = 0.3,
    bias_combined = bias$bias_utility_combined,
    n_total = 200, alpha = 0.025
  )
  
  # Type I error should exceed nominal alpha
  expect_true(t1e$type1_z > 0.025)
  expect_true(t1e$type1_bin > 0.025)
  
  # Z-test should show more inflation than binomial
  expect_true(t1e$type1_z > t1e$type1_bin)
})

test_that("Bias max bound equals response-only selection case", {
  # When U = X (ROSE case), Cov(X,U) = Var(X), so utility bias = max bound
  bias <- calc_analytical_bias(
    p = 0.4, q = 0.8, phi = 0,
    u = c(1, 1, 0, 0),  # ROSE: U = X
    n1 = 60, n2 = 140
  )
  
  expect_equal(bias$bias_utility_combined, bias$bias_response_combined,
               tolerance = 1e-10)
})

test_that("PCS targets are met in approximate design", {
  design <- calc_sample_size_utility_approx(
    pL = 0.3, qL = 0.5, delta = 0.15, d = 0.15, phi = 0,
    alpha_L = 0.8, alpha_H = 0.8
  )
  
  expect_true(design$PCS_L >= 0.8 - 0.01)  # Allow small numerical tolerance
  expect_true(design$PCS_H >= 0.8 - 0.01)
})
