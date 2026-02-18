test_that("utility approx matches ROSE for u=(1,1,0,0)", {
  # Replicate Table 1 of Wang et al. 2025: pL=0.4, delta=0.15, alpha=0.8 => n=58
  res_util <- calc_sample_size_utility_approx(
    pL = 0.4, qL = 0.8, delta = 0.15, d = 0.15,
    phi = 0, alpha_L = 0.8, alpha_H = 0.8,
    u = c(1, 1, 0, 0)
  )
  res_rose <- calc_sample_size_rose_approx(pL = 0.4, delta = 0.15,
                                            alpha_L = 0.8, alpha_H = 0.8)
  expect_equal(res_util$n, res_rose$n)
  expect_equal(res_util$lambda_u, res_rose$lambda_u, tolerance = 1e-6)
})

test_that("utility approx achieves target PCS", {
  res <- calc_sample_size_utility_approx(
    pL = 0.3, qL = 0.7, delta = 0.10, d = 0.15,
    phi = 0, alpha_L = 0.8, alpha_H = 0.8
  )
  expect_gte(res$PCS_L, 0.80 - 0.01)
  expect_gte(res$PCS_H, 0.80 - 0.01)
})

test_that("utility exact >= utility approx in n", {
  res_a <- calc_sample_size_utility_approx(
    pL = 0.3, qL = 0.7, delta = 0.10, d = 0.15,
    phi = 0, alpha_L = 0.8, alpha_H = 0.8
  )
  res_e <- calc_sample_size_utility_exact(
    pL = 0.3, qL = 0.7, delta = 0.10, d = 0.15,
    phi = 0, alpha_L = 0.8, alpha_H = 0.8,
    max_n = 200
  )
  expect_gte(res_e$n, res_a$n - 5)  # exact should be close or larger
  expect_gte(res_e$PCS_L, 0.80)
  expect_gte(res_e$PCS_H, 0.80)
})

test_that("direct approx returns valid result", {
  res <- sample_size_direct_approx(
    p = c(0.3, 0.4), q = c(0.5, 0.5),
    delta = 0.10, d = 0.10
  )
  expect_s3_class(res, "ss_direct")
  expect_true(res$n_design > 0)
  expect_equal(nrow(res$scenarios), 2)
})

test_that("ROSE exact >= ROSE approx", {
  ra <- calc_sample_size_rose_approx(pL = 0.3, delta = 0.10)
  re <- calc_sample_size_rose_exact(pL = 0.3, delta = 0.10, max_n = 300)
  expect_gte(re$n, ra$n - 3)
  expect_gte(re$PCS_L, 0.80)
})

test_that("phi_bounds respected in calc_pi", {
  b <- phi_bounds(0.3, 0.5)
  expect_no_error(calc_pi(0.3, 0.5, phi = b["phi_max"] - 0.01))
  expect_error(calc_pi(0.3, 0.5, phi = b["phi_max"] + 0.1))
})
