test_that("calc_pi returns valid probabilities", {
  pi <- calc_pi(0.4, 0.8, phi = 0)
  expect_length(pi, 4)
  expect_true(all(pi >= 0))
  expect_equal(sum(pi), 1, tolerance = 1e-9)
})

test_that("calc_pi respects Frechet bounds", {
  expect_error(calc_pi(0.4, 0.8, phi = 0.99))
  expect_error(calc_pi(0.4, 0.8, phi = -0.99))
})

test_that("phi_bounds returns valid range", {
  bounds <- phi_bounds(0.4, 0.8)
  expect_true(bounds["phi_min"] < bounds["phi_max"])
  expect_true(bounds["phi_min"] <= 0)
})

test_that("calc_utility is correct for r=1", {
  u <- calc_utility(1)
  expect_equal(u[1], 1); expect_equal(u[4], 0)
  expect_equal(u[2], u[3])
  expect_equal(u[2], 0.5)
})

test_that("calc_utility_moments are non-negative", {
  pi_v <- calc_pi(0.4, 0.8, 0)
  u_v  <- calc_utility(1)
  mom  <- calc_utility_moments(pi_v, u_v)
  expect_true(mom$mu >= 0 && mom$mu <= 1)
  expect_true(mom$sigma2 > 0)
})
