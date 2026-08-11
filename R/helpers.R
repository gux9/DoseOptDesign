# Helper functions for DoseOptDesign

#' Calculate Joint Probabilities from Marginals and Correlation
#' 
#' Computes the joint probabilities for four outcome categories:
#' (1) Response & No AE, (2) Response & AE, (3) No Response & No AE, (4) No Response & AE
#' 
#' @param p Response probability
#' @param q No adverse event probability  
#' @param phi Correlation coefficient between response and no-AE (-1 to 1)
#' @return Vector of joint probabilities (pi1, pi2, pi3, pi4)
#' @export
calc_pi <- function(p, q, phi) {
  if (phi < -1 || phi > 1) stop("phi must be in [-1, 1]")
  
  sqrt_term <- sqrt(p * (1 - p) * q * (1 - q))
  pi1 <- p * q + phi * sqrt_term       # Response & No AE
  pi2 <- p - pi1                        # Response & AE
  pi3 <- q - pi1                        # No Response & No AE  
  pi4 <- 1 - p - q + pi1                # No Response & AE
  
  pi <- c(pi1, pi2, pi3, pi4)
  
  if (any(pi < 0) || any(pi > 1)) {
    stop("Invalid joint probabilities. phi may be incompatible with p and q.")
  }
  
  return(pi)
}

#' Calculate Utility Scores from Trade-off Ratio
#' 
#' Computes utility scores u = (u1, u2, u3, u4) based on the trade-off ratio
#' r = delta/d, which represents the relative importance of efficacy vs safety.
#' 
#' @param r Trade-off ratio (delta/d), must be positive
#' @return Vector of utility scores (u1, u2, u3, u4)
#'         u1 = 1 (response & no AE)
#'         u2 = 1/(1+r) (response & AE)
#'         u3 = r/(1+r) (no response & no AE)
#'         u4 = 0 (no response & AE)
#' @export
calc_utility <- function(r) {
  if (r <= 0) stop("r must be positive")
  c(1, 1/(1+r), r/(1+r), 0)
}

#' Calculate Mean and Variance of Utility Score
#' 
#' Computes the first two moments of the utility distribution
#' 
#' @param pi Joint probabilities (length 4)
#' @param u Utility scores (length 4)
#' @return List with:
#'   - mu: Expected utility E[U]
#'   - sigma2: Variance Var(U)
#' @export
calc_utility_moments <- function(pi, u) {
  mu <- sum(u * pi)
  sigma2 <- sum(u^2 * pi) - mu^2
  list(mu = mu, sigma2 = max(sigma2, 1e-10))
}

#' Compute Binomial PMF Vector
#'
#' @param p Success probability
#' @param n Sample size
#' @return Vector of probabilities P(X=0), P(X=1), ..., P(X=n)
#' @keywords internal
#' Compute PMF of Sum for Binomial (special case of multinomial)
#'
#' @param p Success probability
#' @param n Sample size
#' @return PMF vector of length (n + 1) where pmf[k+1] = P(X = k)
compute_pmf_binomial <- function(p, n) {
  dbinom(0:n, size = n, prob = p)
}
