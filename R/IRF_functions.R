# Last update: 14/06/2017
# User beware: functions not written to check or handle input errors.

# Functions for ML, WML, and MAP estimation of latent traint in 2PL models.
# Results for WML are in Magis & Raiche, 2012, but here they are implemented only for special case of 2PL. In particular, the following equalities are used
  #   P'  = aPQ
  #   P'' = a^2 PQ(1-2P)
  #   Info   = \sum a P'
  #   J   = \sum a P''
  #   Info'  = J
  #   J'  = sum a^4 PQ(1-6P)

#--------------------------------------------------------------------------
#' The IRF for the two parameter logistic model.
#'
#' Computes a matrix of probabilities for correct responses under 2PL model.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities.
#' @export

IRF <-function(parms, theta) {
  Z <- matrix(0, nrow = length(theta), ncol = length(parms[[1]]))
  Z <- Z + theta
  Z <- t(parms$alpha*(t(Z) - parms$beta))
  1/(1 + exp(-Z))
}

#--------------------------------------------------------------------------
#' First derivative of 2PL IRF in theta.
#'
#' Used for obtaining SEs of theta and computing WML.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of first derivatives of IRF.
#' @export

dIRF <-function(parms, theta) {
  t(parms$alpha * t(IRF(parms, theta) * (1-IRF(parms, theta))))
}

#--------------------------------------------------------------------------
#' Second derivative of 2PL IRF in theta.
#'
#' Used for obtaining SEs of theta and computing WML.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of second derivatives of IRF.
#' @export

d2IRF <- function(parms, theta) {
  t(parms$alpha * t(dIRF(parms, theta) * (1 - 2 * IRF(parms, theta))))
}

#--------------------------------------------------------------------------
#' Test information function for 2PL IRF.
#'
#' Computes test information under 2PL model.
#'
#' @param resp a matrix or data.frame containing the binary item responses with \code{NA} for missing responses (used omit missing responses).
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{length(theta)} vector of test information.
#' @export

Info <- function(resp, parms, theta) {
  p <- IRF(parms, theta)*!is.na(resp)
  q <- 1 - p
  temp <-  t(parms$alpha^2 * t(p * q))  # faster with 2PL
  # temp <- dIRF(parms, theta)^2 / p / q  # too slow with 2PL
  apply(temp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Numerator of penalty used in WML.
#'
#' Computes equation 15 of Magis & Raiche (2012). For 2PL, this is also the derivative of the test info.
#'
#' @param resp a matrix or data.frame containing the binary item responses with \code{NA} for missing responses (used omit missing responses).
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{length(theta)} vector of derivatives of test information.
#' @export

J <- function(resp, parms, theta) {
  p <- IRF(parms, theta)*!is.na(resp)
  q <- 1 - p
  temp <- t(parms$alpha^3 * t(p * q * (1 - 2 * p))) # faster with 2PL
  # temp <- dIRF(parms, theta) * d2IRF(parms, theta) / p / q # too slow with 2PL
  apply(temp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Derviative of J.
#'
#' Computes the derivative of J (for 2PL, second derivative of test info).
#'
#' @param resp a matrix or data.frame containing the binary item responses with \code{NA} for missing responses (used omit missing responses).
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{length(theta)} vector of derivatives of J.
#' @export

dJ <- function(resp, parms, theta) {
  p <- IRF(parms, theta)*!is.na(resp)
  q <- 1 - p
  temp <- t(parms$alpha^4 * t(p * q * (1 - 6 * p * q) ))
  apply(temp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Standard error of WML.
#'
#' Computes equation 16 of Magis & Raiche (2012).
#'
#' @param resp a matrix or data.frame containing the binary item responses with \code{NA} for missing responses (used omit missing responses).
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{length(theta)} vector of standard errors for WML estimator.
#' @export

WML_SE <- function(resp, parms, theta) {
    i <- Info(resp, parms, theta)
    temp <- J(resp, parms, theta)^2 + i * dJ(resp, parms, theta)
    temp <- temp / 2 / i^2 + i
    1/sqrt(temp)
  }

#--------------------------------------------------------------------------
#' Log-likelihood of 2PL, conditional on theta.
#'
#' Computes loglikeihood of a response pattern under the 2PL model, given item parms and theta.
#'
#' @param resp a matrix or data.frame containing the binary item responses.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{length(theta)}-vector of log-likelihoods.
#' @export

logL <- function(resp, parms, theta) {
  p <- IRF(parms, theta)
  apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Multiplier for first deriviative of log-likelihood of 2PL.
#'
#' Writing the derivative of the log-likelihood of the 2PL in theta for a single item as M * dIRF, this function computes (a matrix of values of) M, with dIRF being the first deriviative of the 2PL IRF. Called by functions that compute derivative of the log-likelihood.

#' @param resp a matrix or data.frame containing the binary item responses.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{dim(resp)}-matrix of multipliers.
#' @export

M <- function(resp, parms, theta) {
  p <- IRF(parms, theta)
  resp / p - (1 - resp) / (1 - p)
}

#--------------------------------------------------------------------------
#' First derviative of log-likelihood of 2PL, in theta.
#'
#' Used to estimate theta in \code{est_2PL}.
#'
#' @param resp a matrix or data.frame containing the binary item responses.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{length(theta)}- vector of first derivatives of log-likelihood of 2PL, in theta.
#' @export

dlogL <- function(resp, parms, theta) {
  m <- M(resp, parms, theta)
  dp <- dIRF(parms, theta)
  apply(m * dp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Multiplier for second deriviative of log-likelihood of 2PL.
#'
#' Similar to \code{M}, but for the second derivatives. Called by functions that compute second derivative of the log-likelihood.

#' @param resp a matrix or data.frame containing the binary item responses.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @param obs logical: should the observed value be returned? If not, the expected value is returned.
#' @return \code{dim(resp)}-matrix of multipliers.
#' @export

N <- function(resp, parms, theta, obs = F) {
  p <- IRF(parms, theta)
  if (obs) {
    resp / p^2 + (1 - resp) / (1 - p)^2
  } else {
    1 / p / (1 - p) * !is.na(resp)
  }
}

#--------------------------------------------------------------------------
#' Second derviative of log-likelihood of 2PL, in theta.
#'
#' Used to compute standard errors simultaneously estimating individual and group assessments.
#' @param resp a matrix or data.frame containing the binary item responses.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the latent trait.
#' @return \code{length(theta)}-vecor of (possibly weigthed) log-likelihoods.
#' @export

d2logL <- function(resp, parms, theta, obs = T) {
  n <- -1 * N(resp, parms, theta, obs)
  dp <- dIRF(parms, theta)
  d2t <- apply(n * dp * dp, 1, sum, na.rm = T)
  if (obs) {
    m <- M(resp, parms, theta)
    d2p <- d2IRF(parms, theta)
    d2t <- d2t + apply(m * d2p, 1, sum , na.rm = T)
  }
  d2t
}

#--------------------------------------------------------------------------
#' Estimation of latent trait for 2PL model.
#'
#' The latent trait (theta) is estimated by calling \code{uniroot} on q = dlogL + C, where dlogL is the first derivative of the log-likelihood, in theta, and C is a function of theta that depends on the method selected. For \code{method = "ML"}, C = 0. For \code{method = "WML"}, C is given by equation 9 of Warm (1989). For \code{method = "MAP"}, C = - theta. Standard errors (or posterior standard deviations) are computed analytically via the test information of each estimator. The value of logL (not the value of q) at the estimate is also provided. If \code{parallel = T}, the call to \code{uniroot} is parallelized via \code{parallel::mclapply}, but the decrease in runtime is negligible for \code{dim(resp) < c(1000, 100)}
#'
#' ToDo: 1. Remove and index duplicate response patterns. 2. Find a better way of setting the range of x in \code{uniroot}.
#'
#' @param resp a matrix or data.frame containing the binary item responses.
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param method one of \code{"ML", "WML", "MAP"}.
#' @param logical: call \code{parallel:mclapply} instead of looping over \code{nrow(resp)}?
#' @return An \code{nrow(resp)} by 3 matrix of loglikleihoods, theta estimates, analytic SEs for each response pattern in resp.
#' @export

est_2PL <-function(resp, parms, method = "ML", parallel = T) {
  n_items <- length(parms[[1]])
  n_obs <- nrow(resp)
  out <- data.frame(matrix(0, nrow = n_obs, ncol = 3))
  names(out) <- c("logL", "theta", "se")

  if (method == "ML") {
    obj <- function(theta, resp, parms) { dlogL(resp, parms, theta) }
    se <- function(resp, parms, theta) { 1 / sqrt(Info(resp, parms, theta)) }
  }

  if (method == "WML") {
    obj <- function(theta, resp, parms) {
      dlogL(resp, parms, theta) + J(resp, parms, theta) / 2 / Info(resp, parms, theta)
    }
    se <- function(resp, parms, theta) { WML_SE(resp, parms, theta) }
  }

  if (method == "MAP") {
    obj <- function(theta, resp, parms){ dlogL(resp, parms, theta) - theta }
    se <- function(resp, parms, theta) { 1 / sqrt(Info(resp, parms, theta) + 1) }
  }

  fun <- function(resp) {
    if (method == "ML" & sum(resp) %in% c(n_items, 0)) {
      NA
    } else {
      uniroot(obj, c(-10, 10), resp = resp, parms = parms)$root[1]
    }
  }

  if (parallel) {
    out$theta <- parallel::mclapply(data.frame(t(resp)), fun) %>% unlist
  } else {
    for (i in 1:nrow(resp)) {out$theta[i] <- fun(resp[i,])}
  }
  out$se <- se(resp, parms, out$theta)
  out$logL <- logL(resp, parms, out$theta)
  out
}

#--------------------------------------------------------------------------
#' Simulate item responses from the 2PL.

#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta the vector of latent traits.
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of binary responses.
#' @export

sim_2PL <- function(parms, theta) {
  n_row <- length(theta)
  n_col <- nrow(parms)
  r <- array(runif(n_row * n_col), dim = c(n_row, n_col))
  p <- IRF(parms, theta)
  out <- ifelse(p > r, 1, 0)
  colnames(out) <- row.names(parms)
  out
}
