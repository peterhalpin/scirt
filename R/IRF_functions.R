# Last update: 07/12/2016
  # Functions for 2PL models, mainly to implement MLE and WMLE estimation of theta
  # Main results are in Magis & Raiche, 2012, but here they are implemented only for special case of 2PL. In particular, the following equalities are used
  #   P'  = aPQ
  #   P'' = a^2 PQ(1-2P)
  #   I   = \sum a P'
  #   J   = \sum a P''
  #   I'  = J
  #   J'  = sum a^4 PQ(1-6P)
# User beware: functions not written to check or handle input errors.

#--------------------------------------------------------------------------
#' The IRF for the two parameter logistic model
#'
#' Computes a matrix of probabilities for correct responses under 2PL model.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait
#' @param theta2 not used; included for compatabiity with \code{cIRF}
#' @param sorted not used: included for compatbility with \code{cIRF}
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

IRF <-function(parms, theta1, theta2 = NULL, sorted = F) {
  Z <- matrix(0, nrow = length(theta1), ncol = nrow(parms))
  Z <- Z + theta1
  Z <- t(parms$alpha*(t(Z) - parms$beta))
  1/(1 + exp(-Z))
}

#--------------------------------------------------------------------------
#' First deriviate of 2PL in theta
#'
#' Used for obtaining SEs of theta and computing WMLE.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of first derivatives of IRF
#' @export

dIRF <-function(parms, theta) {
  t(parms$alpha * t(IRF(parms, theta) * (1-IRF(parms, theta))))
}

#--------------------------------------------------------------------------
#' Second derivative of 2PL in theta
#'
#' Used for obtaining SEs of theta and computing WMLE.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of second derivatives of IRF
#' @export

d2IRF <- function(parms, theta) {
  t(parms$alpha^2 * t(dIRF(parms, theta) * (1 - 2 * IRF(parms, theta))))
}

#--------------------------------------------------------------------------
#' Test information function for IRF
#'
#' Computes test information under 2PL model.
#'
#' @param resp a matrix or data.frame containing the binary item responses with \code{NA} for missing responses (used omit missing responses)
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)} vector of test information
#' @export

I <- function(resp, parms, theta) {
  p <- IRF(parms, theta)
  q <- 1 - p
  temp <-  t(parms$alpha^2 * t(p * q))  # faster with 2PL
  # temp <- dIRF(parms, theta)^2 / p / q  # too slow with 2PL
  apply(temp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Numerator of penalty used in WMLE
#'
#' Computes equation 15 of Magis & Raiche (2012). For 2PL, this is also the derivative of the test info.
#'
#' @param resp a matrix or data.frame containing the binary item responses with \code{NA} for missing responses (used omit missing responses)
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)} vector of derivatives of test information
#' @export

J <- function(resp, parms, theta) {
  p <- IRF(parms, theta)
  q <- 1 - p
  temp <- t(parms$alpha^3 * t(p * q * (1 - 2 * p))) # faster with 2PL
  # temp <- dIRF(parms, theta) * d2IRF(parms, theta) / p / q # too slow with 2PL
  apply(temp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Derviative of J
#'
#' Computes the derivative of J (for 2PL, second derivative of test info).
#'
#' @param resp a matrix or data.frame containing the binary item responses with \code{NA} for missing responses (used omit missing responses)
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)} vector of derivatives of J
#' @export

dJ <- function(resp, parms, theta) {
  p <- IRF(parms, theta)
  q <- 1 - p
  temp <- t(parms$alpha^4 * t(p * q * (1 - 6 * p * q) ))
  apply(temp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Standard error of WMLE
#'
#' Computes equation 16 of Magis & Raiche (2012).
#'
#' @param resp a matrix or data.frame containing the binary item responses with \code{NA} for missing responses (used omit missing responses)
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)} vector of standard errors for WLME estimator
#' @export

  SE <- function(resp, parms, theta) {
    i <- I(resp, parms, theta)
    temp <- J(resp, parms, theta)^2 + i * dJ(resp, parms, theta)
    temp <- temp / 2 / i^2 + i
    1/sqrt(temp)
  }

#--------------------------------------------------------------------------
#' Loglikelihood conditional on theta
#'
#' Computes loglikeihood of a response pattern under the 2PL model, given item parms and theta.
#'
#' @param resp a matrix or data.frame containing the binary item responses
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)}-vector of loglikelihoods
#' @export

logL <- function(resp, parms, theta) {
  p <- IRF(parms, theta)
  apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' First derviative of loglikelihood, and Warm's weighted loglikelihood, in theta
#'
#' Used to estimate theta in \code{MLE}.
#'
#' @param resp a matrix or data.frame containing the binary item responses
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @param WMLE logical: use weighted ML or not?
#' @return \code{length(theta)}-vecor of (possibly weigthed) loglikelihoods
#' @export

dlogL <- function(theta, resp, parms, WMLE = T) {
  p <- IRF(parms, theta)
  dlog <- apply(parms$alpha * t(resp - p), 2, sum, na.rm = T)
  if (WMLE) {
    dlog <- dlog + J(resp, parms, theta) / 2 / I(resp, parms, theta)
  }
  dlog
}

#--------------------------------------------------------------------------
#' ML or WML estimation of latent trait for 2PL model
#'
#' Theta is estimated by calling \code{uniroot} on \code{dlogL}. SEs are computed analytically. Value of the loglikelihood at the estimate is also provided. If \code{parallel = T}, the call to uniroot is parallelized via \code{parallel::mclapply}.
#'
#' ToDo: 1. Remove duplicate response patterns when computing MLE; 2. Find a better way of setting the range of x in \code{uniroot}.
#'
#' @param resp a matrix or data.frame containing the binary item responses
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param parallel logical: call \code{parallel:mclapply} instead of looping over \code{nrow(resp)}?
#' @return An \code{nrow(resp)} by 3 matrix of loglikleihoods, theta estimates, analytic SEs for each response pattern in resp.
#' @export

MLE <-function(resp, parms, WMLE = T, parallel = T) {
  out <- data.frame(matrix(0, nrow = nrow(resp), ncol = 3))
  names(out) <- c("logL", "theta", "se")

  fun <- function(resp){
    uniroot(dlogL, c(-10, 10), resp = resp, parms = parms, WMLE = WMLE)$root[1]
  }

  if (parallel) {
    out$theta <- parallel::mclapply(data.frame(t(resp)), fun) %>% unlist
  } else {
    for (i in 1:nrow(resp)) {out$theta[i] <- fun(resp[i,])}
  }
  if (WMLE) {
    out$se <- SE(resp, parms, out$theta)
  } else {
    out$se <- 1/sqrt(I(resp, parms, out$theta))
  }
  out$logL <- logL(resp, parms, out$theta)
  out
}
