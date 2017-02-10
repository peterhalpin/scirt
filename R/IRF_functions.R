# Last update: 07/12/2016
  # Functions for 2PL models, mainly to implement MLE and WMLE estimation of theta
  # Main results are in Magis & Raiche, 2012, but only for special case of 2PL.
  # In particular, the following equalities are used
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
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait
#' @param theta2 not used; included for formal compatabiity with collaborative IRFs
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

IRF <-function(parms, theta1, theta2 = NULL) {
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
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of first derivatives of IRF
#' @export

dIRF <-function(parms, theta) {
  t(parms$alpha * t(IRF(parms, theta) * (1-IRF(parms, theta))))
}

#--------------------------------------------------------------------------
#' Second deriviate of 2PL in theta
#'
#' Used for obtaining SEs of theta and computing WMLE.
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of second derivatives of IRF
#' @export

d2IRF <- function(parms, theta) {
  t(parms$alpha^2 * t(dIRF(parms, theta) * (1 - 2 * IRF(parms, theta))))
}

#--------------------------------------------------------------------------
#' Test information function for IRF
#'
#' Computes test information under 2PL model
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)} vector of test information
#' @export

I <- function(parms, theta) {
  p <- IRF(parms, theta)
  q <- 1 - p
  temp <-  t(parms$alpha^2 * t(p * q))
  # temp <- dIRF(parms, theta)^2 / p / q  # too slow with 2PL
  apply(temp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Numerator of penalty used in WMLE
#'
#' Computes equation 15 of Magis & Raiche (2012). For 2PL, this is also the derivtive of the item info
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait
#' @param theta2 not used; included for formal compatabiity with collaborative IRFs
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

J <- function(parms, theta) {
  p <- IRF(parms, theta)
  q <- 1 - p
  temp <- t(parms$alpha^3 * t(p * q * (1 - 2 * p)))
  # temp <- dIRF(parms, theta) * d2IRF(parms, theta) / p / q # too slow with 2PL
  apply(temp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Derviative of J
#'
#' Computes the derivative of J (for 2PL, second derivative of item info)
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait
#' @param theta2 not used; included for formal compatabiity with collaborative IRFs
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

dJ <- function(parms, theta) {
  p <- IRF(parms, theta)
  q <- 1 - p
  temp <- t(parms$alpha^4 * t(p * q * (1 - 6 * p * q) ))
  apply(temp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Standard error of WMLE
#'
#' Computes equation 16 of Magis & Raiche (2012)
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait
#' @param theta2 not used; included for formal compatabiity with collaborative IRFs
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

  SE <- function(parms, theta) {
    temp <- J(parms, theta)^2 + I(parms, theta) * dJ(parms, theta)
    temp <- temp / 2 / I(parms, theta)^2 + I(parms, theta)
    1/sqrt(temp)
  }

#--------------------------------------------------------------------------
#' Loglikelihood of 2PL in theta
#'
#' Computes likeihood of a response pattern under the 2PL model, given item parms and theta
#'
#' @param resp the matrix binary responses
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)}-vecor of (possibly weigthed) loglikelihoods
#' @export

logL <- function(resp, parms, theta) {
  p <- IRF(parms, theta)
  apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Derviative of loglikelihood of 2PL in theta
#'
#' Computes either deriviat of standard MLE or of WMLE
#'
#' @param resp the matrix binary responses
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @param WMLE use Weighted ML or not?
#' @return \code{length(theta)}-vecor of (possibly weigthed) loglikelihoods
#' @export

dlogL <- function(theta, resp, parms, WMLE = T) {
  p <- IRF(parms, theta)
  dlog <- apply(parms$alpha * t(resp - p), 2, sum, na.rm = T)
  if (WMLE) {
    dlog <- dlog + J(parms, theta) / 2 / I(parms, theta)
  }
  dlog
}

#--------------------------------------------------------------------------
#' ML or WML estimation of latent trait for 2PL model
#'
#' Theta is estimated by calling \code{uniroot} on \code{dlogL}. SEs are computed analytically. Value of the loglikelihood at the estimate is also provided. If \coe{parallel = T}, the call to uniroot is parallelized via \code{mclapply}. (ToDo: remove duplicate response patterns; also could do with a better way of setting the range of x in uniroot)
#''
#' @param resp the matrix binary responses
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param parallel call parallel:mclapply instead of looping over \code{nrow(resp)}?
#' @return An \code{nrow(resp)} by 3 matrix of log-likleihoods, theta estimates, analytic SEs for each response pattern in resp.
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
    out$se <- SE(parms, out$theta)
  } else {
    out$se <- 1/sqrt(I(parms, out$theta))
  }
  out$logL <- logL(resp, parms, out$theta)
  out
}


# # XX attempt to deal with duplicate response patterns, to shorten runtine for MLE
#
# get_duplicates <- function(resp) {
#   doubles <- which(duplicated(resp) | duplicated (resp, fromLast = T))
#   if (length(doubles) > 0) {
#     temp <- lapply(data.frame(t(resp[doubles,])), paste0, collapse = "") %>% unlist
#     sort_ind <- order(temp)
#     temp <- temp[sort_ind]
#     doubles <- doubles[sort_ind]
#     key_ind <- cumsum(!duplicated(temp))
#     uniques <- doubles[!duplicated(temp)]
#     originals <- uniques[key_ind]
#     out <- cbind(doubles, originals)
#     out[duplicated(out[,2]),]
#   } else {
#     NA
#   }
# }
