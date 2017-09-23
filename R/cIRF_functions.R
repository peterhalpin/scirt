# Last update: 04/25/2017
# Functions for group assessments. Depends on IRF_functions.R
# User beware: functions do not check or handle input errors.

require(ggplot2)
require(dplyr)

#--------------------------------------------------------------------------
#' The IRF of the one-parameter RSC model.
#'
#' Computes a matrix of probabilities for correct responses under one parameter RSC model, using 2PL model for individual responses.
#'
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @return  \code{length(w)} by \code{nrow(parms)} matrix of response probabilities.
#' @export

RSC <- function(w, parms, theta1, theta2) {
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  W <- w %*% t(rep(1, nrow(parms)))
  W * (p1 + p2) + (1 - 2 * W) * p1 * p2
}



#--------------------------------------------------------------------------
#' First derivatives of the IRF of the one-parameter RSC model.
#'
#' Computes a named list of derivatives of the IRF of the one-parameter RSC model. Uses 2PL model for individual responses. See Appendix of Halpin & Bergner (2017) for details.
#'
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @return a named list of derivatives, with elements \code{c("dw", "dtheta1", "dtheta2")}.
#' @export

d_RSC <- function(w, parms, theta1, theta2) {
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  dp1 <- dIRF(parms, theta1)
  dp2 <- dIRF(parms, theta2)
  W <- w %*% t(rep(1, nrow(parms)))
  dw <- p1 + p2 - 2 * p1 * p2
  dt1 <- (W + (1 - 2 * W) * p2) * dp1
  dt2 <- (W + (1 - 2 * W) * p1) * dp2
  out <- list(dw, dt1, dt2)
  names(out) <- c("dw", "dtheta1", "dtheta2")
  out
}

#--------------------------------------------------------------------------
#' Second derivatives of the IRF of the one-parameter RSC model.
#'
#' Computes a named list of derivatives of the IRF of the one-parameter RSC model. Uses 2PL model for individual responses. See Appendix of Halpin & Bergner (2017) for details.
#'
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @return a named list of second derivatives, with elements \code{c("dw_dtheta1", "dw_dtheta2", "d2theta1", "d2theta2", "dtheta1_dtheta2")}.
#' @export

d2_RSC <- function(w, parms, theta1, theta2) {
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  dp1 <- dIRF(parms, theta1)
  dp2 <- dIRF(parms, theta2)
  d2p1 <- d2IRF(parms, theta1)
  d2p2 <- d2IRF(parms, theta2)
  W <- w %*% t(rep(1, nrow(parms)))

  dwdt1 <- dp1 * (1 - 2 * p2)
  dwdt2 <- dp2 * (1 - 2 * p1)
  d2t1 <- (W + (1 - 2 * W) * p2) * d2p1
  d2t2 <- (W + (1 - 2 * W) * p1) * d2p2
  dt1dt2 <- (1 - 2 * W) * dp1 * dp2
  out <- list(dwdt1, dwdt2, d2t1, d2t2, dt1dt2)
  names(out) <- c("dw_dtheta1", "dw_dtheta2", "d2theta1", "d2theta2", "dtheta1_dtheta2")
  out
}

#--------------------------------------------------------------------------
#' Log-likelihood of the one-parameter RSC model.
#'
#' Computes loglikeihood of a matrix of response patterns under the one-parameter RSC model.
#'
#' @param resp a matrix or data.frame containing the (conjunctively-scored) binary item responses.
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @return  \code{length(w)}-vector of log-likelihood.
#' @export

l_RSC <- function(resp, w, parms, theta1, theta2) {
  p <- RSC(w, parms, theta1, theta2)
  apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Mutiplier for first deriviative of log-likelihood of one-parameter RSC model.
#'
#' Writing the derivative of the log-likelihood of the RSC model a single item as M * grad_IRF, this function computes (a matrix of values of) M, with grad_IRF being the gradient of the RSC IRF. Called by functions that compute derivative of the log-likelihood.

#' @param resp a matrix or data.frame containing the (conjunctively-scored) binary item responses.
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @return \code{dim(resp)}-matrix of multipliers.
#' @export

Mstar <- function(resp, w, parms, theta1, theta2) {
  p <- RSC(w, parms, theta1, theta2)
  resp / p - (1 - resp) / (1 - p)
}



#--------------------------------------------------------------------------
#' Gradient of log-likelihood of one-parameter RSC model.
#'
#' Computes the first derivatives of the log-likelihood, in \code{c(w, theta1, theta2)}. Called by \code{est_RSC}.

#' @param resp a matrix or data.frame containing the (conjunctively-scored) binary item responses.
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @return \code{3 * K} - vector of first derivatives, with \code{K = length(w)}, ordered as \code{rep(c(w_k, theta1_k, theta2_k), times = K)}.
#' @export

dl_RSC <- function(resp, w, parms, theta1, theta2) {
  m <- Mstar(resp, w, parms, theta1, theta2)
  d <- d_RSC(w, parms, theta1, theta2)
  dw <- apply(m * d$dw, 1, sum, na.rm = T)
  dt1 <- apply(m * d$dtheta1, 1, sum, na.rm = T)
  dt2 <- apply(m * d$dtheta2, 1, sum, na.rm = T)
  rbind(dw, dt1, dt2) %>% c
}

#--------------------------------------------------------------------------
#' Multiplier for second deriviative of log-likelihood of one-parameter RSC model.
#'
#' Similar to \code{Mstar}, but for the second derivatives. Called by functions that compute second derivative of the log-likelihood.

#' @param resp a matrix or data.frame containing the (conjunctively-scored) binary item responses.
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @param obs logical: should the observed value be returned? If not, the expected value is returned.
#' @return \code{dim(resp)}-matrix of multipliers.
#' @export

Nstar <- function(resp, w, parms, theta1, theta2, obs = F) {
  p <- RSC(w, parms, theta1, theta2)
  if (obs) {
    resp / p^2 + (1 - resp) / (1 - p)^2
  } else {
    1 / p / (1 - p) #* !is.na(resp)
  }
}

#--------------------------------------------------------------------------
#' Hessian of log-likelihood of one-parameter RSC model.
#'
#' Computes the second derivatives of the log-likelihood, in \code{c(w, theta1, theta2)}. Called by \code{est_RSC}. Calls the function \code{bdiag_m} written by Martin Maechler, ETH Zurich.
#'
#' @param resp a matrix or data.frame containing the (conjunctively-scored) binary item responses.
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @param obs logical: should the observed value be returned? If not, the expected value is returned.
#' @return \code{3 * K} by \code{3 * K} block-diagonal, symmmetrical matrix of second derivatives, with \code{K = length(w)} and rows/cols ordered as \code{rep(c(w_k, theta1_k, theta2_k), times = K)}.
#' @export

d2l_RSC <- function(resp, w, parms, theta1, theta2, obs = T) {
  n <- -1 * Nstar(resp, w, parms, theta1, theta2, obs)
  d <- d_RSC(w, parms, theta1, theta2)
  dwdw <-  apply(n * d$dw * d$dw, 1, sum, na.rm = T)
  dwdt1 <- apply(n * d$dw * d$dtheta1, 1, sum, na.rm = T)
  dwdt2 <- apply(n * d$dw * d$dtheta2, 1, sum, na.rm = T)
  dt1dt1 <- apply(n * d$dtheta1 * d$dtheta1, 1, sum, na.rm = T)
  dt1dt2 <- apply(n * d$dtheta1 * d$dtheta2, 1, sum, na.rm = T)
  dt2dt2 <- apply(n * d$dtheta2 * d$dtheta2, 1, sum, na.rm = T)

  if (obs) {
    m <- Mstar(resp, w, parms, theta1, theta2)
    d2 <- d2_RSC(w, parms, theta1, theta2)
    dwdt1 <- dwdt1 + apply(m * d2$dw_dtheta1, 1, sum, na.rm = T)
    dwdt2 <- dwdt2 + apply(m * d2$dw_dtheta2, 1, sum, na.rm = T)
    dt1dt1 <- dt1dt1 + apply(m * d2$d2theta1, 1, sum, na.rm = T)
    dt1dt2 <- dt1dt2 + apply(m * d2$dtheta1_dtheta2, 1, sum, na.rm = T)
    dt2dt2 <- dt2dt2 + apply(m * d2$d2theta2, 1, sum, na.rm = T)
  }

  fun <- function(i) {
      temp <- c(dwdw[i], dwdt1[i], dwdt2[i], dwdt1[i], dt1dt1[i], dt1dt2[i], dwdt2[i], dt1dt2[i], dt2dt2[i])
      matrix(temp, nrow = 3, ncol = 3)
  }

  #if (parallel) { # This isnt speeding anything up
  #  temp <- parallel::mclapply(1:length(w), fun)
  #} else {
    temp <- vector("list", length(theta1))
    for (i in 1:length(theta1)) {temp[[i]] <- fun(i)}
  #}
  bdiag_m(temp)
}


#--------------------------------------------------------------------------
#' Log-likelihood of a combined assessment.
#'
#' This function computes the log-likelihood for a combined assessment, in which the 2PL model is used for the individual component of the assessment and the one-parameter RSC model is used for the group component of the assessment. The derivatives are taken in \code{c(w, theta1, theta2)}.
#'
#' The response matrix \code{resp} must be formatted to contain one row of binary responses for each respondent (not each dyad). Members of the same dyad must be on adjancent rows, such that \code{resp[odd,]} gives the responses of one member of a dyad and \code{resp[odd + 1, ]} gives the responses of the other member of the dyad, where \code{odd} is any odd integer in \code{c(1, nrow(resp))}. The (column) names for items on the individual assessment must include \code{"IND"}; those on the (conjunctively-scored) group assessment just include \code{"COL"} -- these text-keys are grepped from \code{names(resp)} to obtain the response patterns for the individual assessment and the group assessment. Note that only the odd rows of \code{resp[grep("COL", names(resp))]} are used when computing the log-likelihood for the group component.

#' The order of items (columns) of \code{resp} is assumed to correpond to that of items (rows) of \code{parms}, for each of \code{c("IND", "COL")}. Similarly to the procedure described for \code{names(resp)}, \code{row.names(parms)} is grepped for each of \code{c("IND", "COL")} to obtain the item parameters of the individual assessment and the group assessment.
#'
#' This description is much longer than the source code -- type \code{l_full} for a shorter explanation.

#' @param resp a data.frame containing the binary item responses of both the individual assessment and the (conjunctively scored) group assessment. See details for information on formatting.
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively. See details for information on formatting.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @return \code{length(w)}-vector of log-likelihoods.
#' @export

l_full <- function(resp, w, parms, theta1, theta2)
{
  ind <- grep("IND", names(resp))
  col <- grep("COL", names(resp))
  ind_parms <- parms[grep("IND", row.names(parms)),]
  col_parms <- parms[grep("COL", row.names(parms)),]
  odd <- seq(1, nrow(resp), by = 2)
  (logL(resp[odd, ind], ind_parms, theta1) +
    logL(resp[(odd + 1), ind], ind_parms, theta2) +
    l_RSC(resp[odd, col], w, col_parms, theta1, theta2)) %>% sum
}

#--------------------------------------------------------------------------
#' Log-likelihood of a combined assessment, with sum option.
#'
#' This function is identical to \code{l_full}, except that it has an option for whether or not to sum the log-liklihood over respondents. See \code{help(l_full)} for details on formatting \code{resp} and \code{parms}.
#'
#' @param resp a data.frame containing the binary item responses of both the individual assessment and the (conjunctively scored) group assessment. See details for information on formatting.
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively. See details for information on formatting.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @return \code{length(w)}-vector of log-likelihoods.
#' @export

l_full_sum <- function(resp, w, parms, theta1, theta2, Sum = F)
{
  ind <- grep("IND", names(resp))
  col <- grep("COL", names(resp))
  ind_parms <- parms[grep("IND", row.names(parms)),]
  col_parms <- parms[grep("COL", row.names(parms)),]
  odd <- seq(1, nrow(resp), by = 2)
  temp <- logL(resp[odd, ind], ind_parms, theta1) +
    logL(resp[(odd + 1), ind], ind_parms, theta2) +
    l_RSC(resp[odd, col], w, col_parms, theta1, theta2)
  if (Sum) {sum(temp)} else {temp}
}

#--------------------------------------------------------------------------
#' Log of prior distribution for a combined assessment.
#'
#' This function computes the log of the prior (minus a constant) for a combined assessment, in which the 2PL model is used for the individual component of the assessment and the one-parameter RSC model is used for the group component of the assessment. A standard normal prior is used for individual ability. A two-parameter Beta prior is the parameter of the RSC model, in which both parameters are equal to 1 + epsilon.
#'

#' @param w the weight parameter of the RSC model.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @param epsilon a small positive number, see description for details.
#' @return \code{length(w)}-vector of log-priors (minus a constant).
#' @export

lp <- function(w, theta1, theta2, epsilon = .05)
{
  (epsilon * log(w - w^2) - theta1^2 / 2 - theta2^2 / 2)
}

#--------------------------------------------------------------------------
#' Gradient of the log-likelihood of a combined assessment.
#'
#' This function computes first derivatives of the log-likelihood for a combined assessment, in which the 2PL model is used for the individual component of the assessment and the one-parameter RSC model is used for the group component of the assessment. The derivatives are taken in \code{c(w, theta1, theta2)}. See \code{help(l_full)} for details on formatting \code{resp} and \code{parms}.
#'
#' @param resp a data.frame containing the binary item responses of both the individual assessment and the (conjunctively scored) group assessment. See details for information on formatting.
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively. See details for information on formatting.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @return \code{3 * K} - vector of first derivatives, with \code{K = length(w)}, ordered as \code{rep(c(w_k, theta1_k, theta2_k), times = K)}.
#' @export

dl_full <- function(resp, w, parms, theta1, theta2) {
  ind <- grep("IND", names(resp))
  col <- grep("COL", names(resp))
  ind_parms <- parms[grep("IND", row.names(parms)),]
  col_parms <- parms[grep("COL", row.names(parms)),]
  odd <- seq(1, nrow(resp), by = 2)

  dl_t1 <- dlogL(resp[odd, ind], ind_parms, theta1)
  dl_t2 <- dlogL(resp[(odd+1), ind], ind_parms, theta2)
  temp <- rbind(0, dl_t1, dl_t2) %>% c
  dl_RSC(resp[odd, col], w, col_parms, theta1, theta2) + temp
}

#--------------------------------------------------------------------------
#' Graident of the log of the prior distribution for a combined assessment.
#'
#' This function computes the derivative of the log of the prior for a combined assessment, in which the 2PL model is used for the individual component of the assessment and the one-parameter RSC model is used for the group component of the assessment. A standard normal prior is used for individual ability. A two-parameter Beta prior is the parameter of the RSC model, in which both parameters are equal to 1 + epsilon.
#'
#'
#' @param w the weight parameter of the RSC model.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @param epsilon a small positive number, see description for details.
#' @return \code{3 * K} - vector of first derivatives, with \code{K = length(w)}, ordered as \code{rep(c(w_k, theta1_k, theta2_k), times = K)}.
#' @export

dlp <- function(w, theta1, theta2, epsilon = .05)
{
  dlw <- epsilon * (1 - 2 * w) / (w - w^2)
  rbind(dlw, -theta1, -theta2) %>% c
}

#--------------------------------------------------------------------------
#' Hessian of the log-likelihood for a combined assessment.
#'
#' This function computes second derivatives of the log-likelihood for a combined assessment, in which the 2PL model is used for the individual component of the assessment and the one-parameter RSC model is used for the group component of the assessment. The derivatives are taken in \code{c(w, theta1, theta2)}. See \code{help(l_full)} for details on formatting \code{resp} and \code{parms}.
#'
#' @param resp a data.frame containing the binary item responses of both the individual assessment and the (conjunctively scored) group assessment. See details for information on formatting.
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively. See details for information on formatting.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @param obs logical: should the observed value be returned? If not, the expected value is returned.
#' @return \code{3 * K} by \code{3 * K} block-diagonal, symmmetrical matrix of second derivatives, with \code{K = length(w)} and rows/cols ordered as \code{rep(c(w_k, theta1_k, theta2_k), times = K)}.
#' @export

d2l_full <- function(resp, w, parms, theta1, theta2, obs = T, parallel = F) {
  ind <- grep("IND", names(resp))
  col <- grep("COL", names(resp))
  ind_parms <- parms[grep("IND", row.names(parms)),]
  col_parms <- parms[grep("COL", row.names(parms)),]
  odd <- seq(1, nrow(resp), by = 2)

  d2l_t1 <- d2logL(resp[odd, ind], ind_parms, theta1, obs)
  d2l_t2 <- d2logL(resp[(odd+1), ind], ind_parms, theta2, obs)
  temp <- rbind(0, d2l_t1, d2l_t2) %>% c %>% diag
  d2l_RSC(resp[odd, col], w, col_parms, theta1, theta2, obs) + temp
}

#--------------------------------------------------------------------------
#' Hessian of the log of the prior distribution for a combined assessment.
#'
#' This function computes the Hessian of the log of the prior for a combined assessment, in which the 2PL model is used for the individual component of the assessment and the one-parameter RSC model is used for the group component of the assessment. A standard normal prior is used for individual ability. A two-parameter Beta prior is the parameter of the RSC model, in which both parameters are equal to 1 + epsilon.
#'
#' @param w the weight parameter of the RSC model.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @param epsilon a small positive number, see description for details.
#' @return \code{3 * K} by \code{3 * K} diagonal matrix of second derivatives, with \code{K = length(w)} and rows/cols ordered as \code{rep(c(w_k, theta1_k, theta2_k), times = K)}.
#' @export

d2lp <- function(w, theta1, theta2, epsilon = .05)
{
  w2 <- w - w^2
  d2lw <- -1 * epsilon * (2 / w2 + ((1 - 2 * w) / w2)^2 )
  rbind(d2lw, -1, -1) %>% c %>% diag
}

#--------------------------------------------------------------------------
#' Simultaneous estimation of latent traits and the one-parameter RSC model for a combined assessment.
#'
#' This function calls \code{optim} to estimate the parameter vector \code{c(w, theta1, theta2)} from the repsonses to a combined assessment, in which the 2PL model is used for the individual component of the assessment and the one-parameter RSC model is used for the group component of the assessment.
#'
#' @details Esimation is via either maximum likelihood (ML) or modal a'posteriori (MAP), with the latter being prefered. For MAP, a standard normal prior is used for individual ability. A two-parameter Beta prior is the parameter of the RSC model, in which both parameters are equal to 1 + \code{epsilon}. Standard errors (or posterior standard deviations) are computed by numerically inverting the analytically computed Hessian of the objective function, at the parameter estimates. The value of the objective function at the estimate is is also provided. If \code{parallel = T}, the call to \code{optim} is parallelized via \code{parallel::mclapply}.
#'
#' The response matrix \code{resp} must be formatted to contain one row of binary responses for each respondent (not each dyad). Members of the same dyad must be on adjancent rows, such that \code{resp[odd,]} gives the responses of one member of a dyad and \code{resp[odd + 1, ]} gives the responses of the other member of the dyad, where \code{odd} is any odd integer in \code{c(1, nrow(resp))}. The (column) names for items on the individual assessment must include \code{"IND"}; those on the (conjunctively-scored) group assessment just include \code{"COL"} -- these text-keys are grepped from \code{names(resp)} to obtain the response patterns for the individual assessment and the group assessment. Note that only the odd rows of \code{resp[grep("COL", names(resp))]} are used when computing the log-likelihood for the group component.
#'
#' The order of items (columns) of \code{resp} is assumed to correpond to that of items (rows) of \code{parms}, for each of \code{c("IND", "COL")}. Similarly to the procedure described for \code{names(resp)}, \code{row.names(parms)} is grepped for each of \code{c("IND", "COL")} to obtain the item parameters of the individual assessment and the group assessment.
#' Type \code{l_full} for an illustration of how the formatting calls are made.
#'
#' @param resp a data.frame containing the binary item responses of both the individual assessment and the (conjunctively scored) group assessment. See details for information on formatting.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively. See details for information on formatting.
#' @param starts starting values, ordered as triplets of \code{c(w, theta1, theta2)} for each row or \code{resp} (optional).
#' @param method one of \code{c("ML", "MAP")}. The latter is strongly recommended.
#' @param obs logical: should standard errors be computed using the observed (\code{TRUE}) or expected (\code{FALSE}) Hessian?
#' @param epsilon a small positive number, see description for details.
#' @param parallel logical: call \code{parallel:mclapply} instead of looping over \code{nrow(resp)}?
#' @return An named \code{nrow(resp)} by 7 data.frame containing the estimates, their standard errors, and the value of the objective function at the solution.
#' @export

est_RSC <- function(resp, parms, starts = NULL, method = "MAP", obs = F, epsilon = .05, parallel = F) {

  stopifnot(ncol(resp) == nrow(parms),
          method%in%c("MAP", "ML"),
          is.logical(obs),
          is.logical(parallel))

  n_obs <- nrow(resp)/2
  odd <- seq(from = 1, to = n_obs*2, by = 2)
  parm_index <- seq(from = 1, to = n_obs*3, by = 3)
  out <- data.frame(matrix(0, nrow = n_obs, ncol = 7))
  names(out) <- c("log", "w", "w_se", "theta1", "theta1_se", "theta2", "theta2_se")


  # Starting values
  if(is.null(starts)) {starts <- rep(c(.5, 0, 0), times = n_obs)}
  lower <- rep(c(.00001, -8, -8), times = n_obs)
  upper <- rep(c(.99999, 8, 8), times = n_obs)

  # Select objective function and gradient
  if (method == "ML"){
    obj <- function(par, resp) {
      ind <- parm_index[1:(length(par) / 3)]
      -1 * l_full(resp, par[ind], parms, par[(ind+1)], par[(ind+2)])
    }
    grad <- function(par, resp) {
      ind <- parm_index[1:(length(par) / 3)]
      -1 * dl_full(resp, par[ind], parms, par[(ind+1)], par[(ind+2)])
    }
  }

  if (method == "MAP"){
    obj <- function(par, resp) {
      ind <- parm_index[1:(length(par) / 3)]
      -1 * l_full(resp, par[ind], parms, par[(ind+1)], par[(ind+2)]) -
        sum(lp(par[ind], par[(ind+1)], par[(ind+2)], epsilon))
    }
    grad <- function(par, resp) {
      ind <- parm_index[1:(length(par) / 3)]
      -1 * dl_full(resp, par[ind], parms, par[(ind+1)], par[(ind+2)]) -
        dlp(par[ind], par[(ind+1)], par[(ind+2)], epsilon)
    }
  }

  # Set up parm_indexing for parallel
  fun <- function(i) {
    blocksize <- floor(n_obs/n_cores)
    m <- (i-1) * blocksize + 1
    if (i < n_cores) {n <- i * blocksize} else {n <- n_obs}
    ind1 <- parm_index[m] : (parm_index[n] + 2)
    ind2 <- odd[m] : (odd[n] + 1)
    q <- optim(starts[ind1], obj,
             gr = grad,
             resp = resp[ind2,],
             method = "L-BFGS-B",
             lower = lower[ind1],
             upper = upper[ind1]
             )$par
  }

  # Estimation
  n_cores <- parallel::detectCores()
  if (parallel & n_cores < n_obs) {
    temp <- parallel::mclapply(1:n_cores, fun) %>% unlist
  } else {
    n_cores <- 1
    temp <- fun(1)
  }
  out$w <- temp[parm_index]
  out$theta1 <- temp[parm_index+1]
  out$theta2 <- temp[parm_index+2]

  # Standard errors
  temp_se <- d2l_full(resp, out$w, parms, out$theta1, out$theta2, obs)
  if (method == "MAP") {
    temp_se <- temp_se + d2lp(out$w, out$theta1, out$theta2, epsilon)
  }
  temp_se <- (-1 * temp_se) %>% solve %>% diag %>% sqrt

  out$w_se <- temp_se[parm_index]
  out$theta1_se <- temp_se[parm_index+1]
  out$theta2_se <- temp_se[parm_index+2]

  # Objective function
  out$log <- l_full_sum(resp, out$w, parms, out$theta1, out$theta2)
  if (method == "MAP") {
    out$log <- out$log + lp(out$w, out$theta1, out$theta2, epsilon)
  }

  out
}

#--------------------------------------------------------------------------
#' Estimation of the one-parameter RSC model, with latent traits assumed to be known.
#'
#' This function calls \code{optim} to estimate the one-parameter RSC model from the (conjunctively-scored) repsonses of dyads to a group assessment.
#'
#' @details Estimation is via either maximum likelihood (ML) or modal a'posteriori (MAP), with the latter being prefered. For MAP, a two-parameter Beta prior is used with the parameter of the RSC model, in which both parameters are equal to \code{1 + epsilon}. Standard errors (or posterior standard deviations) are computed via the inverse of the analytically computed second derivatives of the objective function, at the parameter estimates. The value of the objective function at the estimate is is also provided. If \code{parallel = T}, the call to \code{optim} is parallelized via \code{parallel::mclapply}.
#'
#' @param resp a matrix or data.frame containing the (conjunctively-scored) binary item responses.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait of member 1.
#' @param theta2 the latent trait of member 2.
#' @param method one of \code{c("ML", "MAP")}. The latter is strongly recommended.
#' @param obs logical: should standard errors be computed using the observed (\code{TRUE}) or expected (\code{FALSE}) Fisher information?
#' @param epsilon a small positive number, see description for details.
#' @param parallel logical: call \code{parallel:mclapply} instead of looping over \code{nrow(resp)}?
#' @return An named \code{nrow(resp)} by 3 data.frame containing the estimates, their standard errors, and the value of the log-likelihood of the RSC model at the solution (not log posterior with MAP).
#' @export

est_RSC2 <- function(resp, parms, theta1, theta2, method = "MAP", obs = F, epsilon = .05, parallel = F) {

  stopifnot(ncol(resp) == nrow(parms),
            nrow(resp) == length(theta1),
            length(theta1) == length(theta2),
            method%in%c("MAP", "ML"),
            is.logical(obs),
            is.logical(parallel))

  n_obs <- nrow(resp)
  odd <- seq(from = 1, to = n_obs, by = 2)
  parm_index <- seq(from = 1, to = n_obs*3, by = 3)
  out <- data.frame(matrix(0, nrow = n_obs, ncol = 3))
  names(out) <- c("log", "w", "se")
  starts <- rep(.5, times = n_obs)

  # Select objective function and gradient
  fun1 <- function(par, resp, theta1, theta2) {
    -1 * sum(l_RSC(resp, par, parms, theta1, theta2))
  }
  fun2 <- function(par, resp, theta1, theta2) {
    -1 * dl_RSC(resp, par, parms, theta1, theta2)[parm_index[1:length(par)]]
  }

  if (method == "ML") {
    obj <- function(par, resp, theta1, theta2) {fun1(par, resp, theta1, theta2)}
    grad <- function(par, resp, theta1, theta2) {fun2(par, resp, theta1, theta2)}
  }

  if (method == "MAP") {
    obj <- function(par, resp, theta1, theta2) {
      fun1(par, resp, theta1, theta2) - sum(epsilon * log(par - par^2))
    }
    grad <- function(par, resp, theta1, theta2) {
      fun2(par, resp, theta1, theta2) - epsilon * (1 - 2 * par) / (par - par^2)
    }
  }
 
   # Set up parm_indexing for parallel
  fun <- function(i) {
    blocksize <- floor(n_obs/n_cores)
    m <- (i-1) * blocksize + 1
    if (i < n_cores) {n <- i * blocksize} else {n <- n_obs}
    q <- optim(starts[m:n], obj,
             gr = grad,
             resp = resp[m:n, ],
             theta1 = theta1[m:n],
             theta2 = theta2[m:n],
             method = "L-BFGS-B",
             lower = .000001,
             upper = .999999)$par
  }

  # Estimation
  if(parallel) {
    n_cores <- parallel::detectCores()
    out$w <- parallel::mclapply(1:n_cores, fun) %>% unlist
  } else {
    n_cores <- 1
    out$w <- fun(1)
  }

  # Standard errors
  out$se <- diag(d2l_RSC(resp, out$w, parms, theta1, theta2, obs))[parm_index]
  if (method == "MAP") {
    out$se <- out$se + diag(d2lp(out$w, theta1, theta2, epsilon))[parm_index]
  }
  out$se <- sqrt(-1 / out$se)

  # Objective function
  out$log <- l_RSC(resp, out$w, parms, theta1, theta2)
  #if (method == "MAP") {
  #  out$log <- out$log + epsilon * log(out$w - out$w^2)
  #}
  out
}


#--------------------------------------------------------------------------
#' Simulate item responses from the one-parameter RSC model.
#'
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively. See details for information on formatting.
#' @param theta1 the latent trait of member 1.
#' @param theta1 the latent trait of member 2.
#' @return \code{length(theta1)} by \code{nrow(parms)} matrix of binary responses.
#' @export

sim_RSC <- function(w, parms, theta1 = 0, theta2 = 0) {
  n_row <- length(theta1)
  n_col <- nrow(parms)
  r <- array(runif(n_row * n_col), dim = c(n_row, n_col))
  p <- RSC(w, parms, theta1, theta2)
  out <- ifelse(p > r, 1, 0)
  colnames(out) <- row.names(parms)
  out
}

#--------------------------------------------------------------------------
#' Generate data from the one-parameter RSC model.
#'
#' This is a wrapper for \code{sim_RSC} that saves the data generating parms and allows for multple response patterns per group.
#'
#' To generate data from a 2PL model, set \code{theta1 = theta2} and \code{w = 1/2}. See Halpin and Bergner (2017) for discussion.
#'
#' @param n_reps integer indicating how many datasets to generate.
#' @param w the weight parameter of the RSC model.
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait for member 1.
#' @param theta2 the latent trait for member 2.
#' @param theta1_se the standard error of the latent trait for member 1 (Optional). If included, data generation samples \code{n_reps} values from \code{rnorm(theta1, theta1_se)}.
#' @param theta2_se the standard error of the latent trait for member 2 (Optional). If included, data generation samples \code{n_reps} values from \code{rnorm(theta2, theta2_se)}.
#' @param NA_pattern an (optional) \code{length(w)} by \code{nrow(parms)} data.frame with \code{NA} entries denoting missing data. The missing values are preserved in the generated data.

#' @return A data.frame with \code{length(w)} rows containing an id variable for each pair and each sample, the data generating values of \code{w}, \code{theta1}, and \code{theta2}, and the simulated response patterns.
#' @export

data_gen <- function(n_reps, w, parms, theta1, theta2, theta1_se = NULL, theta2_se = NULL, NA_pattern = NULL) {

  # Storage
  n_obs <- length(theta1)
  out <- data.frame(rep(1:n_obs, each = n_reps), rep(1:n_reps, times = n_obs))
  names(out) <- c("pairs", "samples")
  out$w <- rep(w, each = n_reps)
  out$theta1 <- theta_gen(n_reps, theta1, theta1_se)
  out$theta2 <- theta_gen(n_reps, theta2, theta2_se)

  # Simulate data
  data <- data.frame(matrix(NA, nrow = n_reps*n_obs, ncol = nrow(parms)))
  names(data) <- row.names(parms)
  data <- sim_RSC(out$w, parms, out$theta1, out$theta2)
  data <- format_NA(data, NA_pattern)
  cbind(out[], data[])
}



#--------------------------------------------------------------------------
#' Formats a data.frame to have the same variables (column names) as that of another data.frame
#'
#' Drops columns in the target data.frame (resp) not in the variable name list (items). Adds columns (with NA entries) for entries in the items that are not alrady in resp. This is mainly to simplify using item parameters obtained from one sample (e.g., calibration sample) with a second sample that may not use all of the items.

#' @param resp the target data.frame to be re-formatted.
#' @param item string vector of variable names.
#' @param version an optional string used to subset resp via \code{grep(version, names(resp))}
#' @return a data.frame that results from padding and deleting \code{resp} as described.
#' @export

format_resp <- function(resp, items, version = NULL) {
  if (!is.null(version)) {
    temp <- resp[grep(version, names(resp))]
  }
  temp[items[!items%in%names(temp)]] <- NA
  temp <- temp[items]
  temp
}


#--------------------------------------------------------------------------
#' Transforms one matrix/data.frame have the same NA entries as another matrix/data.frame
#'
#' Useful for reproducing observed NA patterns in simulated data. If \code{q = (nrow(data) / nrow(NA_pattern)} is an integer greater than 1, will replicate NA_pattern as \code{kronecker(as.matrix(NA_pattern), rep(1, q)}.
#'
#' @param data the matrix/data.frame into which to write NA entries.
#' @param NA_pattern the matrix/data.frame from which to write NA entries.

#' @return a data.frame/matrix similar to \code{data} that includes all the NA entires as \code{NA_pattern}. If \code{NA_pattern == NULL}, the data.frame/matrix is returned unaltered.
#' @export

format_NA <- function(data, NA_pattern = NULL){
  if (!is.null(NA_pattern)) {
    n_reps = nrow(data) / nrow(NA_pattern)
    temp <- kronecker(as.matrix(NA_pattern), rep(1, n_reps))
    temp[!is.na(temp)] <- 1
    data <- data*temp
 }
 data
}


#--------------------------------------------------------------------------
#' Generates or replicates values of latent variable.
#'
#' If \code{theta_se}' is not null, values are generated using \code{rnorm}. Otherwise, each value is replicated.

#' @param n number of values to generate / replicate for each value of \code{theta}.
#' @param theta the latent trait.
#' @param theta_se the standard error of the latent trait.

#' @return An \code{n * length(theta)} vector in which each value \code{theta} is generated / replicated \code{n} times.
#' @export

theta_gen <- function(n, theta, theta_se = NULL){
  theta_long <- rep(theta, each = n)

  if(!is.null(theta_se)) {
    theta_long <- rnorm(length(theta_long), theta_long, rep(theta_se, each = n))
  }
  theta_long
}


bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
    (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
    all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
    ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
    i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
    p = k * 0L:M,
    x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}
