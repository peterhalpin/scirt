
require(Matrix)

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
#' First deriviate of 2PL IRF in theta
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
#' Second derivative of 2PL IRF in theta
#'
#' Used for obtaining SEs of theta and computing WMLE.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of second derivatives of IRF
#' @export

d2IRF <- function(parms, theta) {
  t(parms$alpha * t(dIRF(parms, theta) * (1 - 2 * IRF(parms, theta))))
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
#' First derviative of loglikelihood in theta
#'
#' @param resp a matrix or data.frame containing the binary item responses
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @param WMLE logical: use weighted ML or not?
#' @return \code{length(theta)}-vecor of (possibly weigthed) loglikelihoods
#' @export

dlogL <- function(resp, parms, theta) {
  m <- M(resp, parms, theta)
  dp <- dIRF(parms, theta)
  apply(m * dp, 1, sum, na.rm = T)
}

#--------------------------------------------------------------------------
#' Second derviative of loglikelihood in theta
#'
#' @param resp a matrix or data.frame containing the binary item responses
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta the latent trait
#' @param WMLE logical: use weighted ML or not?
#' @return \code{length(theta)}-vecor of (possibly weigthed) loglikelihoods
#' @export

d2logL <- function(resp, parms, theta, obs = T) {
  if (obs) {
    n <- -1 * N(resp, parms, theta)
  } else {
    n <- -1 * exp_N(resp, parms, theta)
  }
  dp <- dIRF(parms, theta)
  d2t <- apply(n * dp * dp, 1, sum, na.rm = T)
  if (obs) {
    m <- M(resp, parms, theta)
    d2p <- d2IRF(parms, theta)
    d2t <- d2t + apply(m * d2p, 1, sum , na.rm = T)
  }
  d2t
}


WA <- function(w, parms, theta1, theta2) {
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  W <- w %*% t(rep(1, nrow(parms)))
  W * (p1 + p2) + (1 - 2 * W) * p1 * p2
}

d_WA <- function(w, parms, theta1, theta2) {
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


d2_WA <- function(w, parms, theta1, theta2) {
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

M <- function(resp, parms, theta) {
  p <- IRF(parms, theta)
  resp / p - (1 - resp) / (1 - p)
}

N <- function(resp, parms, theta) {
  p <- IRF(parms, theta)
  resp / p^2 + (1 - resp) / (1 - p)^2
}

exp_N <- function(resp, parms, theta) {
  p <- IRF(parms, theta)
  1 / p / (1- p) # * !is.na(resp)
}

Mstar <- function(resp, w, parms, theta1, theta2) {
  p <- WA(w, parms, theta1, theta2)
  resp / p - (1 - resp) / (1 - p)
}

Nstar <- function(resp, w, parms, theta1, theta2) {
  p <- WA(w, parms, theta1, theta2)
  resp / p^2 + (1 - resp) / (1 - p)^2
}

exp_Nstar <- function(resp, w, parms, theta1, theta2) {
  p <- WA(w, parms, theta1, theta2)
  1 / p / (1 - p) # * !is.na(resp)
}

l_WA <- function(resp, w, parms, theta1, theta2) {
  p <- WA(w, parms, theta1, theta2)
  apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
}

lpost_WA <- function(resp, w, parms, theta1, theta2, epsilon = .05) {
  lp <- epsilon * log(w - w^2)
  l_WA(resp, w, parms, theta1, theta2) + lp
}

dl_WA <- function(resp, w, parms, theta1, theta2) {
  m <- Mstar(resp, w, parms, theta1, theta2)
  d <- d_WA(w, parms, theta1, theta2)
  dw <- apply(m * d$dw, 1, sum, na.rm = T)
  dt1 <- apply(m * d$dtheta1, 1, sum, na.rm = T)
  dt2 <- apply(m * d$dtheta2, 1, sum, na.rm = T)
  rbind(dw, dt1, dt2) %>% c
}

dlpost_WA <- function(resp, w, parms, theta1, theta2, epsilon = .05) {
  dlw <- epsilon * (1 - 2 * w) / (w - w^2)
  dlp <- rbind(dlw, -theta1, -theta2) %>% c
  dl_WA(resp, w, parms, theta1, theta2) + dlp
}


d2l_WA <- function(resp, w, parms, theta1, theta2, obs = T, parallel = F) {
  if (obs) {
    n <- -1 * Nstar(resp, w, parms, theta1, theta2)
  } else {
    n <- -1 * exp_Nstar(resp, w, parms, theta1, theta2)
  }
  d <- d_WA(w, parms, theta1, theta2)
  dwdw <-  apply(n * d$dw * d$dw, 1, sum, na.rm = T)
  dwdt1 <- apply(n * d$dw * d$dtheta1, 1, sum, na.rm = T)
  dwdt2 <- apply(n * d$dw * d$dtheta2, 1, sum, na.rm = T)
  dt1dt1 <- apply(n * d$dtheta1 * d$dtheta1, 1, sum, na.rm = T)
  dt1dt2 <- apply(n * d$dtheta1 * d$dtheta2, 1, sum, na.rm = T)
  dt2dt2 <- apply(n * d$dtheta2 * d$dtheta2, 1, sum, na.rm = T)

  if (obs) {
    m <- Mstar(resp, w, parms, theta1, theta2)
    d2 <- d2_WA(w, parms, theta1, theta2)
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

  if (parallel) {
    temp <- parallel::mclapply(1:length(w), fun)
  } else {
    temp <- vector("list", length(theta1))
    for (i in 1:length(theta1)) {temp[[i]] <- fun(i)}
  }
  bdiag_m(temp)
}


d2lpost_WA <- function(resp, w, parms, theta1, theta2, epsilon = .05, obs = T, parallel = T) {
  w2 <- w - w^2
  d2lw <- -1 * epsilon * (2 / w2 + ((1 - 2 * w) / w2)^2 )
  d2lp <- rbind(d2lw, -1, -1) %>% c %>% diag
  d2l_WA(resp, w, parms, theta1, theta2, obs = obs, parallel = parallel) + d2lp
}


l_full <- function(resp, w, parms, theta1, theta2)
{
  ind <- grep("IND", names(resp))
  col <- grep("COL", names(resp))
  ind_parms <- parms[grep("IND", row.names(parms)),]
  col_parms <- parms[grep("COL", row.names(parms)),]
  odd <- seq(1, nrow(resp), by = 2)
  (logL(resp[odd, ind], ind_parms, theta1) +
    logL(resp[(odd + 1), ind], ind_parms, theta2) +
    l_WA(resp[odd, col], w, col_parms, theta1, theta2)) %>% sum
}


lpost_full <- function(resp, w, parms, theta1, theta2, epsilon = .05)
{
  l_p <- epsilon * log(w - w^2) - theta1^2 / 2 - theta2^2 / 2
  (l_full(resp, w, parms, theta1, theta2) + l_p) %>% sum
}
d2lpost_full(resp, w, parms, theta1, theta2)

dl_full <- function(resp, w, parms, theta1, theta2) {
  ind <- grep("IND", names(resp))
  col <- grep("COL", names(resp))
  ind_parms <- parms[grep("IND", row.names(parms)),]
  col_parms <- parms[grep("COL", row.names(parms)),]
  odd <- seq(1, nrow(resp), by = 2)

  dl_t1 <- dlogL(resp[odd, ind], ind_parms, theta1)
  dl_t2 <- dlogL(resp[(odd+1), ind], ind_parms, theta2)
  temp <- rbind(0, dl_t1, dl_t2) %>% c
  dl_WA(resp[odd, col], w, col_parms, theta1, theta2) + temp
}


dlpost_full <- function(resp, w, parms, theta1, theta2, epsilon = .05) {
  dlw <- epsilon * (1 - 2 * w) / (w - w^2)
  dlp <- rbind(dlw, -theta1, -theta2) %>% c
  dl_full(resp, w, parms, theta1, theta2) + dlp
}


d2l_full <- function(resp, w, parms, theta1, theta2, obs = T, parallel = T) {
  ind <- grep("IND", names(resp))
  col <- grep("COL", names(resp))
  ind_parms <- parms[grep("IND", row.names(parms)),]
  col_parms <- parms[grep("COL", row.names(parms)),]
  odd <- seq(1, nrow(resp), by = 2)

  d2l_t1 <- d2logL(resp[odd, ind], ind_parms, theta1, obs)
  d2l_t2 <- d2logL(resp[(odd+1), ind], ind_parms, theta2, obs)
  temp <- rbind(0, d2l_t1, d2l_t2) %>% c %>% diag
  d2l_WA(resp[odd, col], w, col_parms, theta1, theta2, obs, parallel) + temp
}


d2lpost_full <- function(resp, w, parms, theta1, theta2, obs = T, parallel = T, epsilon = .05) {
  w2 <- w - w^2
  d2lw <- -1 * epsilon * (2 / w2 + ((1 - 2 * w) / w2)^2 )
  d2lp <- rbind(d2lw, -1, -1) %>% c %>% diag
  d2l_full(resp, w, parms, theta1, theta2, obs, parallel) + d2lp
}


est_WA <- function(resp, parms, starts = NULL, SE = "obs", method = "map", epsilon = .05, parallel = T) {
  n_obs <- nrow(resp)/2
  odd <- seq(from = 1, to = n_obs*2, by = 2)
  parm_index <- seq(from = 1, to = n_obs*3, by = 3)
  out <- data.frame(matrix(0, nrow = n_obs, ncol = 7))
  names(out) <- c("log", "w", "w_se", "theta1", "theta1_se", "theta2", "theta2_se")

  if(is.null(starts)) {starts <- rep(c(.5, 0, 0), times = n_obs)}
  lower <- rep(c(.00001, -8, -8), times = n_obs)
  upper <- rep(c(.99999, 8, 8), times = n_obs)

  if (method == "map"){
    obj <- function(par, resp) {
      ind <- parm_index[1:(length(par) / 3)]
      -1 * lpost_full(resp, par[ind], parms, par[(ind+1)], par[(ind+2)], epsilon)
    }
    grad <- function(par, resp) {
      ind <- parm_index[1:(length(par) / 3)]
      -1 * dlpost_full(resp, par[ind], parms, par[(ind+1)], par[(ind+2)], epsilon)
    }
  }

  if (method == "ml"){
    obj <- function(par, resp) {
      ind <- parm_index[1:(length(par) / 3)]
      -1 * l_full(resp, par[ind], parms, par[(ind+1)], par[(ind+2)])
    }
    grad <- function(par, resp) {
      ind <- parm_index[1:(length(par) / 3)]
      -1 * dl_full(resp, par[ind], parms, par[(ind+1)], par[(ind+2)])
    }
  }

  # set up parm_indexing for parallel
  fun <- function(i) {
    blocksize <- floor(n_obs/n_cores)
    m <- (i-1) * blocksize + 1
    if (i < n_cores) {n <- i * blocksize} else {n <- n_obs}
    ind1 <- index[m] : (index[n] + 2)
    ind2 <- odd[m] : (odd[n] + 1)
    q <- optim(starts[ind1], obj,
             gr = grad,
             resp = resp[ind2,],
             method = "L-BFGS-B",
             lower = lower[ind1],
             upper = upper[ind1]
             )$par
  }

  if(parallel) {
    n_cores <- parallel::detectCores()
    temp <- parallel::mclapply(1:n_cores, fun) %>% unlist
  } else {
    n_cores <- 1
    temp <- fun(1)
  }
  out$w <- temp[parm_index]
  out$theta1 <- temp[parm_index+1]
  out$theta2 <- temp[parm_index+2]

  if (SE == "obs") {obs = T} else {obs = F}
  if (method == "map"){
    temp_se <- (-1 * d2lpost_full(resp, out$w, parms, out$theta1, out$theta2, obs, parallel)) %>% solve %>% diag %>% sqrt
  }
  if (method == "ml"){
    temp_se <- (-1 * d2l_full(resp, out$w, parms, out$theta1, out$theta2, obs, parallel)) %>% solve %>% diag %>% sqrt
  }

  out$w_se <- temp_se[parm_index]
  out$theta1_se <- temp_se[parm_index+1]
  out$theta2_se <- temp_se[parm_index+2]
  out$log <- -1 * obj(temp, resp) # un-sum
  out
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


#-------------------------------------------------------------------
#Esimation for just WA
#-------------------------------------------------------------------

m_WA <- function(resp, w, parms, theta1, theta2, Log = T, Sum = F) {
  alpha <- 1.05
  p <- WA(w, parms, theta1, theta2)
  temp <- apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
  temp <- temp + (alpha - 1) * log(w - w^2)
  if (Sum) {temp <- sum(temp)}
  if (Log) {temp} else {exp(temp)}
}


pgrad_WA <- function(resp, w, parms, theta1, theta2, Sum = T) {
  alpha <- 1.05
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  dw <- p1 + p2 - 2 * p1 * p2
  m <- Mstar(resp, w, parms, theta1, theta2)
  lp <- (alpha - 1) / (w - w^2) * (1 - 2 * w)
  apply(m * dw, 1, sum, na.rm = T) + lp
}

# pgrad_WA <- function(resp, w, parms, theta1, theta2, Sum = T) {
#   alpha <- 1.05
#   r <- WA(w, parms, theta1, theta2)
#   p1 <- IRF(parms, theta1)
#   p2 <- IRF(parms, theta2)
#   m <- (alpha - 1) / (w - w^2) * (1 - 2 * w)
#   n <- resp/r - (1 - resp)/(1 - r)
#   temp <- n * (p1 * (1 - p2) + p2 * (1 - p1)) + m
#   if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
# }


pinfo_WA <- function(resp, w, parms, theta1, theta2, type = "obs", Sum = T) {
  alpha <- 1.05
  r <- WA(w, parms, theta1, theta2)
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  m <- (1 - alpha) / (w - w^2) * (2 + (1 - 2 * w)^2 / (w - w ^2))
  if (type == "obs") {
    n <- resp/r^2 + (1 - resp)/(1 - r)^2
  } else {
    n <- 1/r/(1 - r) #*!is.na(resp)
  }
  temp <- n * (p1*(1 - p2) + p2*(1 - p1))^2 - m
  if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
}

map_WA <- function(resp, parms, theta1, theta2, SE = "obs", starts = NULL, parallel = F) {
  n_obs <- nrow(resp)
  out <- data.frame(matrix(0, nrow = n_obs, ncol = 3))
  names(out) <- c("logp", "w", "psd")
  if(is.null(starts)) starts <- rep(.5, n_obs)

  map <- function(w, resp, theta1, theta2) {
    -1*m_WA(resp, w, parms, theta1, theta2, Sum = T)
  }

  dmap <- function(w, resp, theta1, theta2) {
    -1*pgrad_WA(resp, w, parms, theta1, theta2)
  }

  fun <- function(i) {
    blocksize <- floor(n_obs/n_cores)
    m <- (i-1)*blocksize + 1
    if (i < n_cores) {n <- i*blocksize} else {n <- n_obs}
    optim(starts[m:n], map,
             gr = dmap,
             resp = resp[m:n,],
             theta1 = theta1[m:n],
             theta2 = theta2[m:n],
             method = "L-BFGS-B",
             lower = rep(.000001, n-m+1),
             upper = rep(.999999, n-m+1)
             )$par
  }
  if(parallel){
    n_cores <- parallel::detectCores()
    out$w<- parallel::mclapply(1:n_cores, fun) %>% unlist
  } else {
    n_cores <- 1
    out$w <- fun(1)
  }
  out$psd <- 1/sqrt(pinfo_WA(resp, out$w, parms, theta1, theta2, type = SE))
  out$logp <- m_WA(resp, out$w, parms, theta1, theta2)
  out
}

# # ML -------------------------------------------------------------------------
#
# l_WA <- function(resp, w, parms, theta1, theta2, Log = T, Sum = F) {
#   p <- WA(w, parms, theta1, theta2)
#   temp <- apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
#   if (Sum) {temp <- sum(temp)}
#   if (Log) {temp} else {exp(temp)}
# }
#
# grad_WA <- function(resp, w, parms, theta1, theta2, Sum = T) {
#   r <- WA(w, parms, theta1, theta2)
#   p1 <- IRF(parms, theta1)
#   p2 <- IRF(parms, theta2)
#   n <- resp/r - (1 - resp)/(1 - r)
#   temp <- n * (p1*(1 - p2) + p2*(1 - p1))
#   #temp <- n * (p1 + p2 - p1 * p2)
#   if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
# }
#
#
# info_WA <- function(resp, w, parms, theta1, theta2, type = "obs", Sum = T) {
#   r <- WA(w, parms, theta1, theta2)
#   p1 <- IRF(parms, theta1)
#   p2 <- IRF(parms, theta2)
#   if (type == "obs") {
#     n <- resp/r^2 + (1 - resp)/(1 - r)^2
#   } else {
#     n <- 1/r/(1 - r)*!is.na(resp)
#   }
#   temp <- n * (p1*(1 - p2) + p2*(1 - p1))^2
#   if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
# }
#
# mle_WA <- function(resp, parms, theta1, theta2, SE = "obs", starts = NULL, parallel = F) {
#   n_obs <- nrow(resp)
#   out <- data.frame(matrix(0, nrow = n_obs, ncol = 3))
#   names(out) <- c("logL", "w", "se")
#   if(is.null(starts)) starts <- rep(.5, n_obs)
#
#   logl <- function(w, resp, theta1, theta2) {
#     -1*l_WA(resp, w, parms, theta1, theta2, Sum = T)
#   }
#
#   dlogl <- function(w, resp, theta1, theta2) {
#     -1*grad_WA(resp, w, parms, theta1, theta2)
#   }
#
#   fun <- function(i) {
#     blocksize <- floor(n_obs/n_cores)
#     m <- (i-1)*blocksize + 1
#     if (i < n_cores) {n <- i*blocksize} else {n <- n_obs}
#     optim(starts[m:n], logl,
#              gr = dlogl,
#              resp = resp[m:n,],
#              theta1 = theta1[m:n],
#              theta2 = theta2[m:n],
#              method = "L-BFGS-B",
#              lower = rep(.000001, n-m+1),
#              upper = rep(.999999, n-m+1)
#              )$par
#   }
#   if(parallel){
#     n_cores <- parallel::detectCores()
#     out$w<- parallel::mclapply(1:n_cores, fun) %>% unlist
#   } else {
#     n_cores <- 1
#     out$w <- fun(1)
#   }
#   out$se <- 1/sqrt(info_WA(resp, out$w, parms, theta1, theta2, type = SE))
#   out$logL <- l_WA(resp, out$w, parms, theta1, theta2)
#   out
# }
#
#
#
#
# # l_full <- function(resp, parms, theta1, theta2, w)
# # {
# #   ind <- grep("IND", names(resp))
# #   col <- grep("COL", names(resp))
# #   -1*(logL(resp[1, ind], parms, theta1) + logL(resp[2, ind], parms, theta2) + l_WA(resp[1, col], w, parms, theta1, theta2) )
# # }
#
# WA2 <- function(w, parms, theta1, theta2, theta1_se, theta2_se) {
#   p1_lower <- IRF(parms, theta1 - 1.96*theta1_se)
#   p1_upper <- IRF(parms, theta1 + 1.96*theta1_se)
#   p2_lower <- IRF(parms, theta2 - 1.96*theta2_se)
#   p2_upper <- IRF(parms, theta2 + 1.96*theta2_se)
#   W <- w %*% t(rep(1, nrow(parms)))
#   temp <- p1_lower * p2_lower + W * (p1_upper + p2_upper - 2 * p1_upper * p2_upper)
#   #temp[temp >= .99999] <- .99999
#   #temp
# }
#
#
# l_WA2 <- function(resp, w, parms, theta1, theta2, theta1_se, theta2_se, Log = T, Sum = F) {
#   r <- WA2(w, parms, theta1, theta2, theta1_se, theta2_se)
#   temp <- apply(log(r) * (resp) + log(1-r) * (1-resp), 1, sum, na.rm = T)
#   if (Sum) {temp <- sum(temp)}
#   if (Log) {temp} else {exp(temp)}
# }
#
# grad_WA2 <- function(resp, w, parms, theta1, theta2, theta1_se, theta2_se, Sum = T) {
#   r <- WA2(w, parms, theta1, theta2, theta1_se, theta2_se)
#   p1_lower <- IRF(parms, theta1 - 1.96*theta1_se)
#   p1_upper <- IRF(parms, theta1 + 1.96*theta1_se)
#   p2_lower <- IRF(parms, theta2 - 1.96*theta2_se)
#   p2_upper <- IRF(parms, theta2 + 1.96*theta2_se)
#   n <- resp/r - (1 - resp)/(1 - r)
#   temp <- n * (p1_upper * (1 - p2_lower) + (1 - p1_lower) * p2_upper)
#   if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
# }
#
# info_WA2 <- function(resp, w, parms, theta1, theta2, theta1_se, theta2_se, type = "obs", Sum = T) {
#   r <- WA2(w, parms, theta1, theta2, theta1_se, theta2_se)
#   p1_lower <- IRF(parms, theta1 - 1.96*theta1_se)
#   p1_upper <- IRF(parms, theta1 + 1.96*theta1_se)
#   p2_lower <- IRF(parms, theta2 - 1.96*theta2_se)
#   p2_upper <- IRF(parms, theta2 + 1.96*theta2_se)
#   if (type == "obs") {
#     n <- resp/r^2 + (1 - resp)/(1 - r)^2
#   } else {
#     n <- 1/r/(1 - r)*!is.na(resp)
#   }
#   temp <- n * (p1_upper * (1 - p2_lower) + (1 - p1_lower) * p2_upper)^2
#   if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
# }
#
#
# mle_WA2 <- function(resp, parms, theta1, theta2, theta1_se, theta2_se, SE = "obs", starts = NULL, parallel = F) {
#   n_obs <- nrow(resp)
#   out <- data.frame(matrix(0, nrow = n_obs, ncol = 3))
#   names(out) <- c("logL", "w", "se")
#   if(is.null(starts)) starts <- rep(.5, n_obs)
#
#   logl <- function(w, resp, theta1, theta2, theta1_se, theta2_se) {
#     -1*l_WA2(resp, w, parms, theta1, theta2, theta1_se, theta2_se, Sum = T)
#   }
#
#   dlogl <- function(w, resp, theta1, theta2, theta1_se, theta2_se) {
#     -1*grad_WA2(resp, w, parms, theta1, theta2, theta1_se, theta2_se)
#   }
#
#   fun <- function(i) {
#     blocksize <- floor(n_obs/n_cores)
#     m <- (i-1)*blocksize + 1
#     if (i < n_cores) {n <- i*blocksize} else {n <- n_obs}
#     optim(starts[m:n], logl,
#              gr = dlogl,
#              resp = resp[m:n,],
#              theta1 = theta1[m:n],
#              theta2 = theta2[m:n],
#              theta1_se = theta1_se[m:n],
#              theta2_se = theta2_se[m:n],
#              method = "L-BFGS-B",
#              lower = rep(0, n-m+1),
#              upper = rep(1, n-m+1)
#              )$par
#   }
#   if(parallel){
#     n_cores <- parallel::detectCores()
#     out$w<- parallel::mclapply(1:n_cores, fun) %>% unlist
#   } else {
#     n_cores <- 1
#     out$w <- fun(1)
#   }
#   out$se <- 1/sqrt(info_WA2(resp, out$w, parms, theta1, theta2, theta1_se, theta2_se, type = SE))
#   out$logL <- l_WA2(resp, out$w, parms, theta1, theta2, theta1_se, theta2_se)
#   out
# }
#
#
#

#--------------------------------------------------------------------------
#' Simulate item responses from the 2PL or a model of pairwise collaboration obtained from the 2PL. Called by \code{data_gen}.

#' @param model is one of \code{c("IRF", "Ind", "Min", "Max", "AI") }.
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait for member 1.
#' @param theta2 the latent trait for member 2.
#' @return \code{length(theta1)} by \code{nrow(parms)} matrix of binary responses.
#' @export

sim_WA <- function(w, parms, theta1 = 0, theta2 = 0) {
  n_row <- length(theta1)
  n_col <- nrow(parms)
  r <- array(runif(n_row * n_col), dim = c(n_row, n_col))
  p <- WA(w, parms, theta1, theta2)
  out <- ifelse(p > r, 1, 0)
  colnames(out) <- row.names(parms)
  out
}

#--------------------------------------------------------------------------
#' Simulate data from a mixture of models of collaboration.
#'
#' @param n_reps integer indicating how many datasets to generate.
#' @param mix_prop a n_models-vector by \code{length(theta)} matrix of mixing proportions for the models.
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param theta1_se the standard error of the latent trait for member 1 (Optional). If included, data generation incorporates PV for theta1.
#' @param theta2_se the standard error of the latent trait for member 2 (Optional). If included, data generation incorporates PV for theta2.
#' @param NA_pattern an (optional) \code{length(theta)} by \code{nrow(parms)} data.frame with \code{NA} entries for each item a dyad did not answer. The missing values are preserved in the generated data.

#' @return A data.frame with \code{length(theta)} rows containing an id variable for each pair and each sample, the data generating values of theta1, theta2, and mix_prop; the model used to simulate the response pattern; and the simulated response pattern.
#' @export

data_gen <- function(n_reps, w, parms, theta1, theta2, theta1_se = NULL, theta2_se = NULL, NA_pattern = NULL) {

  # Storage
  n_obs <- length(theta1)
  out <- data.frame(rep(1:n_obs, each = n_reps), rep(1:n_reps, times = n_obs))
  names(out) <- c("pairs", "samples")
  out$w <- rep(w, each = n_reps)
  # Generate / replicate thetas
  #temp_theta <- theta_sort(theta1, theta2, theta1_se, theta2_se)
  #out$theta1 <- theta_gen(n_reps, temp_theta$theta_min, temp_theta$se_min)
  #out$theta2 <- theta_gen(n_reps, temp_theta$theta_max, temp_theta$se_max)
  out$theta1 <- theta_gen(n_reps, theta1, theta1_se)
  out$theta2 <- theta_gen(n_reps, theta2, theta2_se)

  # Simulate data
  data <- data.frame(matrix(NA, nrow = n_reps*n_obs, ncol = nrow(parms)))
  names(data) <- row.names(parms)
  data <- sim_WA(out$w, parms, out$theta1, out$theta2)
  data <- format_NA(data, NA_pattern)

  # Return
  cbind(out[], data[])
}


# #--------------------------------------------------------------------------
# #' Generate plausible values.
# #'
# #  Replicates the input data and combines it with draws from the "approximate" (i.e., normal) distribution of theta1 and theta2
# #' @param n_reps integer indicating how many replications to use.
# #' @param resp data.frame of binary, conjuncitvely-scored responses patterns.
# #' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
# #' @param theta1 the latent trait for member 1.
# #' @param theta2 the latent trait for member 2.
# #' @param theta1_se the standard error of the latent trait for member 1.
# #' @param theta2_se the standard error of the latent trait for member 2.
# #' @param true_model an (optional) \code{length(theta)} vector indicating which model is associated with each value of theta. Useful when working with simulated data.
#
# #' @return A data.frame with \code{length(theta) \times n_reps} rows containing an id variable for each pair and for each sample, the plausible values of theta1, theta2, the observed data for each plausible value, and (optionally) the true_model for each value plausible value.
# #' @export
#
# pv_gen <- function(n_reps, resp, parms, theta1, theta2, theta1_se, theta2_se, weights = NULL) {
#
#   # Expand data generating parms
#   n_obs <- length(theta1)
#   n_long <- n_obs * n_reps
#   out <- data.frame(rep(1:n_obs, each = n_reps), rep(1:n_reps, times = n_obs))
#   names(out) <- c("pairs", "samples")
#   if (!is.null(weights)) {out$weights <- rep(weights, each = n_reps) }
#
#   # Sort thetas
#   temp_theta <- theta_sort(theta1, theta2, theta1_se, theta2_se)
#
#   # Generate PVs for theta
#   out$theta1 <- theta_gen(n_reps, temp_theta$theta_min, temp_theta$se_min)
#   out$theta2 <- theta_gen(n_reps, temp_theta$theta_max, temp_theta$se_max)
#
#   # Expand data (replicate each row n_reps times)
#   data <- kronecker(as.matrix(resp), rep(1, n_reps))
#   colnames(data) <- row.names(parms)
#   cbind(out[], data[])
# }
