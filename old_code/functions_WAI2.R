# Last update: 04/12/2016
# Functions for computing bootstrapped LRtests of collaboration models.
# User beware: functions not written to check or handle input errors.

WA <- function(w, parms, theta1, theta2) {
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  W <- w %*% t(rep(1, nrow(parms)))
  W * (p1 + p2) + (1 - 2 * W) * p1 * p2
}

l_WA <- function(resp, w, parms, theta1, theta2, Log = T, Sum = F) {
  p <- WA(w, parms, theta1, theta2)
  temp <- apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
  if (Sum) {temp <- sum(temp)}
  if (Log) {temp} else {exp(temp)}
}

grad_WA <- function(resp, w, parms, theta1, theta2, Sum = T) {
  r <- WA(w, parms, theta1, theta2)
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  n <- resp/r - (1 - resp)/(1 - r)
  temp <- n * (p1*(1 - p2) + p2*(1 - p1))
  if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
}

info_WA <- function(resp, w, parms, theta1, theta2, type = "obs", Sum = T) {
  r <- WA(w, parms, theta1, theta2)
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  if (type == "obs") {
    n <- resp/r^2 + (1 - resp)/(1 - r)^2
  } else {
    n <- 1/r/(1 - r)*!is.na(resp)
  }
  temp <- n * (p1*(1 - p2) + p2*(1 - p1))^2
  if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
}

mle_WA <- function(resp, parms, theta1, theta2, SE = "obs", starts = NULL, parallel = F) {
  n_obs <- nrow(resp)
  out <- data.frame(matrix(0, nrow = n_obs, ncol = 3))
  names(out) <- c("logL", "w", "se")
  if(is.null(starts)) starts <- rep(.5, n_obs)

  logl <- function(w, resp, theta1, theta2) {
    -1*l_WA(resp, w, parms, theta1, theta2, Sum = T)
  }

  dlogl <- function(w, resp, theta1, theta2) {
    -1*grad_WA(resp, w, parms, theta1, theta2)
  }

  fun <- function(i) {
    blocksize <- floor(n_obs/n_cores)
    m <- (i-1)*blocksize + 1
    if (i < n_cores) {n <- i*blocksize} else {n <- n_obs}
    optim(starts[m:n], logl,
             gr = dlogl,
             resp = resp[m:n,],
             theta1 = theta1[m:n],
             theta2 = theta2[m:n],
             method = "L-BFGS-B",
             lower = rep(0, n-m+1),
             upper = rep(1, n-m+1)
             )$par
  }
  if(parallel){
    n_cores <- parallel::detectCores()
    out$w<- parallel::mclapply(1:n_cores, fun) %>% unlist
  } else {
    n_cores <- 1
    out$w <- fun(1)
  }
  out$se <- 1/sqrt(info_WA(resp, out$w, parms, theta1, theta2, type = SE))
  out$logL <- l_WA(resp, out$w, parms, theta1, theta2)
  out
}


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
  temp_theta <- theta_sort(theta1, theta2, theta1_se, theta2_se)
  out$theta1 <- theta_gen(n_reps, temp_theta$theta_min, temp_theta$se_min)
  out$theta2 <- theta_gen(n_reps, temp_theta$theta_max, temp_theta$se_max)

  # Simulate data
  data <- data.frame(matrix(NA, nrow = n_reps*n_obs, ncol = nrow(parms)))
  names(data) <- row.names(parms)
  data <- sim_WA(out$w, parms, out$theta1, out$theta2)
  data <- format_NA(data, NA_pattern)

  # Return
  cbind(out[], data[])
}


#--------------------------------------------------------------------------
#' Generate plausible values.
#'
#  Replicates the input data and combines it with draws from the "approximate" (i.e., normal) distribution of theta1 and theta2
#' @param n_reps integer indicating how many replications to use.
#' @param resp data.frame of binary, conjuncitvely-scored responses patterns.
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait for member 1.
#' @param theta2 the latent trait for member 2.
#' @param theta1_se the standard error of the latent trait for member 1.
#' @param theta2_se the standard error of the latent trait for member 2.
#' @param true_model an (optional) \code{length(theta)} vector indicating which model is associated with each value of theta. Useful when working with simulated data.

#' @return A data.frame with \code{length(theta) \times n_reps} rows containing an id variable for each pair and for each sample, the plausible values of theta1, theta2, the observed data for each plausible value, and (optionally) the true_model for each value plausible value.
#' @export

pv_gen <- function(n_reps, resp, parms, theta1, theta2, theta1_se, theta2_se, weights = NULL) {

  # Expand data generating parms
  n_obs <- length(theta1)
  n_long <- n_obs * n_reps
  out <- data.frame(rep(1:n_obs, each = n_reps), rep(1:n_reps, times = n_obs))
  names(out) <- c("pairs", "samples")
  if (!is.null(weights)) {out$weights <- rep(weights, each = n_reps) }

  # Sort thetas
  temp_theta <- theta_sort(theta1, theta2, theta1_se, theta2_se)

  # Generate PVs for theta
  out$theta1 <- theta_gen(n_reps, temp_theta$theta_min, temp_theta$se_min)
  out$theta2 <- theta_gen(n_reps, temp_theta$theta_max, temp_theta$se_max)

  # Expand data (replicate each row n_reps times)
  data <- kronecker(as.matrix(resp), rep(1, n_reps))
  colnames(data) <- row.names(parms)
  cbind(out[], data[])
}
