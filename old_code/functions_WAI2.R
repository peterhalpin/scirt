# Last update: 04/12/2016
# Functions for computing bootstrapped LRtests of collaboration models.
# User beware: functions not written to check or handle input errors.

WA <- function(w, parms, theta1, theta2) {
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  W <- w %*% t(rep(1, nrow(parms)))
  W * (p1 + p2) + (1 - 2 * W) * p1 * p2
}

l_WA <- function(resp, w, parms, theta1, theta2, Log = T) {
  p <- WA(w, parms, theta1, theta2)
  out <- apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
  if (Log) {out} else {exp(out)}
}

info_WA <- function(resp, w, parms, theta1, theta2, type = "obs", Sum = T) {
  r <- WA(w, parms, theta1, theta2)
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  if (type == "obs") {
    n <- resp/r^2 + (1 - resp)/(1 - r)^2
  } else {
    n <- 1/r/(1 - r)
  }
  temp <- n * (p1*(1 - p2) + p2*(1 - p1))^2
  if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
}

ML_WA <- function(resp, parms, theta1, theta2, SE = "obs", parallel = T) {
  n_obs <- nrow(resp)
  out <- data.frame(matrix(0, nrow = n_obs, ncol = 3))
  names(out) <- c("logL", "w", "se")

  fun <- function(w, resp, theta1, theta2)  {
    -l_WA(resp, w, parms, theta1, theta2)
  }

  ml <- function(i) {
    optimize(fun, c(0, 1), resp = resp[i,], theta1 = theta1[i], theta2 = theta2[i])
  }
w
  if (parallel) {
    q <-  parallel::mclapply(1:n_obs, ml) %>% unlist
    ind <- seq(1, n_obs*2, by = 2)
    out$w <- q[ind]
    out$logL <- -1*q[ind+1]
  } else {
    for (i in 1:nrow(resp)) {
      temp <- ml(i)
      out$w[i] <- temp$minimum
      out$logL <- -1*temp$objective
    }
  }
  out$se <- 1/sqrt(info_WA(resp, out$w, parms, theta1, theta2))
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
