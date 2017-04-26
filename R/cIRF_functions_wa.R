# Last update: 04/25/2017
# Functions for testing collaboration models. Depends on IRF_functions.R
# User beware: functions do not check or handle input errors.

require(ggplot2)
require(dplyr)


#--------------------------------------------------------------------------
#' Item response function for ``the Independence model"
#'
#' Computes a matrix of probabilities for correct responses using the independence model.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param sorted logical used for compatbility with \code{cIRF}
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities.
#' @export

Ind <- function(parms, theta1, theta2, sorted = F) {
  IRF(parms, theta1)*IRF(parms, theta2)
}


#--------------------------------------------------------------------------
#' Item response function for ``the Minimum model"
#'
#' Computes a matrix of probabilities for correct responses using the minimum model.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param sorted logical indicating whether to compute Min using \code{theta1}, regardless of the value of the \code{theta2} (useful for data simulation, where the simulated values of theta may not reflect the ordering of the data generating values)
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities.
#' @export

Min <- function(parms, theta1, theta2, sorted = F) {
  if (sorted) {
    IRF(parms, theta1)
  } else {
    theta <- apply(cbind(theta1, theta2), 1, min, na.rm = T)
    IRF(parms, theta)
  }
}


#--------------------------------------------------------------------------
#' Item response function for ``the Maximum model"
#'
#' Computes a matrix of probabilities for correct responses using the maximum model.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param sorted logical indicating whether to compute Min using \code{theta1}, regardless of the value of the \code{theta2} (useful for data simulation, where the simulated values of theta may not reflect the ordering of the data generating values)
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

Max <- function(parms, theta1, theta2, sorted = F) {
  if (sorted) {
    IRF(parms, theta2)
  } else {
    theta <- apply(cbind(theta1, theta2), 1, max, na.rm = T)
    IRF(parms, theta)
  }
}


#--------------------------------------------------------------------------
#' Item response function for ``the Additive model"
#'
#' Computes a matrix of probabilities for correct responses using the additive model.
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param sorted logical used for compatbility with \code{cIRF}
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

Add <- function(parms, theta1, theta2, sorted = F) {
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  p1 + p2 - p1 * p2
}

#--------------------------------------------------------------------------
#' Wrapper for collaborative item response functions
#'
#' Computes a matrix of probabilities for correct responses using the named \code{model} for collaboration, and the 2PL model for the items.
#'
#' @param model one of \code{c("Ind", "Min", "Max", "Add", "IRF")}
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param sorted logical indicating whether to compute Min using \code{theta1}, regardless of the value of the \code{theta2} (useful for data simulation, where the simulated values of theta may not reflect the ordering of the data generating values)
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

cIRF <- function(model, parms, theta1, theta2, sorted = F) {
  if (model %in% c("Ind", "Min", "Max", "Add", "IRF")) {
    fun <- match.fun(model)
    fun(parms, theta1, theta2, sorted)
  } else {
    cat("\'model\' must be one of c(\"Ind\", \"Min\", \"Max\", \"Add\", \"IRF\")")
  }
}

#--------------------------------------------------------------------------
#' Likelihood of a matrix of binary responses for one or more models of collaboration, conditional on theta.
#'
#' \code{logL} is faster for 2PL.
#'
#' @param models one or more of \code{c("IRF", "Ind", "Min", "Max", "Add") }
#' @param resp a matrix or data.frame containing the binary item responses
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param sorted logical indicating whether to compute Min using \code{theta1}, regardless of the value of the \code{theta2} (useful for data simulation, where the simulated values of theta may not reflect the ordering of the data generating values)
#' @param Log logical indicating whether to return loglikelihood or likelihood
#' @return An \code{nrow(resp)} by \code{length(models)} matrix of (log-) likleihoods for each response pattern and each model
#' @export

likelihood <- function(models, resp, parms, theta1, theta2 = NULL, sorted = F, Log = T) {
  n_models <- length(models)
  out <- array(0, dim = c(nrow(resp), n_models))
  for (i in 1:n_models) {
    p <- cIRF(models[i], parms, theta1, theta2, sorted)
    out[,i] <- apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
  }
  if (n_models == 1) {out <- c(out)} # un-matrix
  if (Log) {out} else {exp(out)}
}

#--------------------------------------------------------------------------
#' Item response function for ``the Weighted Additive model"
#'
#' Computes a matrix of probabilities for correct responses using the weighted additive model.
#'
#' @param w a numeric vector containing the weights for each group
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

WA <- function(w, parms, theta1, theta2) {
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  W <- w %*% t(rep(1, nrow(parms)))
  W * (p1 + p2) + (1 - 2 * W) * p1 * p2
}

#--------------------------------------------------------------------------
#' (Log-) likelihood for ``the Weighted Additive model"
#'
#' @param resp a matrix or data.frame containing the binary item responses
#' @param w a numeric vector containing the weights for each group
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param Log logical: return the loglikelihood or the likelihood?
#' @param Sum logical: return the sum over items, for each group, or a matrix for each group and each item?
#'
#' @return If \code{Sum = T), a \code{length(theta)} vector of (log-) likelihoods for each group. If \code{Sum = T), a \code{length(theta)} by \code{nrow(parms)} matrix of (log-) likelihoods for each group and each item
#' @export

l_WA <- function(resp, w, parms, theta1, theta2, Log = T, Sum = F) {
  p <- WA(w, parms, theta1, theta2)
  temp <- apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
  if (Sum) {temp <- sum(temp)}
  if (Log) {temp} else {exp(temp)}
}

#--------------------------------------------------------------------------
#' First derivative (in w) of loglikelihood for ``the Weighted Additive model"
#'
#' @param resp a matrix or data.frame containing the binary item responses
#' @param w a numeric vector containing the weights for each group
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param Sum logical: return the sum over items, for each group, or a matrix for each group and each item?
#'
#' @return If \code{Sum = T), a \code{length(theta)} vector of (log-) likelihoods for each group. If \code{Sum = T), a  \code{length(theta)} by \code{nrow(parms)} matrix of (log-) likelihoods for each group and each item
#' @export

grad_WA <- function(resp, w, parms, theta1, theta2, Sum = T) {
  r <- WA(w, parms, theta1, theta2)
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  n <- resp/r - (1 - resp)/(1 - r)
  temp <- n * (p1*(1 - p2) + p2*(1 - p1))
  if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
}


#--------------------------------------------------------------------------
#' Second derivative (in w) of loglikelihood for ``the Weighted Additive model"
#'
#' @param resp a matrix or data.frame containing the binary item responses
#' @param w a numeric vector containing the weights for each group
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param type one of \code{"obs", "exp"}, indicating whether to return the observed or expected Fisher information
#' @param Sum logical: return the sum over items, for each group, or a matrix for each group and each item?
#'
#' @return If \code{Sum = T), a \code{length(theta)} vector of (log-) likelihoods for each group. If \code{Sum = T), a  \code{length(theta)} by \code{nrow(parms)} matrix of (log-) likelihoods for each group and each item
#' @export

info_WA <- function(resp, w, parms, theta1, theta2, type = "obs", Sum = T) {
  r <- WA(w, parms, theta1, theta2)
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  if (type == "obs") {
    n <- resp/r^2 + (1 - resp)/(1 - r)^2
  } else {
    n <- 1/r/(1 - r)#*!is.na(resp)
  }
  temp <- n * (p1*(1 - p2) + p2*(1 - p1))^2
  if (Sum) {apply(temp, 1, sum, na.rm = T)} else {temp}
}


#--------------------------------------------------------------------------
#' Maximum likelihood estimation of weights for ``the Weighted Additive model"
#'
#' @param resp a matrix or data.frame containing the binary item responses
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param SE one of \code{"obs", "exp"}, indicating whether to return the observed or expected Fisher information
#' @param starts a numertic vector of \code{length(theta)} starting values for the weights (optional)
#' @param parallel logical: call \code{parallel:mclapply} instead of looping over \code{nrow(resp)}?
#'
#' @return An \code{nrow(resp)} by 3 matrix of loglikleihoods, theta estimates, and SEs for each response pattern in resp.
#' @export

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
#' Simulate item responses from ``the Weighted Additive model"

#' Called by \code{data_gen}.

#' @param w a numeric vector containing the weights for each group
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of binary responses
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
#' Generate data from ``the Weighted Additive model"

#' @param n_reps integer indicating how many response vectors to generate for group
#' @param w a numeric vector containing the weights for each group
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param theta1_se the standard error of the latent trait for member 1 (optional -- if included, data generation incorporates plausible values for theta1)
#' @param theta1_se the standard error of the latent trait for member 2 (optional -- if included, data generation incorporates plausible values for theta2)
#' @param NA_pattern a \code{length(theta)} by \code{nrow(parms)} data.frame with \code{NA} entries denoting missing values to be preserved in the generated data (optional -- non-\code{NA} entires are ignored, see \code{format_NA})

#' @return A data.frame with \code{length(theta) * n_reps} rows and columns denoting: an id variable for each group; an id variable for each replication; the data generating values of theta1, theta2, and w; and the simulated response patterns.
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
#' Generate plausible values for group response data
#'
#  Replicates the input data and combines it with draws from the "approximate" (i.e., normal) distribution of \code{theta1} and code{theta2}.

#' @param n_reps integer indicating how many replication to use for group
#' @param resp a matrix or data.frame containing the binary item responses for each group
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param theta1_se the standard error of the latent trait for member 1
#' @param theta2_se the standard error of the latent trait for member 2
#' @param weights a \code{length(theta)} vector indicating weights (from the weighted additive model) associated with each group (optional)

#' @return A data.frame with \code{length(theta) * n_reps} rows and columns denoting: an id variable for each group; an id variable for each replication; the generated values of \code{theta1} and \code{theta2}; the \code{weights} for each group; and the replicated response patterns.
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


#--------------------------------------------------------------------------
#' Computes a martix of item deltas (redundancy measure)
#'
#' @param parms a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param sorted logical indicating whether to compute Min using \code{theta1}, regardless of the value of the \code{theta2} (useful for data simulation, where the simulated values of theta may not reflect the ordering of the data generating values)
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of item deltas.
#' @export

item_delta <- function(parms, theta1, theta2, sorted = F, NA_pattern = NULL) {
  temp <- Min(parms, theta1, theta2, sorted) * (1 - Max(parms, theta1, theta2, sorted))
  format_NA(temp, NA_pattern)
}




# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Helper functions----------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------


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
    resp <- resp[grep(version, names(resp))]
  }
  names(resp) <- substr(names(resp), 1, 5)
  resp <- resp[names(resp)%in%items]
  resp[items[!items%in%names(resp)]] <- NA
  resp <- resp[items]
  resp
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
#' Generates or replicates values of latent variable
#'
#' If \code{theta_se}' is no null, values are generated using \code{rnorm}. Otherwise, each value is replicated.

#' @param n number of values to generate / replicate for each value of \code{theta}.
#' @param theta the latent trait.
#' @param theta_se the standard error of the latent trait.

#' @return An \code{n \times length(theta)} vector in which each value \code{theta} is generated / replicated \code{n} times.
#' @export

theta_gen <- function(n, theta, theta_se = NULL){
  theta_long <- rep(theta, each = n)

  if(!is.null(theta_se)) {
    theta_long <- rnorm(length(theta_long), theta_long, rep(theta_se, each = n))
  }
  theta_long
}


#--------------------------------------------------------------------------
#' Uses element-wise comparison of two vectors to create two new vectors, one with all the smaller values and the other with all the larger values.
#''
#' Optionally sorts an additional pair of vectors based on the element-wise comparison of the first two vectors.

#' @param theta1 the latent trait for member 1.
#' @param theta2 the latent trait for member 2.
#' @param theta1_se the standard error of the latent trait for member 1 (optional).
#' @param theta2_se the standard error of the latent trait for member 2 (optional).

#' @return A named list of length 4, containin \code{theta_min, theta_max, se_min, se_max}. The last two elements are null of the corresponding input is.
#' @export

theta_sort <- function(theta1, theta2, theta1_se = NULL, theta2_se = NULL) {
  n <- length(theta1)
  temp_theta <- cbind(theta1, theta2)
  min_ind <- cbind(1:n, apply(temp_theta, 1, which.min))
  max_ind <- cbind(1:n, 3 - min_ind[,2])
  out <- list(theta_min = temp_theta[min_ind], theta_max = temp_theta[max_ind])
  if (!is.null(theta1_se)) {
    temp_se <- cbind(theta1_se, theta2_se)
    out_se <- list(se_min = temp_se[min_ind], se_max = temp_se[max_ind])
  } else {
    out_se <- list(se_min = NULL, se_max = NULL)
  }
  c(out, out_se)
}




# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Old Plots ----------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------


#--------------------------------------------------------------------------
#' Plots individual versus collaborative peformance.
#'
#' Wrapper on ggplot to make barbell plots for pariwise collboration.
#'
#' @param ind_theta vector of test scores on a individual assessment.
#' @param col_theta corresponding vector of test scores on collaborative assessment.
#' @param group_score optional numeric vector used to color barbells; if omitted, each pair  has its own color.
#' @param legend passed to ggplot2 \code{legend.position}.
#' @return A barbell plot.
#' @export

barbell_plot <- function(ind_theta, col_theta, group_score = NULL, legend = "none") {
  data <- data.frame(ind_theta, col_theta)
  lim <- c(min(data)-.2, max(data)+.2)
  data$pairs <- factor(rep(1:(length(ind_theta)/2), each = 2))
  if (is.null(group_score)) {
    data$group_score <- data$pairs
    legend_title <- "pairs"
  } else {
    data$group_score <- group_score
    legend_title <- "group_score"
  }
  ggplot(data = data, aes(x = ind_theta, y = col_theta, group = pairs)) +
    geom_line(aes(color = group_score)) +
    geom_point(aes(color = group_score), size = 4) +
    scale_x_continuous(limits = lim) +
    scale_y_continuous(limits = lim) +
    geom_abline(intercept = 0, slope = 1, col = "grey") +
    theme(legend.position = legend) +
    labs(color = legend_title) +
    ggtitle("Collaborative vs Individual Performance") +
    xlab("Individual Theta")+
    ylab("Collaborative Theta")+
    theme(axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13)
   )
}

#--------------------------------------------------------------------------
#' Plots posterior probabilities of each model for each observation.
#'
#' Wrapper on ggplot.
#'
#' @param mix_prop a matrix of mixing proportions, where rows are obsevrations and columns are latent classes (e.g, the output of \code{em$poster})
#' @param sort logical indicating whether the observations should be sorted in expectation.
#' @param grey_scale logical indicating whether to use grey scale.

#' @return A raster plot.
#' @export

raster_plot <-function(mix_prop, sort = F, grey_scale = F) {
  mix_prop <- as.matrix(mix_prop)

  if (sort) {
    #u <- mix_prop[,4]
    u <- mix_prop%*%1:4
    mix_prop <- mix_prop[order(u, decreasing = F),]
  }

  temp <- data.frame(cbind(1:nrow(mix_prop), mix_prop))
  names(temp) <- c("pair", "Ind", "Min", "Max", "Add")

  gg <- reshape(temp,
    varying = names(temp)[-1],
    v.names = "prob",
    timevar = "model",
    times = names(temp)[-1],
    direction = "long"
  )
  gg$model <- ordered(gg$model, c("Ind", "Min", "Max", "Add"))

  #NYU <- rgb(87, 6, 140, maxColorValue = 255)
  # scale_fill_gradient2( high=muted('NYU'))

  p <- ggplot(gg, aes(pair, model, fill = prob)) +
    geom_raster() +
    theme_bw() +
    xlab("Group") +
    ylab("Model")

  if (grey_scale) {
    p <-  p + scale_fill_gradient(low="grey10", high="grey80")
  }
  p
}
