# Last update: 02/14/2017
# Functions for testing collaboration models. Depends on IRF_functions.R
# User beware: functions do not check or handle input errors.

source("~/github/cirt/R/IRF_functions.R")
require(ggplot2)
require(dplyr)


#--------------------------------------------------------------------------
#' Item response function for ``the Independence model"
#'
#' Computes a matrix of probabilities for correct responses using the independence model of collaboration for the members and the 2PL model for the items.
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
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
#' Computes a matrix of probabilities for correct responses using the minimum model of collaboration for the members and the 2PL model for the items.
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param sorted logical indicating whether to compute Min with theta1, regardless of the value of the theta2. Useful for data simulation, where variability in the simulated values of theta may not reflect the ordering of the data generating values.
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
#' Computes a matrix of probabilities for correct responses using the maximum model of collaboration for the members and the 2PL model for the items.
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait for member 1.
#' @param theta2 the latent trait for member 2.
#' @param sorted logical indicating whether to compute Max with theta2, regardless of the value of the theta2. Useful for data simulation, where variability in the simulated values of theta may not reflect the ordering of the data generating values.
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities.
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
#' Item response function for ``the Additive Independence model"
#'
#' Computes a matrix of probabilities for correct responses using the additive. independence model of collaboration for the members and the 2PL model for the items.
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait for member 1.
#' @param theta2 the latent trait for member 2.
#' @param sorted logical used for compatbility with \code{cIRF}.
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities.
#' @export

AI <- function(parms, theta1, theta2, sorted = F) {
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  p1 + p2 - p1*p2
}


#--------------------------------------------------------------------------
#' Wrapper for collaborative item response functions.
#'
#' Computes a matrix of probabilities for correct responses using the named \code{model} for collaboration, and the 2PL model for the items.
#'
#' @param model one of \code{c("Ind", "Min", "Max", "AI", "IRF")}
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait for member 1.
#' @param theta2 the latent trait for member 2.
#' @param sorted logical indicating whether to compute Min with theta1  and / or Max with theta2, regardless of the value of the theta1 and theta2. Useful for data simulation, where variability in the simulated values of theta may not reflect the ordering of the data generating values.
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities.
#' @export

cIRF <- function(model, parms, theta1, theta2, sorted = F) {
  if (model %in% c("Ind", "Min", "Max", "AI", "IRF")) {
    fun <- match.fun(model)
    fun(parms, theta1, theta2, sorted)
  } else {
    cat("\'model\' must be one of c(\"Ind\", \"Min\", \"Max\", \"AI\", \"IRF\")")
  }
}


#--------------------------------------------------------------------------
#' Computes a martix of item deltas.
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait for member 1.
#' @param theta2 the latent trait for member 2.
#' @param sorted logical indicating whether to compute Min with theta1  and / or Max with theta2, regardless of the value of the theta1 and theta2. Useful for data simulation, where variability in the simulated values of theta may not reflect the ordering of the data generating values.
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities.
#' @export

item_delta <- function(parms, theta1, theta2, sorted = F, NA_pattern = NULL) {
  temp <- Min(parms, theta1, theta2, sorted) * (1 - Max(parms, theta1, theta2, sorted))
  format_NA(temp, NA_pattern)
}


#--------------------------------------------------------------------------
#' Likelihood of a matrix of binary responses for one or more models of collaboration, conditional on theta.
#'
#' Mainly used to provide the component likelihoods for the finite mixture approach to model selection. Note that \code{logL} is faster for 2PL.
#'
#' @param resp the matrix of binary responses.
#' @param models is one or more of \code{c("IRF", "Ind", "Min", "Max", "AI") }.
#' @param parms a named list or data.frame with parms$alpha and parms$beta corresponding. to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait for member 1.
#' @param theta2 the latent trait for member 2.
#' @return An \code{nrow(resp)} by \code{length(models)} matrix of (log-) likleihoods for each response pattern and each model.
#' @export

likelihood <- function(models, resp, parms, theta1, theta2 = NULL, sorted = F, Log = T) {
  n_models <- length(models)
  out <- array(0, dim = c(nrow(resp), n_models))
  for (i in 1:n_models) {
    p <- cIRF(models[i], parms, theta1, theta2, sorted)
    if (Log) {
      out[,i] <- apply(log(p) * resp + log(1-p) * (1-resp), 1, sum, na.rm = T)
    } else {
      out[,i] <- apply(p * resp + (1-p) * (1-resp), 1, prod, na.rm = T)
    }
  }

  # for (i in 1:n_models) {
  #   p <- cIRF(models[i], parms, theta1, theta2, sorted)
  #
  #   out[,i] <- apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
  # }
  # if (Log) {out} else {exp(out)}

  if (n_models == 1) {out <- c(out)} # un-matrix
  out
}


#--------------------------------------------------------------------------
#' Incomplete data logliklihood for a mixture of collaboration models
#'
#' @param likelihood n_resp by n_models matrix of likelihoods (\strong{not loglikelihoods}) for each response pattern and each model (e.g., the output of \code{likelihood} with \code{Log = F}).
#' @param mix_prop the mixing proporitions for the models. Can be either an n_resp by n_models matrix (useful for computing the "posterior predicted" loglikelihood for each response pattern); or a n_models-vector, which is applied to each row of \code{likelihood} (useful for EM).
#' @param Sum logical indicating whether the output be summer over rows of \code{likelihood}.
#' @return A scalar (if \code{Sum = T}) or a n_resp-vector of incomplete data loglikelihoods.
#' @export

incomplete_data <- function(components, mix_prop, Sum = T) {
  if (!is.null(dim(mix_prop))) {
    temp <- components * mix_prop
  } else {
    temp <- t(t(components) * mix_prop)
  }
  out <- log(apply(temp, 1, sum))
  if(Sum) {sum(out)} else {out}
}


#--------------------------------------------------------------------------
#' Posterior probabilities of components in a mixture of collaboration models.
#'
#' This is the E-step of the EM algorithm for estimating the mixing proportions.
#'
#' @param components n_resp by n_models matrix of componentss (\strong{not logcomponentss}) for each response pattern and each model (e.g., the output of \code{components} with \code{Log = F}).
#' @param mix_prop a \code{length(models)}-vector of mixing proporitions for the models.
#' @return An n_resp by n_models matrix of posterior proabilities for each response pattern and each component.
#' @export

posterior <- function(components, mix_prop) {
  temp <- t(t(components) * mix_prop)
  temp / apply(temp, 1, sum)
}


#--------------------------------------------------------------------------
#' Computes updated mixing proportions based on posterior probabilities of components in a  mixture of collaboration models.
#'
#' This is the M-step of the EM algorithm for estimating the mixing proportions.
#''
#' @param post is output from \code{posterior}.
#' @return An n_models- vector mixing proportions each component.

prior <- function(post) {
   apply(post, 2, sum) / nrow(post)
}


#--------------------------------------------------------------------------
#' Computes standard errors of mixing proportions via observed Hessian.
#'
# Checked for accuracy again numDeriv::hessian.
#' @param components n_resp by n_models matrix of likehooods (\strong{not loglikelihoods}) for each response pattern and each model (e.g., the output of \code{likelihood} with \code{Log = F}).
#' @param mix_prop a n_models-vector of mixing proporitions for the models.

#' @return An n_models-vector standard errors for the mixing proportions in a mixture model.

prior_se <- function(components, mix_prop) {
  if (sum(mix_prop < .001) == 0) {
    temp <- apply(t(t(components) * mix_prop), 1, sum)
    temp <- components / temp
    H <- t(temp) %*% temp
    I <- solve(H)
    sqrt(diag(I))
  } else {
    rep(NA, times = length(mix_prop))
  }
}

#--------------------------------------------------------------------------
#' Runs EM algorithm for mixing proportions of a mixture of collaboration models
#'
#' Only the mixing proportions are estimated.
#'
#' @param models one or more of \code{c("Ind", "Min", "Max", "AI") }
#' @param resp data.frame of binary, conjunctively scored \strong{collaborative responses}
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param sorted logical indicating whether to compute Min with theta1  and / or Max with theta2, regardless of the value of the theta1 and theta2. Useful for data simulation, where variability in the simulated values of theta may not reflect the ordering of the data generating values.
#' @param max_iter maximum number of iterations for EM
#' @param conv convegence criterional applied to difference of incomplete data loglikelihood
#' @return A list of length 3 containing the optimization trace, priors, and posteriors
#' @export

EM <- function(models, resp, parms, theta1, theta2, sorted = F, max_iter = 100, conv = 1e-5) {
  n_models <- length(models)
  p <- rep(1/n_models, n_models)
  l <- likelihood(models, resp, parms, theta1, theta2, sorted, Log = F)
  trace <- incomplete_data(l, p)
  i <- 1
  delta <- 1

  while(i <= max_iter & delta > conv) {
    post <- posterior(l, p) # E
    p <- prior(post) # M
    trace <- c(trace, incomplete_data(l, p))
    delta <- trace[i+1] - trace[i]
    i <- i + 1
  }
  se <- prior_se(l, p)
  out <- list(trace, p, se, post)
  names(out) <- c("trace", "prior", "se", "posterior")
  out
}


#--------------------------------------------------------------------------
#' Formats a data.frame to have the same variables (columns) as that of another data.frame.
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
#' Simulate item responses from the 2PL or a model of pairwise collaboration obtained from the 2PL. Called by \code{data_gen}.

#' @param model is one of \code{c("IRF", "Ind", "Min", "Max", "AI") }.
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.
#' @param theta1 the latent trait for member 1.
#' @param theta2 the latent trait for member 2.
#' @return \code{length(theta1)} by \code{nrow(parms)} matrix of binary responses.
#' @export

sim_model <- function(model, parms, theta1 = 0, theta2 = 0, sorted = F, expected = F) {
  n_row <- length(theta1)
  n_col <- nrow(parms)
  r <- array(runif(n_row * n_col), dim = c(n_row, n_col))
  p <- cIRF(model, parms, theta1, theta2, sorted)
  if (expected) {
    out <- p
  } else {
    out <- ifelse(p > r, 1, 0)
  }
  colnames(out) <- row.names(parms)
  out
}


#--------------------------------------------------------------------------
#' Transforms one matrix/data.frame have the same NA entries as another matrix/data.frame.
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
#' Computes average classification probabilities for posterior distribution of mixture model
#'
#' @param mix_prop an n_obs by n_model matrix of probabilities, (e.g., \code{em$posterior}).
#' @param known_model an integer vector indicating the known class membership for each row of mix_prop. If omitted, it is internally replaced by \code{apply(em$posterior, 1, which.max)}

#' @return An n_model by n_model matrix with each row giving the mean posterior probabilities conditional on known_class.

#' @export
class_probs <-function(mix_prop, known_model = NULL){
  n_models <- ncol(mix_prop)
  out <- matrix(0, nrow = n_models, ncol = n_models)

  if (is.null(known_model)) {
    known_model <- apply(mix_prop, 1, which.max)
  }
  for(i in 1:n_models) {
    out[i, ] <-  apply(mix_prop[known_model == i, ], 2, mean)
  }
  row.names(out) <- paste0("model", 1:4)
  colnames(out) <- paste0("prob", 1:4)
  out
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

pv_gen <- function(n_reps, resp, parms, theta1, theta2, theta1_se, theta2_se, true_model = NULL) {

  # Expand data generating parms
  n_obs <- length(theta1)
  n_long <- n_obs * n_reps
  out <- data.frame(rep(1:n_obs, each = n_reps), rep(1:n_reps, times = n_obs))
  names(out) <- c("pairs", "samples")
  if (!is.null(true_model)) {out$model <- rep(true_model, each = n_reps) }

  # Sort thetas (for using sorted = T in post processing)
  temp_theta <- theta_sort(theta1, theta2, theta1_se, theta2_se)
  s_theta1 <- temp_theta$theta_min
  s_theta1_se <- temp_theta$se_min
  s_theta2 <- temp_theta$theta_max
  s_theta2_se <- temp_theta$se_max

  # Generate PVs for theta
  out$theta1 <- theta_gen(n_reps, s_theta1, s_theta1_se)
  out$theta2 <- theta_gen(n_reps, s_theta2, s_theta2_se)

  # Expand data (replicate each row n_reps times)
  data <- kronecker(as.matrix(resp), rep(1, n_reps))
  colnames(data) <- row.names(parms)
  cbind(out[], data[])
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

data_gen <- function(n_reps, mix_prop, parms, theta1, theta2, theta1_se = NULL, theta2_se = NULL, NA_pattern = NULL, expected = F, fixed_class = F) {

  # Expand data generating parms
  models <- c("Ind", "Min", "Max", "AI")
  n_obs <- length(theta1)
  n_long <- n_obs * n_reps
  mix_prop_long <- kronecker(as.matrix(mix_prop) , rep(1, n_reps))

  # Storage
  out <- data.frame(rep(1:n_obs, each = n_reps), rep(1:n_reps, times = n_obs))
  names(out) <- c("pairs", "samples")

  # Sort thetas
  temp_theta <- theta_sort(theta1, theta2, theta1_se, theta2_se)
  s_theta1 <- temp_theta$theta_min
  s_theta1_se <- temp_theta$se_min
  s_theta2 <- temp_theta$theta_max
  s_theta2_se <- temp_theta$se_max

  # Theta gen / replicate
  out$theta1 <- theta_gen(n_reps, s_theta1, s_theta1_se)
  out$theta2 <- theta_gen(n_reps, s_theta2, s_theta2_se)

  # Get model indices
  out[models] <- mix_prop_long

  # Get DGM
  if (fixed_class) {
    out$model <- apply(mix_prop, 1,
      function(x) sample.int(length(models), 1, prob = x)) %>% rep(each = n_reps)
  } else {
    out$model <- apply(mix_prop_long, 1,
      function(x) sample.int(length(models), 1, prob = x))
  }

  # Simulate data
  data <- data.frame(matrix(NA, nrow = n_long, ncol = nrow(parms)))
  names(data) <- row.names(parms)

  for (i in 1:length(models)) {
    temp <- out$model == i
    data[temp, ] <- sim_model(models[i], parms, out$theta1[temp], out$theta2[temp], sorted = T, expected = expected)
  }
  data <- format_NA(data, NA_pattern)
  cbind(out[], data[])
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
  names(temp) <- c("pair", "Ind", "Min", "Max", "AI")

  gg <- reshape(temp,
    varying = names(temp)[-1],
    v.names = "prob",
    timevar = "model",
    times = names(temp)[-1],
    direction = "long"
  )
  gg$model <- ordered(gg$model, c("Ind", "Min", "Max", "AI"))

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



# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# OLD CODE -----------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------


# #' Likelihood ratio tests for various models of collaboration.
# #'
# #' Computes a likelihood ratio test for one or more models of pairwise collaboration, given ``assumed to be known" item and person parameters (i.e., neither estimation error in item parameters nor prediction error in latent variables is accounted for by this procedure).
# #'
# #' @param resp the matrix binary data from the conjunctively scored \strong{collaborative responses}
# #' @param model is one or more of \code{c("Ind", "Min", "Max", "AI") }
# #' @param alpha the item discriminations of (only) the resp items
# #' @param beta the item difficulties of (only) the resp items
# #' @param ind_theta the \code{nrow(resp)*2}-dimensional vector of latent traits for each member, as estimated from a non-collaborative form
# #' @param col_theta the \code{nrow(resp)}-dimensional vector of latent traits for each pair, as estimated from a (conjunctively scored) collaborative form
# #' @param n_boot number of bootstraps to use for testing the likelihood ratio.
# #' @return A list of length \code{length(model)}, each element of which is a data frame with \code{nrow(resp)} rwos containing the output for lr_tests for each pair.
# #' @export
#
# lr_test_old <-function(resp, model, alpha, beta, ind_theta, col_theta, n_boot = 0) {
#
#   odd <- seq(from = 1, to = length(ind_theta), by = 2)
#   theta1 <- ind_theta[odd]
#   theta2 <- ind_theta[odd+1]
#   n_pair <- length(odd)
#   n_model <- length(model)
#   mod <- matrix(0, nrow = n_pair, ncol = n_model)
#   out <- vector("list", n_model)
#   names(out) <- model
#
#   # Helper function for computing P(x > obs)
#   pobs <- function(cdf, obs) {
#     1 - environment(cdf)$y[which.min(abs(environment(cdf)$x-obs))]
#   }
#
#   # logL for collaboration models
#   for (i in 1:n_model) {
#     mod[, i] <-
#       logL(resp, model = model[i], alpha, beta, theta1, theta2)
#   }
#
#   # logL for reference model
#   ref <- logL(resp, "IRF", alpha, beta, col_theta)
#   lr <- -2*(mod - ref)
#
#   # Bootstrapping (could fancy this up...)
#   if (n_boot > 0) {
#     theta1_long <- rep(theta1, each = n_boot)
#     theta2_long <- rep(theta2, each = n_boot)
#     boot_ind <- rep(1:n_pair, each = n_boot)
#
#     for (i in 1:n_model) {
#       message(cat("Running bootstraps for model", i, "..."),"\r", appendLF = FALSE)
#       flush.console()
#
#       boot_data <- sim_data(model[i], alpha, beta, theta1_long, theta2_long)
#       boot_mod <- logL(boot_data, model[i], alpha, beta, theta1_long, theta2_long)
#       boot_ref <- ml_IRF(boot_data, alpha, beta)
#       boot_lr <- -2*(boot_mod - boot_ref[,1])
#
#       # 95% CIs
#       temp <- tapply(boot_lr, boot_ind, function(x) quantile(x, p = c(.025, .975)))
#       boot_ci <- t(matrix(unlist(temp), nrow = 2, ncol = n_pair))
#
#       # P(lr > obs)
#       boot_cdf <- tapply(boot_lr, boot_ind, ecdf)
#       boot_p <-mapply(function(x,y) pobs(x, y), boot_cdf, lr[,i])
#
#       # Storage
#       temp <- data.frame(cbind(lr[,i], boot_ci[], boot_p))
#       names(temp) <- c("lr", "ci_lower", "ci_upper", "p_obs")
#       out[[i]] <- temp
#
#     }
#   }
#   out
# }
#
#
#
#
# #--------------------------------------------------------------------------
# #' Simulates data from an averaged model of collaboration resulting from application of \code{EM}.
# #'
# #' Generates data from an averaged model of pairwise collaboration, for one or more dyads indexed by theta1 and theta2. For each dyad, arppoximately n_i = \code{mix_prop[i] * n_boot} response patterns are generated from each of the i = 1,..4 models of collaboration. If \code{theta_se} are included, data generation uses a plausible values approach in which \code{n_boot} values of theta are generated using \code{rnorm(n_boot, theta, theta_se)}, for each dayd.
# #'
# #' @param n_boot number of samples to generate for each dyad
# #' @param mix_prop is em$posterior
# #' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
# #' @param theta1 the latent trait for member 1
# #' @param theta2 the latent trait for member 2
# #' @param theta1_se the standard error of the latent trait for member 1 (optional)
# #' @param theta2_se the standard errr of the latent trait for member 2 (optional)
# #' @param NA_pattner an (optional) \code{length(theta)} by \code{nrow(parms)} data.frame with \code{NA} entries for each item a dyad did not answer. The missing values are preserved in the generated data. This would usually be the original response matrix.
#
# #' @return A data.frame with \code{length(theta)*n_boot} rows containing an id variable for each pair, the data generating values of theta1, theta2, and mix_prop; the model used to simulate the response pattern; and the simulated response pattern.
# #' @export
#
# theta_gen <- function(n, theta, theta_se = NULL){
#   theta_long <- rep(theta, each = n)
#
#   if(!is.null(theta_se)) {
#     theta_long <- rnorm(n, theta_long, rep(theta_se, each = n))
#   }
#   theta_long
# }
#
# sim_data_old <- function(n_obs, mix_prop, parms, theta1 = 0, theta2 = 0, theta1_se = NULL, theta2_se = NULL, NA_pattern = NULL) {
#
#   # Expand data generating parms
#   models <- c("Ind", "Min", "Max", "AI")
#   n_long <- length(theta1)*n_obs
#   pairs_long <- rep(1:length(theta1), each = n_obs)
#   mix_prop_long <- kronecker(mix_prop, rep(1, n_obs))
#
#   # Step 1. Use PV for theta if SE given
#   theta1_long <- rep(theta1, each = n_obs)
#   theta2_long <- rep(theta2, each = n_obs)
#
#   if(!is.null(theta1_se)) {
#     theta1_long <- rnorm(n_long, theta1_long, rep(theta1_se, each = n_obs))
#   }
#   if(!is.null(theta2_se)) {
#     theta2_long <- rnorm(n_long, theta2_long, rep(theta2_se, each = n_obs))
#   }
#
#   # Set up storage
#   out <- data.frame(pairs_long, theta1_long, theta2_long, mix_prop_long)
#   names(out) <- c("pairs", "theta1", "theta2", models)
#
#   # Step 2. Get model for each rep
#   out$model <- model_indices(n_long, mix_prop_long)
#
#   # Step 3. Simulate data
#   data <- data.frame(matrix(NA, nrow = n_long, ncol = nrow(parms)))
#   names(data) <- row.names(parms)
#
#   for (i in 1:length(models)) {
#     temp <- out$model == i
#     data[temp, ] <- sim_model(models[i], parms, out$theta1[temp], out$theta2[temp])
#   }
#   NA_pattern_long <- kronecker(as.matrix(NA_pattern), rep(1, n_boot))
#   data <- format_NA(data, NA_pattern_long)
#   cbind(out[], data[])
# }
#
#
# #--------------------------------------------------------------------------
# #' Bootstrapped likelihood ratio test for averaged model of collaboration resulting from application of \code{EM}.
# #'
# #' Computes a likelihood ratio test for the averaged model of pairwise collaboration, given ``assumed to be known" item parameters (i.e., estimation error in item parameters is not accounted for by this procedure). If SEs are included for the thetas, data generation uses a plausible values approach in which \code{n_obs} values of theta are generated using \code{rnorm(n_boot, theta, theta_se)}, for each dayd (i.e., prediction error in theta1 and theta2 can be accounted for using this procedure).
# #'
# #' @param resp the matrix binary data from the conjunctively scored \strong{collaborative responses}
# #' @param mix_prop is the em$posterior
# #' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
# #' @param theta1 the latent trait for member 1
# #' @param theta2 the latent trait for member 2
# #' @param theta1_se the standard error of the latent trait for member 1
# #' @param theta2_se the standard errr of the latent trait for member 2
# #' @param n_boot number of bootstraps to use for testing the likelihood ratio.
# #' @return A data.frame with \code{nrow(resp)} rows containing the output for lr_tests for each pair.
# #' @export
#
# lr_test <- function(resp, mix_prop, parms, theta1 = 0, theta2 = 0, theta1_se = NULL, theta2_se = NULL, n_boot = 0) {
#
#   models <- c("Ind", "Min", "Max", "AI")
#   items <- row.names(parms)
#
#   # Helper function for computing P(x > obs)
#   pobs <- function(cdf, obs) {
#     1 - environment(cdf)$y[which.min(abs(environment(cdf)$x-obs))]
#   }
#
#   # Observed likelihood ratios for each dyad
#   components <- likelihood(models, resp, parms, theta1, theta2, Log = F)
#   log_mix <- incomplete_data(components, mix_prop, Sum = F)
#   mle_ref <- MLE(resp, parms)
#   lr_obs <- -2*(log_mix - mle_ref$logL)
#
#   # Bootstrapping
#   if (n_boot == 0) {
#     return(lr_obs)
#   } else {
#     boot <- sim_mix(100, mix_prop, parms, theta1, theta2, theta1_se, theta2_se)
#     temp <- likelihood(models, boot[items], parms, boot$theta1, boot$theta2, Log = F)
#     boot$log_mix <- incomplete_data(temp, boot[models], Sum = F)
#     temp <- MLE(boot[items], parms)
#     boot <- cbind(boot, temp)
#     boot$lr <- -2*(boot$log_mix - boot$logL)
#
#     # 95% CIs
#     temp <- tapply(boot$lr, boot$pairs, function(x) quantile(x, p = c(.025, .975)))
#     boot_ci <- t(matrix(unlist(temp), nrow = 2, ncol = length(theta1)))
#
#     # P(lr > obs)
#     boot_cdf <- tapply(boot$lr, boot$pairs, ecdf)
#     boot_p <- mapply(function(x,y) pobs(x, y), boot_cdf, lr_obs)
#
#     # Storage
#     out <- data.frame(cbind(lr_obs, boot_ci, boot_p))
#     names(out) <- c("lr", "ci_lower", "ci_upper", "p_obs")
#   }
#   out
# }
