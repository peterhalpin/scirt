# Last update: 07/12/2016
# Functions for estimating and testing collaboration models. Depends on IRF_functions.R
# User beware: functions not written to check or handle input errors.
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
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities.
#' @export

Ind <- function(parms, theta1, theta2) {
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
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities.
#' @export

Min <- function(parms, theta1, theta2) {
  theta <- apply(cbind(theta1, theta2), 1, min, na.rm = T)
  IRF(parms, theta)
}


#--------------------------------------------------------------------------
#' Item response function for ``the Maximum model"
#'
#' Computes a matrix of probabilities for correct responses using the maximum model of collaboration for the members and the 2PL model for the items.
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

Max <- function(parms, theta1, theta2) {
  theta <- apply(cbind(theta1, theta2), 1, max, na.rm = T)
  IRF(parms, theta)
}


#--------------------------------------------------------------------------
#' Item response function for ``the Additive Independence model"
#'
#' Computes a matrix of probabilities for correct responses using the additive independence model of collaboration for the members and the 2PL model for the items.
#'
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

AI <- function(parms, theta1, theta2) {
  p1 <- IRF(parms, theta1)
  p2 <- IRF(parms, theta2)
  p1 + p2 - p1*p2
}


#--------------------------------------------------------------------------
#' Wrapper for collaborative item response functions
#'
#' Computes a matrix of probabilities for correct responses using the named \code{model} for collaboration, and the 2PL model for the items.
#' @param model
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta)} by \code{nrow(parms)} matrix of response probabilities
#' @export

cIRF <- function(model, parms, theta1, theta2) {
  if (model %in% c("Ind", "Min", "Max", "AI", "IRF")) {
    fun <- match.fun(model)
    fun(parms, theta1, theta2)
  } else {
    cat("\'model\' must be one of c(\"Ind\", \"Min\", \"Max\", \"AI\", \"IRF\")")
  }
}


#--------------------------------------------------------------------------
#' Likelihood of a matrix of binary responses for one or more models, conditional on theta.
#'
#' \code{logL} is faster for 2PL. Mainly used to provide the component likelihoods for the finite mixture approach to selection / averaging models of collaboration.
#'
#' @param resp the matrix binary responses
#' @param models is one or more of \code{c("IRF", "Ind", "Min", "Max", "AI") }
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return An \code{nrow(resp)} by \code{length(models)} matrix of log-likleihoods for each response pattern and each model
#' @export

likelihood <- function(models, resp, parms, theta1, theta2 = NULL, Log = T) {
  n_models <- length(models)
  out <- array(0, dim = c(nrow(resp), n_models))
  for (i in 1:n_models) {
    p <- cIRF(models[i], parms, theta1, theta2)
    out[,i] <- apply(log(p) * (resp) + log(1-p) * (1-resp), 1, sum, na.rm = T)
  }
  if (n_models == 1) {out <- c(out)} # un-matrix
  if (Log) {out} else {exp(out)}
}


#--------------------------------------------------------------------------
#' Incomplete data logliklihood for a mixture of collaboration models
#'
#' @param components n_resp by n_models matrix of likelihoods (\strong{not loglikelihoods}) for each response pattern and each model (e.g., the output of \code{likelihood} with \code{Log = F})
#' @param mix_prop the mixing proporitions for the models. Can be either an n_resp by n_models matrix (useful for computing the "posterior predicted" loglikelihood for each response pattern); or a n_models-vector, which is applied to each row of \code{components} (useful for EM).
#' @param Sum should the output be summer over rows of \code{components} ?
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
#' @param components n_resp by n_models matrix of likelihoods (\strong{not loglikelihoods}) for each response pattern and each model (e.g., the output of \code{likelihood} with \code{Log = F})
#' @param prior the mixing proporitions for the models. Must be a \code{length(models)}-vector, which is applied to each row of \code{components}
#' @return An n_resp by n_models matrix of posterior proabilities for each response pattern and each component.

#' @export

posterior <- function(components, prior) {
  temp <- t(t(components) * prior)
  temp / apply(temp, 1, sum)
}


#--------------------------------------------------------------------------
#' Computes updated mixing proportions based on posterior probabilities of components in a  mixture of collaboration models
#'
#' This is the M-step of the EM algorithm for estimating the mixing proportions.
#''
#' @param posterior is output from \code{posterior}
#' @param prior the mixing proporitions for the models. Must be a \code{length(models)}-vector, which is applied to each row of \code{components}
#' @return An n_models- vector mixing proportions each component.

prior <- function(posterior) {
   apply(posterior, 2, sum) / nrow(posterior)
}


#--------------------------------------------------------------------------
#' Runs EM algorithm for mixing proportions for a mixture of collaboration models
#'
#' Only the mixing proportions are estimated. Multiple starts are not used becase the Q function is concave in mixing proportions.
#'
#' @param models one or more of \code{c("Ind", "Min", "Max", "AI") }
#' @param resp data.frame of binary, conjunctively scored \strong{collaborative responses}
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param max_iter maximum number of iterations for EM
#' @param conv convegence criterional applied to difference of incomplete data loglikelihood
#' @return A list of length 3 containing the optimization trace, priors, and posteriors
#' @export

EM <- function(models, resp, parms, theta1, theta2, max_iter = 100, conv = 1e-5) {
  n_models <- length(models)
  p <- rep(1/n_models, n_models)
  l <- likelihood(models, resp, parms, theta1, theta2, Log = F)

  trace <- incomplete_data(l, p)
  i <- 1
  delta <- 1

  while(i <= max_iter & delta > conv) {
    post <- posterior(l, p) # Eeeeee
    p <- prior(post) # Mmmmmm
    trace <- c(trace, incomplete_data(l, p))
    delta <- trace[i+1] - trace[i]
    i <- i + 1
  }
  out <- list(trace, p, post)
  names(out) <- c("trace", "prior", "posterior")
  out
}


#--------------------------------------------------------------------------
#' Formats responses a response data.frame to match the calibration data.frame
#'
#' Drops items in the target df not in the calibration sample; adds items (with NA entries) in the calibration df not in the target df. This simplifies using item parms obtained from the calibration df with the target df

#' @param resp the target df to be formatted
#' @param item string vector of names of items in the calibration sample
#' @param version an optional string used to subset resp via \code{grep(version, names(resp))}
#' @return a version of resp that has the same names as calib
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
#' Simulate data from the 2PL or a model of pairwise collaboration.
#''
#' Simulate data using either the 2PL or a model for pariwise collaboration obtained from the 2PL.
#' @param model is one of \code{c("IRF", "Ind", "Min", "Max", "AI") }
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta1)} by \code{nrow(parms)} matrix of binary response patterns
#' @export

sim_data <- function(model, parms, theta1 = 0, theta2 = 0) {
  n_row <- length(theta1)
  n_col <- nrow(parms)
  fun <- match.fun(model)
  r <- array(runif(n_row * n_col), dim = c(n_row, n_col))
  p <- cIRF(model, parms, theta1, theta2)
  out <- ifelse(p > r, 1, 0)
  colnames(out) <- row.names(parms)
  out
}


#--------------------------------------------------------------------------
#' Used by sim_mix to write out the appropriate number of model labels for each dyad.
#'
#' The desried output is to replicate each of \code{c("Ind", "Min", "Max", "AI")} n_i = mix_prop[i] * n_boot times. This function (badly) handles rounding error when computing the n_i.
#'
#' @param mix_prop is em$posterior
#' @param n_boot is the desired number of replications of each row of mix prop

#' @return An \code{nrow(miz_prop) * n_boot} vector that replicates \code{c("Ind", "Min", "Max", "AI")} according to the number simulated data sets desired for each model.
#' @export

model_indices <- function(mix_prop, n_boot){
  indices <-  mix_prop * n_boot
  dif <- apply(round(indices), 1, sum) - n_boot
  temp <- apply((indices %% 1)*10, 1, order, decreasing = T)
  add_ind <- cbind(1:nrow(mix_prop), temp[2,])
  sub_ind <- cbind(1:nrow(mix_prop), temp[1,])
  indices[sub_ind[dif > 0,]] <- indices[sub_ind[dif > 0,]] - dif[dif > 0]
  indices[add_ind[dif < 0,]] <- indices[add_ind[dif < 0,]] - dif[dif < 0]
  models_long <- rep(models, times = length(theta1))
  rep(models_long, c(round(t(indices))))
}


#--------------------------------------------------------------------------
#' Simulates data from an averaged model of collaboration resulting from application of \code{EM}.
#'
#' Generates data from n averaged model of pairwise collaboration, for one or more dyads indexed by theta1 and theta2. For each dyad, n_i \code{mix_prop[i] * n_boot} response patterns are generated from each of the i = 1,..4 models of collaboration. If SEs are included, data generation uses a plausible values approach in which \code{n_boot} values of theta are generated using \code{rnorm(n_boot, theta, theta_se)}, for each dayd.
#'
#' @param n_boot number of samples to generate for each dyad
#' @param mix_prop is em$posterior
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param theta1_se the standard error of the latent trait for member 1
#' @param theta2_se the standard errr of the latent trait for member 2

#' @return A data.frame with \code{length(theta)*n_boot} rows containing an id variable for for each pair, the data generating values of theta1, theta2, and mix_prop; the model used to simulate the response pattern; and the simulated response pattern.
#' @export

sim_mix <- function(n_boot, mix_prop, parms, theta1 = 0, theta2 = 0, theta1_se = NULL, theta2_se = NULL) {

  # Expand data generating parms
  models <- c("Ind", "Min", "Max", "AI")
  pairs_long <- rep(1:length(theta1), each = n_boot)
  theta1_long <- rep(theta1, each = n_boot)
  theta2_long <- rep(theta2, each = n_boot)
  mix_prop_long <- kronecker(mix_prop, rep(1, n_boot))

  # Use PV for theta if SE given
  if(!is.null(theta1_se)) {
    theta1_long <- rnorm(length(theta1_long), theta1_long, rep(theta1_se, each = n_boot))
  }
  if(is.null(theta2_se)) {
    theta2_long <- rnorm(length(theta2_long), theta2_long, rep(theta2_se, each = n_boot))
  }

  # Set up output
  out <- data.frame(pairs_long, theta1_long, theta2_long, mix_prop_long)
  names(out) <- c("pairs", "theta1", "theta2", models)
  out$model <- model_indices(mix_prop, n_boot)

  # Simulate data
  data <- data.frame(matrix(NA, nrow = nrow(out), ncol = nrow(parms)))
  names(data) <- row.names(parms)

  for (i in models) {
    temp <- out$model == i
    data[temp, ] <- sim_data(i, parms, out$theta1[temp], out$theta2[temp])
  }
  cbind(out[], data[])
}


#--------------------------------------------------------------------------
#' Bootstrapped likelihood ratio test for averaged model of collaboration resulting from application of \code{EM}.
#'
#' Computes a likelihood ratio test for the averaged model of pairwise collaboration, given ``assumed to be known" item parameters (i.e., estimation error in item parameters is not accounted for by this procedure). If SEs are included for the thetas, data generation uses a plausible values approach in which \code{n_boot} values of theta are generated using \code{rnorm(n_boot, theta, theta_se)}, for each dayd (i.e., prediction error in theta1 and theta2 can be accounted for using this procedure).
#'
#' @param resp the matrix binary data from the conjunctively scored \strong{collaborative responses}
#' @param mix_prop is the em$posterior
#' @param parms a list or data.frame with elements parms$alpha and parms$beta corresponding to the discrimination and difficulty parameters of the 2PL model, respectively
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @param theta1_se the standard error of the latent trait for member 1
#' @param theta2_se the standard errr of the latent trait for member 2
#' @param n_boot number of bootstraps to use for testing the likelihood ratio.
#' @return A data.frame with \code{nrow(resp)} rows containing the output for lr_tests for each pair.
#' @export

lr_test <- function(resp, mix_prop, parms, theta1 = 0, theta2 = 0, theta1_se = NULL, theta2_se = NULL, n_boot = 0) {

  models <- c("Ind", "Min", "Max", "AI")
  items <- row.names(parms)

  # Helper function for computing P(x > obs)
  pobs <- function(cdf, obs) {
    1 - environment(cdf)$y[which.min(abs(environment(cdf)$x-obs))]
  }

  # Observed likelihood ratios for each dyad
  components <- likelihood(models, resp, parms, theta1, theta2, Log = F)
  log_mix <- incomplete_data(components, mix_prop, Sum = F)
  mle_ref <- MLE(resp, parms)
  lr_obs <- -2*(log_mix - mle_ref$logL)

  # Bootstrapping
  if (n_boot == 0) {
    return(lr_obs)
  } else {
    boot <- sim_mix(100, mix_prop, parms, theta1, theta2, theta1_se, theta2_se)
    temp <- likelihood(models, boot[items], parms, boot$theta1, boot$theta2, Log = F)
    boot$log_mix <- incomplete_data(temp, boot[models], Sum = F)
    temp <- MLE(boot[items], parms)
    boot <- cbind(boot, temp)
    boot$lr <- -2*(boot$log_mix - boot$logL)

    # 95% CIs
    temp <- tapply(boot$lr, boot$pairs, function(x) quantile(x, p = c(.025, .975)))
    boot_ci <- t(matrix(unlist(temp), nrow = 2, ncol = length(theta1)))

    # P(lr > obs)
    boot_cdf <- tapply(boot$lr, boot$pairs, ecdf)
    boot_p <- mapply(function(x,y) pobs(x, y), boot_cdf, lr_obs)

    # Storage
    out <- data.frame(cbind(lr_obs, boot_ci, boot_p))
    names(out) <- c("lr", "ci_lower", "ci_upper", "p_obs")
  }
  out
}





# --- under construction / undocumented functions ----


delta <- function(alpha, beta, theta1, theta2) {
  Min(alpha, beta, theta1, theta2) * (1 - Max(alpha, beta, theta1, theta2))
}

screen <- function(cutoff, alpha, beta, theta1, theta2) {
  screen <- delta(alpha, beta, theta1, theta2)
  screen[screen < cutoff] <- NA
  screen[!is.na(screen)] <- 1
  screen
}

raster_plot <-function(em, sort = F) {
  temp <- em$posterior
  u <- temp%*%1:4
  if (sort) {
      temp <- temp[order(u, decreasing = F),]
  }
  temp <- data.frame(cbind(1:nrow(temp), temp))
  names(temp) <- c("pair", "Ind", "Min", "Max", "AI")


  q <- reshape(temp,
    varying = names(temp)[-1],
    v.names = "prob",
    timevar = "model",
    times = names(temp)[-1],
    direction = "long"
    )
  q$model <- ordered(q$model, c("Ind", "Min", "Max", "AI"))

  #NYU <- rgb(87, 6, 140, maxColorValue = 255)
  # scale_fill_gradient2( high=muted('NYU'))
  ggplot(q, aes(pair, model, fill = prob)) + geom_raster() + theme_bw()
  #   scale_fill_gradient2(high = NYU) +theme_bw()
}

class_accuracy <- function(em) {
  ind <- apply(em$posterior, 1, which.max)
  arr_ind <- cbind(1:nrow(em$posterior), ind)
  cp <- em$posterior[arr_ind]
  tapply(cp, ind, mean)
}


#--------------------------------------------------------------------------
#' Plots individual versus collaborative peformance
#'
#' Wrapper on ggplot to make barbell plots for pariwise collboration.
#'
#' @param ind_theta vector of test scores on a individual assessment
#' @param col_theta corresponding vector of test scores on collaborative assessment
#' @param group_score optional numeric vector used to color barbells; if omitted, each pair  has its own color
#' @param legend passed to ggplot2 \code{legend.position}
#' @return A barbell plot
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






sim_em <- function(n_obs, n_items, prior = NULL, alpha = NULL, beta = NULL, sort = T) {
  temp1 <- rnorm(n_obs)
  temp2 <- rnorm(n_obs)
  theta1 <- apply(cbind(temp1, temp2), 1, max)
  theta2 <- apply(cbind(temp1, temp2), 1, min)
  if (is.null(alpha)) { alpha <- rep(1, times = n_items) }
  if (is.null(beta)) { beta <- sort(runif(n_items, -3, 3)) }
  if (is.null(prior)) { prior <- c(.25, .25, .25, .25) }
  parms <- data.frame(alpha, beta)
  out <- matrix(NA, nrow = n_obs, ncol = n_items)
  n <- c(0, round(prior * n_obs))
  models <- c("Ind", "Min", "Max", "AI")

  for(i in 1:length(models)) {
    j <- sum(n[1:i]) + 1
    k <- sum(n[1:(i+1)])
    if (j < k) {
      out[j:k, ] <-  sim_data(models[i], parms, theta1[j:k], theta2[j:k])
    }
  }
  r <- 1:nrow(out)
  if (sort) r <- sample(1:nrow(out))
  EM(models, out[r,], parms, theta1[r], theta2[r])
}





#--------------------------------------------------------------------------
#' Likelihood ratio tests for various models of collaboration.
#'
#' Computes a likelihood ratio test for one or more models of pairwise collaboration, given ``assumed to be known" item and person parameters (i.e., neither estimation error in item parameters nor prediction error in latent variables is accounted for by this procedure).
#'
#' @param resp the matrix binary data from the conjunctively scored \strong{collaborative responses}
#' @param model is one or more of \code{c("Ind", "Min", "Max", "AI") }
#' @param alpha the item discriminations of (only) the resp items
#' @param beta the item difficulties of (only) the resp items
#' @param ind_theta the \code{nrow(resp)*2}-dimensional vector of latent traits for each member, as estimated from a non-collaborative form
#' @param col_theta the \code{nrow(resp)}-dimensional vector of latent traits for each pair, as estimated from a (conjunctively scored) collaborative form
#' @param n_boot number of bootstraps to use for testing the likelihood ratio.
#' @return A list of length \code{length(model)}, each element of which is a data frame with \code{nrow(resp)} rwos containing the output for lr_tests for each pair.
#' @export

lr_test_old <-function(resp, model, alpha, beta, ind_theta, col_theta, n_boot = 0) {

  odd <- seq(from = 1, to = length(ind_theta), by = 2)
  theta1 <- ind_theta[odd]
  theta2 <- ind_theta[odd+1]
  n_pair <- length(odd)
  n_model <- length(model)
  mod <- matrix(0, nrow = n_pair, ncol = n_model)
  out <- vector("list", n_model)
  names(out) <- model

  # Helper function for computing P(x > obs)
  pobs <- function(cdf, obs) {
    1 - environment(cdf)$y[which.min(abs(environment(cdf)$x-obs))]
  }

  # logL for collaboration models
  for (i in 1:n_model) {
    mod[, i] <-
      logL(resp, model = model[i], alpha, beta, theta1, theta2)
  }

  # logL for reference model
  ref <- logL(resp, "IRF", alpha, beta, col_theta)
  lr <- -2*(mod - ref)

  # Bootstrapping (could fancy this up...)
  if (n_boot > 0) {
    theta1_long <- rep(theta1, each = n_boot)
    theta2_long <- rep(theta2, each = n_boot)
    boot_ind <- rep(1:n_pair, each = n_boot)

    for (i in 1:n_model) {
      message(cat("Running bootstraps for model", i, "..."),"\r", appendLF = FALSE)
      flush.console()

      boot_data <- sim_data(model[i], alpha, beta, theta1_long, theta2_long)
      boot_mod <- logL(boot_data, model[i], alpha, beta, theta1_long, theta2_long)
      boot_ref <- ml_IRF(boot_data, alpha, beta)
      boot_lr <- -2*(boot_mod - boot_ref[,1])

      # 95% CIs
      temp <- tapply(boot_lr, boot_ind, function(x) quantile(x, p = c(.025, .975)))
      boot_ci <- t(matrix(unlist(temp), nrow = 2, ncol = n_pair))

      # P(lr > obs)
      boot_cdf <- tapply(boot_lr, boot_ind, ecdf)
      boot_p <-mapply(function(x,y) pobs(x, y), boot_cdf, lr[,i])

      # Storage
      temp <- data.frame(cbind(lr[,i], boot_ci[], boot_p))
      names(temp) <- c("lr", "ci_lower", "ci_upper", "p_obs")
      out[[i]] <- temp

    }
  }
  out
}
