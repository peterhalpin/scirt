# Last update: 04/12/2016
# Functions for computing bootstrapped LR tests of collaboration models.
# User beware: functions not written to check or handle input errors.

require(stats4)
require(ggplot2)
require(dplyr)

#--------------------------------------------------------------------------
#' Item response function for 2PL
#'
#' Computes a matrix of probabilities for correct responses under 2PL model
#'
#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta the latent trait
#' @return \code{length(theta)} by \code{length(alpha)} matrix of response probabilities
#' @export

twoPL <-function(alpha, beta, theta){
  Z <- matrix(0, nrow = length(theta), ncol = length(alpha))
  Z <- Z + theta
  Z <- t(alpha*(t(Z) - beta))
  1/(1 + exp(-Z))
}


#--------------------------------------------------------------------------
#' Item response function for ``the Independence model"
#'
#' Computes a matrix of probabilities for correct responses using the independence model of collaboration for the members and the 2PL model for the items.
#'
#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta)} by \code{length(alpha)} matrix of response probabilities.
#' @export

Ind <- function(alpha, beta, theta1, theta2){
  twoPL(alpha, beta, theta1)*twoPL(alpha, beta, theta2)
}


#--------------------------------------------------------------------------
#' Item response function for ``the Minimum model"
#'
#' Computes a matrix of probabilities for correct responses using the minimum model of collaboration for the members and the 2PL model for the items.
#'
#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta)} by \code{length(alpha)} matrix of response probabilities.
#' @export

Min <- function(alpha, beta, theta1, theta2){
  theta <- apply(cbind(theta1, theta2), 1, min, na.rm = T)
  twoPL(alpha, beta, theta)
}


#--------------------------------------------------------------------------
#' Item response function for ``the Maximum model"
#'
#' Computes a matrix of probabilities for correct responses using the maximum model of collaboration for the members and the 2PL model for the items.
#'
#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta)} by \code{length(alpha)} matrix of response probabilities
#' @export

Max <- function(alpha, beta, theta1, theta2){
  theta <- apply(cbind(theta1, theta2), 1, max, na.rm = T)
  twoPL(alpha, beta, theta)
}


#--------------------------------------------------------------------------
#' Item response function for ``the Additive Independence model"
#'
#' Computes a matrix of probabilities for correct responses using the additive independence model of collaboration for the members and the 2PL model for the items.
#'
#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta)} by \code{length(alpha)} matrix of response probabilities
#' @export

AI <- function(alpha, beta, theta1, theta2){
  twoPL(alpha, beta, theta1) + twoPL(alpha, beta, theta2) -  twoPL(alpha, beta, theta1) * twoPL(alpha, beta, theta2)
}


#--------------------------------------------------------------------------
#' Simulate data from a the 2PL or a model of pairwise collaboration.
#'
#' Simulate data using either the 2PL or a model for pariwise collaboration obtained from the 2PL.
#' @param model is one of \code{c("twoPL", "Ind", "Min", "Max", "AI") }
#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta1)} by \code{length(alpha)} matrix of binary response patterns
#' @export

sim_data <- function(model, alpha, beta, theta1 = 0, theta2 = 0){
  n_row <- length(theta1)
  n_col <- length(alpha)
  fun <- match.fun(model)
  Q <- array(runif (n_row * n_col), dim = c(n_row, n_col))

  if (model == "twoPL"){
    P <- fun(alpha, beta, theta1)
  } else {
    P <- fun(alpha, beta, theta1, theta2)
  }

  OUT <- ifelse (P > Q, 1, 0)
  colnames(OUT) <- names(alpha)
  OUT
}

sim_data <- function(model, alpha, beta, theta1 = 0, theta2 = 0){
  n_row <- length(theta1)
  n_col <- length(alpha)
  fun <- match.fun(model)
  Q <- array(runif (n_row * n_col), dim = c(n_row, n_col))

  if (model == "twoPL"){
    P <- fun(alpha, beta, theta1)
  } else {
    P <- fun(alpha, beta, theta1, theta2)
  }

  OUT <- ifelse (P > Q, 1, 0)
  colnames(OUT) <- names(alpha)
  OUT
}

#--------------------------------------------------------------------------
#' Log-ikelihood of given matrix of binary responses for a stated model, conditional on theta.
#'
#' Computes a vector of likelihoods for a given matrix binary response patterns, for a selection of models for pairwise collaboration, conditional on theta
#'
#' @param resp the matrix binary responses
#' @param model is one of \code{c("twoPL", "Ind", "Min", "Max", "AI") }
#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return An \code{nrow(resp)}-vector of log-likleihoods for each response pattern.
#' @export


logL <- function(resp, model, alpha, beta, theta1, theta2 = NULL){
  if (model == "twoPL"){
    p <- twoPL(alpha, beta, theta1)
  } else {
    fun <- match.fun(model)
    p <- fun(alpha, beta, theta1, theta2)
  }
  apply(log(p)*(resp) + log(1-p)*(1-(resp)), 1, sum, na.rm = T)
}


#--------------------------------------------------------------------------
#' Internal function used in lr_test; under development.
#' @export

neg_logL <- function(theta, resp, alpha, beta){
  -1*logL(resp, "twoPL", alpha, beta, theta)
}



#--------------------------------------------------------------------------
#' Internal function used in lr_test; under development.
#' @export

ml_twoPL<-function(resp, alpha, beta, method = "ML"){
  OUT <- matrix(0, nrow = nrow(resp), ncol = 3)
  colnames(OUT) <- c("logL", "theta", "se")

  for (i in 1:nrow(resp))
  {
    temp <- mle(neg_logL,
               start = list(theta = 0),
               fixed = list(resp = resp[i,], alpha = alpha, beta = beta),
               method = "Brent",
               lower = -4,
               upper = 4)

    OUT[i,] <- c(logLik(temp), coef(temp)[1], vcov(temp))
  }
  OUT[,3] <- sqrt(OUT[,3])
  OUT
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

lr_test <-function(resp, model, alpha, beta, ind_theta, col_theta, n_boot = 0){

  odd <- seq(from = 1, to = length(ind_theta), by = 2)
  theta1 <- ind_theta[odd]
  theta2 <- ind_theta[odd+1]
  n_pair <- length(odd)
  n_model <- length(model)
  mod <- matrix(0, nrow = n_pair, ncol = n_model)
  OUT <- vector("list", n_model)
  names(OUT) <- model

  # Helper function for computing P(x > obs)
  pobs <- function(cdf, obs){
    1 - environment(cdf)$y[which.min(abs(environment(cdf)$x-obs))]
  }

  # logL for collaboration models
  for (i in 1:n_model){
    mod[, i] <-
      logL(resp, model = model[i], alpha, beta, theta1, theta2)
  }

  # logL for reference model
  ref <- logL(resp, "twoPL", alpha, beta, col_theta)
  lr <- -2*(mod - ref)

  # Bootstrapping (could fancy this up...)
  if (n_boot > 0){
    theta1_long <- rep(theta1, each = n_boot)
    theta2_long <- rep(theta2, each = n_boot)
    boot_ind <- rep(1:n_pair, each = n_boot)

    for (i in 1:n_model){
      message(cat("Running bootstraps for model", i, "..."),"\r",appendLF=FALSE)
      flush.console()

      boot_data <- sim_data(model[i], alpha, beta, theta1_long, theta2_long)
      boot_mod <- logL(boot_data, model[i], alpha, beta, theta1_long, theta2_long)
      boot_ref <- ml_twoPL(boot_data, alpha, beta)
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
      OUT[[i]] <- temp

    }
  }
  OUT
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

barbell_plot <- function(ind_theta, col_theta, group_score = NULL, legend = "none"){
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
#' Formats responses a response df to match the calibration df
#'
#' Drops items in the target df not in the calibration sample; adds items (with NA entries) in the calibration sample not in the target df. Used for obtaining factor scores from calibration sample \code{ltm}.

#' @param resp the target df to be formatted
#' @param calib the calibration df
#' @param version an optional string used to subset resp via \code{grep(version, names(resp))}
#' @return a version of resp that has the same names as calib
#' @export


format_resp <- function(resp, calib, version = NULL){
  if (!is.null(version)) {
    resp <- resp[grep(version, names(resp))]
  }
  names(resp) <- substr(names(resp), 1, 5)
  resp <- resp[names(resp)%in%names(calib)]
  resp[names(calib)[!names(calib)%in%names(resp)]] <- NA
  resp <- resp[names(calib)]
  resp
}


# EM Functions -------------------------------------------------------------------

EM <- function(models, resp, theta1, theta2, parms, weights = NULL, n_reps = 100, conv = 1e-3) {
  n_models <- length(models)
  p <- rep(1/n_models, n_models)
  l <- likelihood(models, resp, theta1, theta2, parms, weights, Log = F)

  trace <- incomplete_data(l, p)
  i <- 1
  delta <- 1

  while(i <= n_reps & delta > conv) {
    post <- posterior(l, p)
    p <- prior(post)
    trace <- c(trace, incomplete_data(l, p))
    delta <- trace[i+1] - trace[i]
    i <- i + 1
  }
  out <- list(trace, p, post)
  names(out) <- c("trace", "prior", "posterior")
  out
}

likelihood <- function(models, resp, theta1, theta2, parms, weights = NULL, Log = T){
  n_models <- length(models)
  out <- array(0, dim = c(n_models, nrow(resp)))
  if (is.null(weights)) weights <- array(1, dim = dim(resp))

  for (i in 1:n_models) {
    fun <- match.fun(models[i])
    p <- fun(parms$alpha, parms$beta, theta1, theta2)
    out[i,] <- apply(weights * (log(p) * (resp) + log(1-p) * (1-(resp))), 1, sum, na.rm = T)
  }
  if (Log) {out} else {exp(out)}
}

posterior <- function(l, p) {
    temp <- l * p
    t(temp) / apply(temp, 2, sum)
}

prior <- function(post) {
   apply(post, 2, sum) / nrow(post)
}

incomplete_data <- function(l, p) {
  sum(log(apply(l * p, 2, sum)))
}

delta <- function(alpha, beta, theta1, theta2) {
  Min(alpha, beta, theta1, theta2) * (1 - Max(alpha, beta, theta1, theta2))
}

screen <- function(cutoff, alpha, beta, theta1, theta2) {
  screen <- delta(alpha, beta, theta1, theta2)
  screen[screen < cutoff] <- NA
  screen[!is.na(screen)] <- 1
  screen
}

raster_plot <-function(em, sort = F){
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

class_accuracy <- function(em){
  ind <- apply(em$posterior, 1, which.max)
  arr_ind <- cbind(1:nrow(em$posterior), ind)
  cp <- em$posterior[arr_ind]
  tapply(cp, ind, mean)
}

sim_mix <- function(n_obs, n_items, prior = NULL, alpha = NULL, beta = NULL) {
  temp1 <- rnorm(n_obs)
  temp2 <- rnorm(n_obs)
  theta1 <- apply(cbind(temp1, temp2), 1, max)
  theta2 <- apply(cbind(temp1, temp2), 1, min)
  if (is.null(alpha)) { alpha <- rep(1, times = n_items) }
  if (is.null(beta)) { beta <- sort(runif(n_items, -3, 3)) }
  if (is.null(prior)) { prior <- c(.25, .25, .25, .25) }

  out <- matrix(NA, nrow = n_obs, ncol = n_items)
  n <- c(0, round(prior * n_obs))
  models <- c("Ind", "Min", "Max", "AI")

  for(i in 1:length(models)) {
    j <- sum(n[1:i]) + 1
    k <- sum(n[1:(i+1)])
    if (j < k) {
      out[j:k, ] <-  sim_data(models[i], alpha, beta, theta1[j:k], theta2[j:k])
    }
  }
  r <- sample(1:nrow(out))
  parms <- data.frame(beta, alpha)
  EM(models, out[r,], theta1[r], theta2[r], parms)
}
