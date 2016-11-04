# Last update: 04/12/2016
# Functions for computing bootstrapped LR tests of collaboration models.
# User beware: functions not written to check or handle input errors.

require(stats4)
require(ggplot2)

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
#' @return A list of lenght \code{length(model)}, each element of which is a data frame with \code{nrow(resp)} rwos containing the output for lr_tests for each pair.
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
#' @param passed to ggplot \code{legend.position}
#' @return A barbell plot
#' @export

barbell_plot <- function(ind_theta, col_theta, legend = "none"){

  data <- data.frame(ind_theta, col_theta)
  lim <- c(min(data)-.2, max(data)+.2)
  data$pairs <- factor(rep(1:(length(ind_theta)/2), each = 2))

  ggplot(data = data, aes(x = ind_theta, y = col_theta, group = pairs)) +
    geom_line(aes(color = pairs)) +
    geom_point(aes(color = pairs), size = 4) +
    scale_x_continuous(limits = lim) +
    scale_y_continuous(limits = lim) +
    geom_abline(intercept = 0, slope = 1, col = "grey") +
    theme(legend.position = legend) +
    ggtitle("Collaborative vs Individual Performance") +
    xlab("Individual Theta")+
    ylab("Collaborative Theta")+
    theme(axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13)
    )
}

























#Scraps --------------------------------------------------------------------------


#' #' Integrand of likelihood of a single response pattern and stated model
#' #'
#' #' Function is passed to \code{likelihood} to evalute likelihood of a model; see \code{likelihood} for details of args.
#' #'
#' #' @return likelihood*prior
#' #'
#' f_cubature <- function(x, resp, model, alpha, beta, theta, se){
#'   if (model == "twoPL"){
#'     irf <-  twoPL(alpha, beta, x)
#'     prior <- dnorm(x, theta, se)
#'   } else {
#'     fun <- match.fun(model)
#'     irf <- fun(alpha, beta, x[1], x[2])
#'     prior <- dnorm(x[1], theta[1], se[1])*dnorm(x[2], theta[2], se[2])
#'   }
#'   lik <- exp(apply(log(irf)*(resp) + log(1-irf)*(1-(resp)), 1, sum, na.rm = T))
#'   lik*prior
#' }
#'
#' #--------------------------------------------------------------------------
#' #' Compute the likelihood of a single response pattern for a stated model
#' #'
#' #' This function is a wrapper that passes \code{f_cubature} to \code{adaptIntegrate}. Requires package \code{cubature}.
#' #'
#' #' @param resp is the response pattern whose likelihood is desired
#' #' @param model is one of \code{c("twoPL", "Ind", "Min", "Max", "AI") }
#' #' @param alpha the item discriminations
#' #' @param beta the item difficulties
#' #' @param theta is a vector of length 1 or 2 depending on the model
#' #' @param se is the standard error(s) for theta (required!)
#' #' @param lim is \code{c(lowerLimit, upperLimit)} passed to adaptIntegrate
#' #' @returns the output from evaluating adpatIntegrate on f_cubature with the stated args (a named list).
#'
#' likelihood <- function(resp, model, alpha, beta, theta, se){
#'   n_theta <- length(theta)
#'   lim <- c(min(qnorm(.01, theta, se)), max(qnorm(.99, theta, se)))
#'
#'   adaptIntegrate(f_cubature,
#'                  lowerLimit = rep(lim[1], n_theta),
#'                  upperLimit = rep(lim[2], n_theta),
#'                  resp = resp,
#'                  model = model,
#'                  alpha = alpha,
#'                  beta = beta,
#'                  theta = theta,
#'                  se = se,
#'                  maxEval = 500
#'                  )
#' }
#'
#'

#' #--------------------------------------------------------------------------
#' #' Computes likelihood ratio test for stated model(s) of pairwise collaboration.
#' #'


#' #' @param data a matrix of binary response patterns
#' #' @param model is one or more of \code{c("Ind", "Min", "Max", "AI") }
#' #' @param alpha the item discriminations of the resp items
#' #' @param beta the item difficulties of the resp items
#' #' @param ind_theta the latent trait for each member, as estimated from a non-collaborative form
#' #' @param col_theta the latent trait for each member, as estimated from a collaborative form
#' #' @return A \code{length(ind_theta)/2} by \code{length(model)} matrix, each column of which contains the lr_tests for each model.
#'
#' lr_test <-function(data, model, alpha, beta, ind_theta, ind_se, col_theta, col_se){
#'   odd <- seq(from = 1, to = length(ind_theta), by = 2)
#'   even <- odd +1
#'   n_pair <- length(odd)
#'   n_model <- length(model)
#'   mod <- matrix(0, nrow = n_pair, ncol = n_model)
#'   ref <- mod[,1]
#'   # likelihoods for each model and each resp pattern -- yick!
#'   for (i in 1:n_pair){
#'
#'     for(j in 1:n_model){
#'       mod[i, j] <- likelihood(data[i,],
#'                     model = "Min",
#'                     alpha,
#'                     beta,
#'                     theta = ind_theta[c(odd[i], even[i])],
#'                     se = ind_se[c(odd[i], even[i])]
#'                     )$integral
#'     }
#'
#'     ref[i] <- likelihood(data[i,],
#'                 model = "twoPL",
#'                 alpha,
#'                 beta,
#'                 theta = col_theta[i],
#'                 se = col_se[i]
#'                 )$integral
#'
#'   }
#'
#'   OUT <- -2*(log(mod) - log(ref))
#'   colnames(OUT) <- model
#'   OUT
#' }
