# Last update: 04/12/2016
# Functions for computing bootstrapped LR tests of collaboration models.
# User beware: functiosn not written to check or handle input errors.


#--------------------------------------------------------------------------
#' Item response function for 2PL
#'
#' Computes a matrix of probabilities for correct responses under 2PL model
#'

#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta the latent trait
#' @return \code{length(theta)} by \code{length(alpha)} matrix of response probabilities

twoPL <-function(alpha, beta, theta){
  Z <- matrix(0, nrow = length(theta), ncol = length(alpha))
  Z <- Z + theta
  Z <- t(alpha*(t(Z) - beta))
  1/(1+exp(-Z))
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
#' @return \code{length(theta)} by \code{length(alpha)} matrix of response probabilities

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
#' @return \code{length(theta)} by \code{length(alpha)} matrix of response probabilities

Min <- function(alpha, beta, theta1, theta2){
  theta <- apply(cbind(theta1, theta2), 1, min, na.rm = T)
  twoPL(alpha, beta, theta)
}


#--------------------------------------------------------------------------
#' Item response function for ``the Maximum model"
#'
#' Computes a matrix of probabilities for correct responses using the maximum model of collaboration for the members and the 2PL model for the items
#'
#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return \code{length(theta)} by \code{length(alpha)} matrix of response probabilities

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

#--------------------------------------------------------------------------
AI <- function(alpha, beta, theta1, theta2){
  twoPL(alpha, beta, theta1) + twoPL(alpha, beta, theta2) -  twoPL(alpha, beta, theta1) * twoPL(alpha, beta, theta2)
}


#--------------------------------------------------------------------------
#' Log-likelihood of given matrix of binary responses for a stated model
#'
#' Computes a vector of log-likelihoods for a given matrix binary response patterns.
#'
#' @param resp the matrix binary responses
#' @param model is one of \code{c("twoPL", "Ind", "Min", "Max", "AI") }
#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta1 the latent trait for member 1
#' @param theta2 the latent trait for member 2
#' @return An \code{nrow(resp)}-vector of log-likleihoods for each response pattern.

logL <- function(resp, model, alpha, beta, theta1, theta2 = NULL){
  if(model == "twoPL"){
    p <- twoPL(alpha, beta, theta1)
  }else{
    fun <- match.fun(model)
    p <- fun(alpha, beta, theta1, theta2)
  }
   apply(log(p)*(resp) + log(1-p)*(1-(resp)), 1, sum, na.rm = T)
}


#--------------------------------------------------------------------------
#' Simulate data from 2PL model.
#'
#' Simulate data using the 2PL model.
#' @param alpha the item discriminations
#' @param beta the item difficulties
#' @param theta the values of latent trait to simulate data for
#' @return \code{length(theta)} by \code{length(alpha)} matrix of binary response patterns

sim_twoPL <- function(alpha, beta, theta = 0){
  n_row <- length(theta)
  n_col <- length(alpha)
  Q <- array(runif(n_row * n_col), dim = c(n_row, n_col))
	P <- twoPL(alpha, beta, theta)
	OUT <- ifelse(P > Q, 1, 0)
	colnames(OUT) <- names(alpha)
	OUT
}


#--------------------------------------------------------------------------
#' Likelihood ratio tests for various models of collaboration.
#'
#' Computes a boostrapped likelihood ratio test for one or more models of pairwise collaboration, given ``assumed to known" item and person parameters (i.e., neither estimation error in item parameters nor prediction error in latent variables is accounted for by this procedure).
#'
#' @param resp the matrix binary data from the \strong{collaborative responses}
#' @param model is one or more of \code{c("twoPL", "Ind", "Min", "Max", "AI") }
#' @param alpha the item discriminations of the resp items
#' @param beta the item difficulties of the resp items
#' @param ind_theta the latent trait for each member, as estimated from a non-collaborative form
#' @param col_theta the latent trait for each member, as estimated from a collaborative form
#' @return A \code{length(ind_theta)/2} by \code{length(model)} matrix, each column of which contains the lr_tests for each model.
#'

lr_test <-function(resp, model, alpha, beta, ind_theta, col_theta){

  odd <- seq(from = 1, to = length(ind_theta), by = 2)
	theta1 <- ind_theta[odd]
	theta2 <- ind_theta[odd+1]
	n_pairs <- length(odd)
	n_model <- length(model)
  model_logL <- data.frame(matrix(0, nrow = n_pairs, ncol = n_model))

	for(i in 1:n_model){
   model_logL[, i] <-
     logL(resp, model = model[i], alpha, beta, theta1, theta2)
  }

  ref_logL <- logL(resp, "twoPL", alpha, beta, col_theta[odd])

  out <- -2*(model_logL - ref_logL)
  colnames(out) <- model
  out
 }


