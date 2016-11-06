
# Bear e.g. ---------------------------------------------------------------------------

#devtools::install_github("peterhalpin/BearShare")
library("BearShare")
library("ltm")
library("ggplot2")
library("stats4")

# for all responsdents and each model
# a) compute likelihood
# b) compute posterior
# c) compute prior (using exponential parameterization)
# iterate b and c; beware of underflow

# data set up -------------------------------------------------------------------
calib_ltm <- ltm(calibration ~ z1)
parms <- coef(calib_ltm) %>% data.frame
names(parms) <- c("beta", "alpha")

ind_form <- col_form <- friends[,-1]
ind_form[, grep("C", names(ind_form))] <- NA
col_form[, grep("I", names(col_form))] <- NA

odd <- seq(1, nrow(col_form), by = 2)
col_resp <- col_form[odd,]*col_form[odd+1,]
ind_resp <- ind_form[odd,]*ind_form[odd+1,]

theta <- factor.scores(calib_ltm, ind_form, type = "EB", prior = F)$score.dat$z1
theta1 <- theta[odd]
theta2 <- theta[odd+1]

# EM -------------------------------------------------------------------

models <- c("Ind", "Min", "Max", "AI")
weights <- delta(parms$alpha, parms$beta, theta1, theta2) / .25
#S <- screen(.025, parms$alpha, parms$beta, theta1, theta2)

ind1 <- EM(models, ind_resp, theta1, theta2, parms, weights)
ind2 <- EM(models, ind_resp, theta1, theta2, parms)
ind2

col1 <- EM(models, col_resp, theta1, theta2, parms, weights)
col2 <- EM(models, col_resp, theta1, theta2, parms)
col1

EM <- function(models, resp, theta1, theta2, parms, weights = NULL, n_reps = 100, conv = 1e-5) {
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
    out <- t(temp) / apply(temp, 2, sum)
}

prior <- function(post) {
   apply(post, 2, sum) / nrow(post)
}

incomplete_data <- function(l, p) {
  sum(log(apply(l * p, 2, sum)))
}

delta <- function(alpha, beta, theta1, theta2){
  Min(alpha, beta, theta1, theta2) * (1 - Max(alpha, beta, theta1, theta2))
}

screen <- function(cutoff, alpha, beta, theta1, theta2){
  screen <- delta(alpha, beta, theta1, theta2)
  screen[screen < cutoff] <- NA
  screen[!is.na(screen)] <- 1
  screen
}
