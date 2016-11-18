
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

# Calibration set up -------------------------------------------------------------------
calib_ltm <- ltm(calibration ~ z1)
parms <- coef(calib_ltm) %>% data.frame
names(parms) <- c("beta", "alpha")

# Friends or Pilot data ------------------------------------------------
data <- friends[,-1]

data <- pilot[,-1]

data <- read.csv("~/Dropbox/Academic/Projects/CA/Data/Response_Matrices/collaboration_sample.csv")
  data <- data[,-c(1:3)]
  dim(data)
  head(data)
# Set up forms and thetas ------------------------------------------------

ind_form <- col_form <- data
ind_form[, grep("C", names(ind_form))] <- NA
col_form[, grep("I", names(col_form))] <- NA

odd <- seq(1, nrow(data), by = 2)
col_resp <- col_form[odd,]*col_form[odd+1,]
# sanity <- ind_form[odd,]*ind_form[odd+1,]

theta <- factor.scores(calib_ltm, ind_form, type = "EAP", prior = F)$score.dat$z1
theta[51] <- NA #wtf
theta1 <- theta[odd]
theta2 <- theta[odd+1]

# EM -------------------------------------------------------------------
models <- c("Ind", "Min", "Max", "AI")
weights <- delta(parms$alpha, parms$beta, theta1, theta2)
s <- screen(.05, parms$alpha, parms$beta, theta1, theta2)

#S <- screen(.025, parms$alpha, parms$beta, theta1, theta2)
col1 <- EM(models, col_resp, theta1, theta2, parms)
col1$prior

# Plots -------------------------------------------------------------------
l <- likelihood(models, col_resp*s, theta1, theta2, parms, Log = F)
  temp <- posterior(l, col1$prior)

temp <- col1$posterior


table(apply(temp, 1, which.max))
hist(temp%*% 1:4, breaks = 20)

temp <- data.frame(cbind(1:nrow(temp), temp))
names(temp) <- c("pair", "Ind", "Min", "Max", "AI")
head(temp)

q <- reshape(temp,
  varying = names(temp)[-1],
  v.names = "prob",
  timevar = "model",
  times = names(temp)[-1],
  direction = "long"
  )

q$model <- ordered(q$model, c("Ind", "Min", "Max", "AI"))
ggplot(q, aes(pair, model, fill = prob)) + geom_raster()

# barbells study -------------------------------------------------------------------

col_resp <- col_form[odd, ]*col_form[odd+1,]
col_theta <- factor.scores(calib_ltm, col_resp, type = "EAP", prior = F)$score.dat$z1
col_theta[col_theta > 18] <- NA


q <- cbind(theta, rep(col_theta, each = 2), rep(temp%*%1:4, each = 2))
q <- data.frame(q)
q$delta <- q[,2] - q[,1]
plot(q[,3], q$deltaA)
ind1 <- q[,3] < 2
ind2 <- q[,3] > 2 & q[,3] < 3
ind3 <- q[,3] > 3.2

barbell_plot(q[ind1,1], q[ind1,2])
barbell_plot(q[ind2,1], q[ind2,2])
barbell_plot(q[ind3,1], q[ind3,2])


length(theta2)



# Simuation study -------------------------------------------------------------------
# Conclusions - dont use the weights?

n_obs <- 1000
n_items <- 30
temp1 <- rnorm(n_obs)
temp2 <- rnorm(n_obs)
theta1 <- apply(cbind(temp1, temp2), 1, max)
theta2 <- apply(cbind(temp1, temp2), 1, min)
beta <- sort(runif(n_items, -3, 3))
alpha <- rep(1, time = n_items)
out <- matrix(NA, nrow = n_obs, ncol = n_items)

prop <- c(.1, .3, .5, .1)
n <- c(0, round(prop * n_obs))
models <- c("Ind", "Min", "Max", "AI")

for(i in 1:length(models)) {
  j <- sum(n[1:i]) + 1
  k <- sum(n[1:(i+1)])
  if (j < k) {
    out[j:k, ] <-  sim_data(models[i], alpha, beta, theta1[j:k], theta2[j:k])
  }
}

parms <- data.frame(beta, alpha)
weights <- delta(parms$alpha, parms$beta, theta1, theta2) / apply(weights, 1, sum)

s <- screen(.05, parms$alpha, parms$beta, theta1, theta2)

col1 <- EM(models, out, theta1, theta2, parms)
col2 <- EM(models, out, theta1, theta2, parms, weights)
col3 <- EM(models, out*s, theta1, theta2, parms)

col1$prior; col2$prior; col3$prior

par(mfrow = c(1, 3))
  plot(col1$posterior[,3])
  points(col1$posterior[,4], col = 2)

  plot(col2$posterior[,3])
  points(col2$posterior[,4], col = 2)

  plot(col3$posterior[,3])
  points(col3$posterior[,4], col = 2)


table(apply(col1$posterior, 1, which.max))
table(apply(col2$posterior, 1, which.max))
table(apply(col3$posterior, 1, which.max))

col2$prior

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
