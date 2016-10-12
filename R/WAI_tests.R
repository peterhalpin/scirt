# Bear e.g. ---------------------------------------------------------------------------

#devtools::install_github("peterhalpin/BearShare")
library("BearShare")
library("ltm")
library("ggplot2")
library("stats4")


calib_ltm <- ltm(calibration ~ z1)
beta <- coef(calib_ltm)[,1]
alpha <- coef(calib_ltm)[,2]

ind_form <- col_form <- friends[,-1]
ind_form[, grep("C", names(ind_form))] <- NA
col_form[, grep("I", names(col_form))] <- NA

odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]

ind_theta <- factor.scores(calib_ltm, ind_form, type = "EB", prior = F)$score.dat$z1
col_theta <- factor.scores(calib_ltm, col_form, type = "EB", prior = F)$score.dat$z1
barbell_plot(ind_theta, col_theta, legend = "left")
barbell_plot

parms <- coef(calib_ltm)
beta_C <- parms[grep("C", row.names(parms)), 1]
alpha_C <- parms[grep("C", row.names(parms)), 2]

# Set up arges for GAI ---------------------------------------------------------------------------

resp <- col_form[odd, grep("C", names(col_form))]
theta1 <- ind_theta[odd]
theta2 <- ind_theta[odd+1]
alpha <- alpha_C
beta <- beta_C
a <- c(.25, .25, .25, .25)

# Functions for GAI -----------------------------------------------------------------------------
#install.packages("Rsolnp")
library("Rsolnp")

GAI <- function(w, alpha, beta, theta1, theta2){
  (w[1] + w[2]) * twoPL(alpha, beta, theta1) + (w[1] + w[3]) * twoPL(alpha, beta, theta2) + (w[4] - w[1]) * twoPL(alpha, beta, theta1) * twoPL(alpha, beta, theta2)
}

neglogGAI <- function(a, resp, alpha, beta, theta1, theta2){
  w <- exp(a)/(1+exp(a))
  p <- GAI(w, alpha, beta, theta1, theta2)
  p[p>1] <- .999
  - apply(log(p) * (resp) + log(1-p) * (1-(resp)), 1, sum, na.rm = T)
}

# sum to 1
lambda1 <- function(a, resp, alpha, beta, theta1, theta2){
  sum(exp(a)/(1+exp(a)))
}

# between 1 and zero
lambda2 <- function(a, resp, alpha, beta, theta1, theta2){
  exp(a)/(1+exp(a))
}

# ML for GAI ------------------------------------------------------------------

x0 <- c(0,0,0,0)
i = 7

q <- solnp(x0, neglogGAI,
   eqfun = lambda1,
   eqB = 1,
   ineqfun = lambda2,
   ineqUB = rep(.9999, 4),
   ineqLB = rep(.0001, 4),
   resp = resp[i,],
   alpha = alpha,
   beta = beta,
   theta1 = theta1[i],
   theta2 = theta2[i])
q

temp <- round(exp(q$pars)/(1+exp(q$pars)),3)
ind <- which(temp > 0 & temp < 1) + 4
temp
sqrt(solve(q$hessian[ind, ind]))

# WAI  ------------------------------------------------------------------

gradGAI <- function(w, resp, alpha, beta, theta1, theta2){
  p <- GAI(w, alpha, beta, theta1, theta2)
  m <- (resp / p - (1-resp) / (1-p))
  p1 <- m * Ind(alpha, beta, theta1, theta2)
  p2 <- m * Min(alpha, beta, theta1, theta2)
  p3 <- m * Max(alpha, beta, theta1, theta2)
  p4 <- m * AI(alpha, beta, theta1, theta2)
  1 * c(sum(p1, na.rm = T), sum(p2, na.rm = T), sum(p3, na.rm = T), sum(p4, na.rm = T))
}

# Plots for WAI

WAI <- function(a, alpha, beta, theta1, theta2){
  w <- exp(a)/(1+exp(a))
  w[1] * twoPL(alpha, beta, theta1) + w[2] * twoPL(alpha, beta, theta2) + (1 - w[1] - w[2]) * twoPL(alpha, beta, theta1) * twoPL(alpha, beta, theta2)
}

neglogWAI <- function(a, resp, alpha, beta, theta1, theta2){
  p <- WAI(a, alpha, beta, theta1, theta2)
  1 * apply(log(p) * (resp) + log(1-p) * (1-(resp)), 1, sum, na.rm = T)
}

# Why does this take forever??
plot_neglogWAI <- function(a1, a2, resp, alpha, beta, theta1, theta2){
  temp <- function(w1, w2){
      neglogWAI(c(w1, w2), resp, alpha, beta, theta1, theta2)
    }
  mapply(temp, a1, a2)
}

gradWAI <- function(a, resp, alpha, beta, theta1, theta2){
	p <- WAI(a, alpha, beta, theta1, theta2)
  r1 <- twoPL(alpha, beta, theta1)
  r2 <- twoPL(alpha, beta, theta2)
	w1 <- exp(a[1])/(1+exp(a[1]))
	w2 <- exp(a[2])/(1+exp(a[2]))

  m <- (resp / p - (1-resp) / (1-p))
  g1 <- w1 * (1-w1) * r1 * (1-r2)
  g2 <- w2 * (1-w2) * r2 * (1-r1)
  -1 * c(sum(m*g1, na.rm = T), sum(m*g2, na.rm = T))
}

# Plot WAI likelihood -----------------------------------------------------------------------

temp_log(a1, a2)
a1 <- a2 <- seq(-10, 10, by = .5)

plot_neglogWAI(a1, a2, resp[i,], alpha, beta, theta1[i], theta2[i])

l <- outer(a1, a2, plot_neglogWAI, resp[i,], alpha, beta, theta1[i], theta2[i])
persp(a1, a2, l, theta = 55, phi = 20, main = "Min Model", expand = .5, ylab = "a2")

hist(l, breaks = 200)
