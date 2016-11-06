
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

parms <- coef(calib_ltm)
beta_C <- parms[grep("C", row.names(parms)), 1]
alpha_C <- parms[grep("C", row.names(parms)), 2]

# Set up arges for screen -----------------------------------------------------------------------

delta <- function(alpha, beta, theta1, theta2){
  Min(alpha, beta, theta1, theta2) * (1 - Max(alpha, beta, theta1, theta2))
}

screen <- function(cutoff, alpha, beta, theta1, theta2){
  screen <- delta(alpha, beta, theta1, theta2)
  screen[screen < cutoff] <- NA
  screen[!is.na(screen)] <- 1
  screen
}

resp <- col_form[odd,grep("C", names(col_form))]
models <- c("Ind", "Min", "Max", "AI")
d <- delta(alpha_C, beta_C, ind_theta[odd], ind_theta[odd+1])
S <- screen(.05, alpha_C, beta_C, ind_theta[odd], ind_theta[odd+1])

lr1 <- lr_test(resp, models, alpha_C, beta_C, ind_theta, col_theta[odd], n_boot = 500)
lr2 <- lr_test(resp*S, models, alpha_C, beta_C, ind_theta, col_theta[odd], n_boot = 500)
i = 4
p1 <- p2 <- matrix(0, nrow = 7, ncol = 4)
p1[] <- unlist(lapply(lr1, function(x) x[,1]))
p2[] <- unlist(lapply(lr2, function(x) x[,1]))
p1;p2
apply(S, 1, sum, na.rm = T)
R

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

WAI <- function(w, alpha, beta, theta1, theta2){
  w[1] * twoPL(alpha, beta, theta1) + w[2] * twoPL(alpha, beta, theta2) + (1 - w[1] - w[2]) * twoPL(alpha, beta, theta1) * twoPL(alpha, beta, theta2)
}

neglogGAI <- function(w, resp, alpha, beta, theta1, theta2){
  #w <- exp(a)/(1+exp(a))
  p <- GAI(w, alpha, beta, theta1, theta2)
  p[p>1] <- .999
  - apply(log(p) * (resp) + log(1-p) * (1-(resp)), 1, sum, na.rm = T)
}

# sum to 1
lambda1 <- function(w, resp, alpha, beta, theta1, theta2){
  #sum(exp(a)/(1+exp(a)))
  sum(w)
}

# between 1 and zero
lambda2 <- function(w, resp, alpha, beta, theta1, theta2){
  #exp(a)/(1+exp(a))
  w
}

# ML for GAI ------------------------------------------------------------------

x0 <- c(.25,.25,.25,.25)
i = 1

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
P <- round(q$par, 4)
sqrt(P*(1-P))

ind <- which(P > .0001 & P < .9999)
SE <- (solve(q$hessian[ind+4, ind+4]))
P; sqrt(diag(SE))

# Plotting SE ------------------------------------------------------------------

SE <- function(w1, w2, alpha, beta, theta1, theta2){
  R <- WAI(c(w1, w2), alpha, beta, theta1, theta2)
  m <-  1 / sqrt(R * (1-R))
  P1 <- twoPL(alpha, beta, theta1)
  P2 <- twoPL(alpha, beta, theta2)
  Q1 <- 1 - P1
  Q2 <- 1 - P2
  U <- sum((m * P1 * Q2)^2)
  V <- sum((m * P2 * Q1)^2)
  W <- sum(m^2 * P1 * P2 * Q1 * Q2)^2
  sqrt(V / (U * V - W))
}

plotSE <- function(w1, w2, alpha, beta, theta1, theta2){
   temp <- function(w1, w2){
      SE(w1, w2, alpha, beta, theta1, theta2)
   }
   mapply(temp, w1, w2)
}

plotSE2 <- function(theta1, theta2, alpha, beta, w1, w2){
   temp <- function(theta1, theta2){
      SE(w1, w2, alpha, beta, theta1, theta2)
   }
   mapply(temp, theta1, theta2)
}

w1 <- w2 <- seq(0, 1, by = .02)
l <- outer(w1, w2, plotSE, alpha, beta, -1, 1)
persp(w1, w2, l, theta = -75, phi = 20, main = "SE", expand = .5, ticktype = "detailed", nticks = 5)

theta1 <- theta2 <- seq(-2, 2, by = .1)
theta1
m <- outer(theta1, theta2, plotSE2, alpha, beta, 1, 1)
m[m>10] <- NA


persp(theta1, theta2, m, theta = -45, phi = 20, main = "SE", expand = .5, ticktype = "detailed", nticks = 5)

hist(l, breaks = 200)




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

WAI <- function(w, alpha, beta, theta1, theta2){
  #w <- exp(a)/(1+exp(a))
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


# More SEs -----------------------------------------------------------------------
alpha = 1
beta = 0
theta1 = 1
theta2 = 0
x = 1
w <- seq(0, 1, by = .01)
info_w <- function(w, alpha, beta, theta1, theta2){
  R <- lapply(w, function (x) WAI(c(x, 0), alpha, beta, theta1, theta2)) %>% unlist
  P1 <- twoPL(alpha, beta, theta1)
  P2 <- twoPL(alpha, beta, theta2)
  1 / sqrt((P1*(1-P2))^2 / (1-R) / R)
}

plot(w, info_w(w, 1, 0, 0, 0))

#  -----------------------------------------------------------------------

alpha <- seq(0, 3, by = .01)
beta = 0
theta1 = 0
theta2 = 0
x = 1
w = 1
info_alpha <- function(alpha, w, beta, theta1, theta2){
  R <-lapply(alpha, function (x) WAI(c(w, 0), x, beta, theta1, theta2)) %>% unlist
  P1 <- twoPL(alpha, beta, theta1)
  P2 <- twoPL(alpha, beta, theta2)
  1/ sqrt(1/(P1*(1-P2))^2 / (1-R) / R)
}

plot(alpha, info_alpha(alpha, 1, .5, 0, 1))


theta1 <- theta2 <- seq(-3, 3, by = .1)
beta = 0
alpha = 1
x = 1
w = 1

info_theta <- function(theta1, theta2, w, alpha, beta){
  R <- mapply(function(x, y) WAI(c(w, 0), alpha, beta, x, y), theta1, theta2) %>% unlist
  P1 <- twoPL(alpha, beta, theta1)
  P2 <- twoPL(alpha, beta, theta2)
  1 / sqrt((P1*(1-P2))^2 / (1-R) / R)
}

l <- outer(theta1, theta2, info_theta, w, alpha, beta)
l[l > 10] = NA
persp(theta1, theta2, l, theta = -25, phi = 20, ticktype = "detailed", nticks = 5)
