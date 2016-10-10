devtools::install_github("peterhalpin/BearShare")
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
barbell_plot(ind_theta, col_theta)


parms <- coef(calib_ltm)
beta_C <- parms[grep("C", row.names(parms)), 1]
alpha_C <- parms[grep("C", row.names(parms)), 2]

resp <- col_form[odd, grep("C", names(col_form))]
theta1 <- ind_theta[odd]
theta2 <- ind_theta[odd+1]
alpha <- alpha_C
beta <- beta_C
a <- c(.5, .5)

WAI(alpha_C, beta_C, theta1[1], theta2[1], .5, .5)

neglogWAI(c(.5, .5), resp[1,], alpha_C, beta_C, theta1[1], theta2[1])
gradWAI(c(.5, .5), resp[1,], alpha_C, beta_C, theta1[1], theta2[1])
temp_grad(.5, .5)

WAI <- function(alpha, beta, theta1, theta2, a1, a2){
  w1 <- exp(a1)/(1+exp(a1))
  w2 <- exp(a2)/(1+exp(a2))

  w1 * twoPL(alpha, beta, theta1) + w2 * twoPL(alpha, beta, theta2) + (1-w1 -w2) * twoPL(alpha, beta, theta1) * twoPL(alpha, beta, theta2)
}

neglogWAI <- function(a, resp, alpha, beta, theta1, theta2){
  p <- WAI(alpha, beta, theta1, theta2, a[1], a[2])
  -1*apply(log(p)*(resp) + log(1-p)*(1-(resp)), 1, sum, na.rm = T)
}

gradWAI <- function(a, resp, alpha, beta, theta1, theta2){
	p <- WAI(alpha, beta, theta1, theta2, a[1], a[2])
  r1 <- twoPL(alpha, beta, theta1)
  r2 <- twoPL(alpha, beta, theta2)
	w1 <- exp(a[1])/(1+exp(a[1]))
	w2 <- exp(a[2])/(1+exp(a[2]))

  m <- (resp / p - (1-resp) / (1-p))
  g1 <- w1 * (1-w1) * r1 * (1-r2)
  g2 <- w2 * (1-w2) * r2 * (1-r1)
  -1*c(sum(m*g1, na.rm = T), sum(m*g2, na.rm = T))
}

i <- 1


q <- optim(par = c(1, 1),
      fn = neglogWAI,
			gr = gradWAI,
		  resp = resp[i,],
			alpha = alpha,
			beta = beta,
			theta1 = theta1[i],
			theta2 = theta2[i],
			method = "L-BFGS-B",
			lower = c(-10, -10),
			upper = c(10, 10),
			hessian = T
			)
q
sqrt(diag(solve(q$hessian)))


temp_neglog <- function(a1, a2){
	vectorize(neglogWAI(a, resp, alpha, beta, theta1, theta2){
}

a1 <- a2 <- seq(-10, 10, by = .1)
outer(a1, a2, temp_neglog)


rgb(87, 6, 140, maxColorValue = 255)

persp(theta1, theta2, mix, theta = -25, phi = 15, col = rgb(87, 6, 140, maxColorValue = 255), main = "Min Model", expand = .5, ylab = "")
