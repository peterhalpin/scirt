# Analayses for "Pscychometric Models for Small Group Collaborations"

#devtools::install_github("peterhalpin/BearShare")
library("BearShare")
library("ltm")
library("ggplot2")
library("stats4")

# Stuff from Collab_outcomes.Rmd-----------------------------------------------------------

calib_ltm <- ltm(calibration ~ z1)
beta <- coef(calib_ltm)[,1]
alpha <- coef(calib_ltm)[,2]

ind_form <- col_form <- friends[,-1]
ind_form[, grep("C", names(ind_form))] <- NA
col_form[, grep("I", names(col_form))] <- NA

odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]

ind_theta <- factor.scores(calib_ltm, ind_form, type = "EB", prior = F)$score.dat$z1
col_theta <- factor.scores(calib_ltm, col_, type = "EB", prior = F)$score.dat$z1
barbell_plot(ind_theta, col_theta, legend = "left")

parms <- coef(calib_ltm)
beta_C <- parms[grep("C", row.names(parms)), 1]
alpha_C <- parms[grep("C", row.names(parms)), 2]

# New data -----------------------------------------------------------------------

machine <- "peterfrancis"
machine <- "peterhalpin"
setwd(paste0("/Users/", machine, "/Dropbox/Academic/Projects/CA/Data/response_matrices"))
calib <- read.csv("calibration_2016.csv", check.names = F)
names(calib)

# The whoe thing?
calib_ltm <- ltm(calib[,-1] ~ z1)
parms <- coef(calib_ltm) %>% data.frame
names(parms) <- c("beta", "alpha")
hist(parms$alpha)

collab <- read.csv("collaboration_2016.csv", check.names = F)
names(collab)

format_resp <- function(collab, calib, version){
  temp <- collab[grep(version, names(collab))]
  names(temp) <- substr(names(temp), 1, 5)
  temp[names(calib)[!names(calib)%in%names(temp)]] <- NA
  temp <- temp[names(calib)]
  temp
}

col_form <- format_resp(collab, calib[,-1], "COL")
ind_form <- format_resp(collab, calib[,-1], "IND")

odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]

ss <- apply(ind_form, 1, sum, na.rm = T)
plot(ss[odd], ss[odd+1])

ind_theta <- factor.scores(calib_ltm, ind_form, type = "EB", prior = F)$score.dat$z1
col_theta <- factor.scores(calib_ltm, col_form, type = "EB", prior = F)$score.dat$z1

hist(col_theta)
hist(ind_theta)

# temp fix for shite
ind_theta[abs(ind_theta) > 4] <- NA
col_theta[abs(col_theta) > 4] <- NA

drop_groups <- collab$group_id[is.na(ind_theta) | is.na(col_theta)]
col_form <- col_form[!collab$group_id%in%drop_groups, ]
ind_theta <- ind_theta[!collab$group_id%in%drop_groups]
col_theta <- col_theta[!collab$group_id%in%drop_groups]
odd <- seq(1, nrow(col_form), by = 2)

barbell_plot(ind_theta[1:60], col_theta[1:60])


# EM -------------------------------------------------------------------
models <- c("Ind", "Min", "Max", "AI")
theta1 <- ind_theta[odd]
theta2 <- ind_theta[odd+1]

cor(theta1, theta2)

weights <- delta(parms$alpha, parms$beta, theta1, theta2)
s <- screen(.05, parms$alpha, parms$beta, theta1, theta2)


sanity <- EM(models, ind_form[odd,]*ind_form[odd+1,]*s, theta1, theta2, parms)
sanity$prior

col <- EM(models, col_form[odd,]*s, theta1, theta2, parms)
col$prior

temp <- col$posterior

temp_ind <- which(temp%*% 1:4 > 3 & temp%*% 1:4 < 3.5)
temp_ind <- sort(c(temp_ind*2, (temp_ind*2)-1))

barbell_plot(ind_theta[temp_ind], col_theta[temp_ind])
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















# Functions for WAI -----------------------------------------------------------------

WAI <- function(w1, w2, alpha, beta, theta1, theta2){
  w1 * twoPL(alpha, beta, theta1) + w2 * twoPL(alpha, beta, theta2) + (1 - w1 - w2) * twoPL(alpha, beta, theta1) * twoPL(alpha, beta, theta2)
}


# Plot WAI in weightspace -------------------------------------------------------------------

w1 <- w2 <- seq(0, 1, by = .05)
alpha <- 1
beta <- 0
theta1 <- .5
theta2 <- -.5
R1 <- outer(w1, w2, WAI, alpha, beta, theta1, theta2)

par(mfrow = c(1, 2))

persp(w1, w2, R1, theta = -55, phi = 20, main = " ", col = "grey", zlim = c(0, 1), expand = .5, ticktype = "detailed", nticks = 2, shade = .25)

# colors?
col <- c("grey", "black")
ncr <- ncol(R)
nrr <- nrow(R)
Rfacet <- (R[-1, -1] + R[-1, -ncr] + R[-nrr, -1] + R[-nrr, -ncr])/4
cutpoints <- c(0, twoPL(1, 0, theta1), 1)
facetcol <- cut(Rfacet, cutpoints)
fill[R > as.numeric(twoPL(1, 0, theta1))] <- "blue"

# Plot WAI in thetaspace -------------------------------------------------------------------

theta1 <- theta2 <- seq(-4, 4, by = .25)
alpha <- 1
beta <- 0
w1 <- .5
w2 <- .7
WAI2 <- function(theta1, theta2, alpha, beta, w1, w2) {
  WAI(w1, w2, alpha, beta, theta1, theta2)
}

R2 <- outer(theta1, theta2, WAI2, alpha, beta, w1, w2)
persp(theta1, theta2, R2, theta = -55, phi = 20, main = " ", col = "grey", expand = .5, ticktype = "detailed", nticks = 2)


# plots IRFS
theta <- seq(-4, 4, by = .25)
nIRF <- 4
a <- 2
alpha <- c(a, a, a, a/3)
beta <- c(-1, 1, 0, 0)
temp <- cbind(rep(1:nIRF, each = length(theta)), rep(theta, times = length(nIRF)))
q <- c()
for(i in 1:nIRF) q <- c(q, twoPL(alpha[i], beta[i], theta))
temp <- data.frame(cbind(temp, q))
names(temp) <- c("item", "theta", "p")
head(temp)
temp$item <- as.factor(temp$item)
ggplot(temp, aes(x = theta, y = p, group = item)) + geom_line(aes(col = item)) + geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 3) + annotate("text", x = .5, y = 1.02, label = "theta_1") +  annotate("text", x = -1.5, y = 1.02, label = "theta_2")


twoPL(alpha, beta, -1) * (1-twoPL(alpha, beta, 1))
