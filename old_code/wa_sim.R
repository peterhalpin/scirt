# devtools::install_github("peterhalpin/cirt")
# library("cirt")

library("ggplot2")
library("dplyr")
source("~/github/cirt/R/cIRF_functions.R")
source("~/github/cirt/R/IRF_functions.R")
source("~/github/cirt/R/bootstrapping.R")
source("~/github/cirt/old_code/functions_WAI2.R")
# ------------------------------------------------------------
# Data simulation
# Parameter recovery with and without PV
# Var increase as a function of item information
# Test Redundancy?
# ------------------------------------------------------------

# Constants
set.seed(101)

# Data generating parameters
n_obs <- 100
n_items <- 25
n_reps <- 100
alpha <- runif(n_items, .65, 2.5)
beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
parms <- data.frame(alpha, beta)
theta <- rnorm(n_obs*2)
theta_se <- SE(parms[sample.int(n_items, 25), ], theta)
a <- rnorm(n_obs)
w <- exp(a)/(1+exp(a))

# Other stuff
item_names <- paste0("item", 1:n_items)
row.names(parms) <- item_names
odd <- seq(1, n_obs*2, by = 2)
rep_order <- order(rep(1:n_obs, times = n_reps))
test_length <- 25

# # ------------------------------------------------------------
# # Simulated example A: Similar abililty
# # ------------------------------------------------------------
# ind <- order(theta)
# a1 <- ind[odd]
# a2 <- ind[odd + 1]
#
# A <- data_gen(1, w, parms, theta[odd], theta[odd+1])
# ml_a <- ML_WA(A[item_names], parms, A$theta1, A$theta2)
#
# plot(w, ml_a$w)
# plot(ml_a$w, ml_a$se)
#
# pv_a <- pv_gen(n_reps, A[item_names], parms, theta[a1], theta[a2], theta_se[a1], theta_se[a2], model = A$model)
#
#
# # ------------------------------------------------------------
# # Simulated example B: Disparate ability
# # ------------------------------------------------------------
# b1 <- ind[1:n_obs]
# b2 <- ind[(n_obs + 1):(n_obs * 2)]
# B <- data_gen(1, mix_prop, parms, theta[b1], theta[b2])
# em_b <- EM(models, B[item_names], parms, theta[b1], theta[b2])
#
# pv_b <- pv_gen(n_reps, B[item_names], parms, theta[b1], theta[b2], theta_se[b1], theta_se[b2], model = B$model)
# em_pv_b <- parallel::mclapply(1:n_reps, em_parallel, sim_data = pv_b, parms = parms)

# ------------------------------------------------------------
# Simulated example C: Randomly selected partners
# w = {0,1} is problematic
#
# ------------------------------------------------------------

theta1 <- theta[odd]
theta2 <- theta[odd + 1]
data <- data_gen(1, w, parms, theta1, theta2)
resp <- data[item_names]
ml <- mle_WA(resp, parms, data$theta1, data$theta2, SE = "exp", starts = w)
plot(ml$w, ml$se)

pv_data <- pv_gen(n_reps, resp, parms, theta1, theta2, theta_se[odd], theta_se[odd+1], weights = w)
ml_pv <- mle_WA(pv_data[item_names], parms, pv_data$theta1, pv_data$theta2, starts = pv_data$w, parallel = T)

pv <- cbind(pv_data, ml_pv)
head(pv)

# ------------------------------------------------------------
# Parameter recovery
# ------------------------------------------------------------
pv_w <- tapply(pv$w, pv$pairs, mean)
plot(ml$w, pv_w)
var_w <- tapply(pv$se^2, pv$pairs, mean)
var_b <- tapply(pv$w, pv$pairs, var)
pv_se <- sqrt(var_w + (1 + 1/n_reps) * var_b)

ind <- var_w!= 0 & var_b!= 0
var_increase <- (1 + 1/n_reps) * var_b/var_w
var_prop <- (1 + 1/n_reps) * var_b/pv_se^2

summary(var_prop)
summary(var_increase)



# ------------------------------------------------------------
# Goodness of fit
# ------------------------------------------------------------

logL <- ml$logL
temp1 <- data_gen(n_reps, pv_w, parms, data$theta1, data$theta2)
temp2 <- mle_WA(temp1[item_names], parms, temp1$theta1, temp1$theta2, SE = "exp", starts = temp1$w)
sim <- cbind(temp1, temp2)
head(sim)

quant <- tapply(sim$logL, sim$pairs, function(x) quantile(x, p = c(.5, .025))) %>% unlist %>% matrix(, nrow = 2, ncol = length(theta1)) %>% t()

# Visual key for misfit
fit <- rep("<.95", times = length(theta1))
fit[logL < quant[,2]] <- ">.95"
fit <- ordered(fit, c(">.95", "<.95"))

# Set up and plot
gg <- data.frame(-2*sim$logL, -2*rep(logL, each = n_reps), rep(fit, each = n_reps), rep(-2*quant[,1], each = n_reps))

names(gg) <- c("l_dist", "l_obs", "fit", "median")
gg <- gg[order(gg$median), ]
gg$pair <- rep(1:length(theta1), each = n_reps)

ggplot(gg, aes(x = pair, y = l_dist, group = pair)) +
  geom_boxplot(outlier.size = 0, outlier.color = "grey90", aes(fill = fit)) +
  geom_point(aes(x = pair, y = l_obs, pch = fit, size = fit)) +
  scale_shape_manual(values = c(4, 20)) +
  scale_fill_grey(start = 0.1, end = 0.8) +
  xlab("Groups") +
  ylab("-2 * loglikelihood") +
  scale_size_manual(values = c(4, 1)) +
  theme(legend.title=element_blank())




# ------------------------------------------------------------
# Parameter recovery
# ------------------------------------------------------------
em_priors_a <-  rbind(em_a$prior, em_a$se) %>% round(3)
em_priors_b <-  rbind(em_b$prior, em_b$se) %>% round(3)
em_priors_c <-  rbind(em_c$prior, em_c$se) %>% round(3)

pv_priors_a <- pv_priors(em_pv_a) %>% round(3)
pv_priors_b <- pv_priors(em_pv_b) %>% round(3)
pv_priors_c <- pv_priors(em_pv_c) %>% round(3)

xtable::xtable(rbind(
  paste0(em_priors_a[1,], " (", em_priors_a[2,], ")"),
  paste0(pv_priors_a[1,], " (", pv_priors_a[4,], ")"),
  paste0(pv_priors_a[5,]),
  paste0(em_priors_b[1,], " (", em_priors_b[2,], ")"),
  paste0(pv_priors_b[1,], " (", pv_priors_b[4,], ")"),
  paste0(pv_priors_b[5,]),
  paste0(em_priors_c[1,], " (", em_priors_c[2,], ")"),
  paste0(pv_priors_c[1,], " (", pv_priors_c[4,], ")"),
  paste0(pv_priors_c[5,])))

# ------------------------------------------------------------
# Process loss with pv
# ------------------------------------------------------------
pv_a[models] <- pv_posteriors(em_pv_a, rep_order)
pv_b[models] <- pv_posteriors(em_pv_b, rep_order)
pv_c[models] <- pv_posteriors(em_pv_c, rep_order)

pv_a$e_mix <- e_mix(pv_a[models], parms, pv_a$theta1, pv_a$theta2, sorted = T)
pv_b$e_mix <- e_mix(pv_b[models], parms, pv_b$theta1, pv_b$theta2, sorted = T)
pv_c$e_mix <- e_mix(pv_c[models], parms, pv_c$theta1, pv_c$theta2, sorted = T)

pv <- rbind(pv_a, pv_b, pv_c)
pv$pl <-  1 - PL(pv$e_mix, parms, pv$theta1, pv$theta2, normed = T)
pv$pl2 <-  PL(pv$e_mix, parms, pv$theta1, pv$theta2, normed = F)
pv$Condition <- rep(Conditions, each = n_obs*n_reps)
pv$Model <- Models[pv$model]
pv$int <- interaction(pv$Condition, pv$Model)

ggplot(pv, aes(x = int, y = pl)) +
  geom_boxplot(outlier.shape = NA, outlier.color = "white", size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  ylab("1 - Process Loss") +
  xlab("Condition by Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(pv, aes(x = int, y = pl2)) +
  geom_boxplot(outlier.shape = NA, size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  ylab("Process Loss") +
  xlab("Condition by Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(c(0,120))

# ------------------------------------------------------------
# Process loss with original data
# ------------------------------------------------------------

A$e_mix <- e_mix(em_a$posterior, parms, theta[a1], theta[a2])
B$e_mix <- e_mix(em_b$posterior, parms, theta[b1], theta[b2])
C$e_mix <- e_mix(em_c$posterior, parms, theta[c1], theta[c2])

D <- rbind(A, B, C)
D$pl <-  1 - PL(D$e_mix, parms, D$theta1, D$theta2, normed = T)
D$pl2 <-  PL(D$e_mix, parms, D$theta1, D$theta2, normed = F)
D$Condition <- rep(Conditions, each = n_obs)
D$Model <- Models[D$model]
D$int <- interaction(D$Condition, D$Model)

ggplot(D, aes(x = int, y = pl)) +
  geom_boxplot(outlier.shape = NA, outlier.color = "white", size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  #  scale_fill_manual(values = c("grey90", "grey60", "grey30")) +
  ylab("1 - Process Loss") +
  xlab("Condition by Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(D, aes(x = int, y = pl2)) +
  geom_boxplot(outlier.shape = NA, size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  ylab("Process Loss") +
  xlab("Condition by Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(c(0,120))

# ------------------------------------------------------------
# Test redundancy
# Use pv or orginal data?
  # data <- D
    data <- pv
# ------------------------------------------------------------

deltas <- item_delta(parms, data$theta1, data$theta2)
delta_order <- apply(deltas, 1, order) %>% t
n_tests <- n_items - test_length

fun <- function(j, normed) {
  screen <- deltas*NA
  i <- j - test_length + 1
  ind <- cbind(rep(1:nrow(data), times = test_length), c(delta_order[,i:j]))
  screen[ind] <- 1
  redundancy(parms, data$theta1, data$theta2, normed = normed, NA_pattern = screen)
}

temp1 <- parallel::mclapply(test_length:(n_items-1), fun, normed = F)
temp2 <- parallel::mclapply(test_length:(n_items-1), fun, normed = T)
gg <- unlist(temp1) %>% data.frame
names(gg) <- "redun"
gg$pl <- 1- unlist(temp2)
gg$test <- as.ordered(rep(1:n_tests, each = nrow(data)))
gg$test2 <- as.ordered(rep(1:(n_tests/5)*5, each = nrow(data)*5))
gg$Condition <- rep(data$Condition, times = n_tests)
gg$Model <- data$Model

ggplot(gg, aes(x = test2, y = redun)) +
  geom_boxplot(outlier.shape = NA, size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(breaks=seq(0,75,5)) +
  ylab("Process Loss") +
  xlab("Test")

ggplot(gg, aes(x = test2, y = pl)) +
  geom_boxplot(outlier.shape = NA, size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(breaks=seq(0,75,5)) +
  ylab("1 - Process Loss") +
  xlab("Test")
