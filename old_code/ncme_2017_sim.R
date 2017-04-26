# ------------------------------------------------------------
# Data simulation for NCME 2017 presentation: IRT-based models for online tasks that involve student collaboration By Peter Halpin and Yoav Bergner.
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")

library("ggplot2")
library("dplyr")
source("~/github/cirt/R/cIRF_functions.R")
source("~/github/cirt/R/IRF_functions.R")
NYU <- rgb(87, 6, 140, maxColorValue = 255)
# ------------------------------------------------------------
# Data simulation for weighted addtive model
# ------------------------------------------------------------
set.seed(101) # Use to replicate NCME results

# Data generating parameters
n_obs <- 100
n_items <- 100
n_reps <- 250
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

# Simulate data
theta1 <- theta[odd]
theta2 <- theta[odd + 1]
theta1_se <- theta_se[odd]
theta2_se <- theta_se[odd+1]
data <- data_gen(1, w, parms, theta1, theta2)
resp <- data[item_names]

# Estimate weights for generated data
ml <- mle_WA(resp, parms, theta1, theta2, SE = "exp", starts = w)

# Generate plausible values from simuluated data
pv_data <- pv_gen(n_reps, resp, parms, theta1, theta2, theta1_se, theta2_se, weights = w)

# Esimate weights using replicated data and PV for data
ml_pv <- mle_WA(pv_data[item_names], parms, pv_data$theta1, pv_data$theta2, starts = pv_data$w, SE = "exp", parallel = T)

pv <- cbind(pv_data, ml_pv)


# ------------------------------------------------------------
# Parameter recovery
# ------------------------------------------------------------


# Means: Figure 1 ----------
pv_w <- tapply(pv$w, pv$pairs, mean)
gg_w <- data.frame(w, ml$w, pv_w)
names(gg_w) <- c("true_values", "ML", "PV")


w1 <- ggplot(gg_w, aes(x = true_values, y = ML)) + geom_point(pch = 1, size = 2) +
  ylim(c(0,1)) + xlim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_set(theme_grey(base_size = 15))

w2 <-ggplot(gg_w, aes(x = true_values, y = PV)) + geom_point(pch = 1, size = 2) +
  ylim(c(0,1)) + xlim(c(0,1)) +
    geom_abline(intercept = 0, slope = 1) +
    theme_set(theme_grey(base_size = 15))

w3 <-ggplot(gg_w, aes(x = ML, y = PV)) + geom_point(pch = 1, size = 2) +
  ylim(c(0,1)) + xlim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_set(theme_grey(base_size = 15))


t <- grid::textGrob("")
grid.arrange(w1, t, w2, w3, ncol=2)


# SEs: Figure 2 ---------------------
var_w <- tapply(pv$se^2, pv$pairs, mean)
var_b <- tapply(pv$w, pv$pairs, var)
pv_se <- sqrt(var_w + (1 + 1/n_reps) * var_b)
var_increase <- (1 + 1/n_reps) * var_b/var_w
var_prop <- (1 + 1/n_reps) * var_b/pv_se^2

gg_v <- data.frame(c(ml$se, pv_se), rep(pv_w, times = n_obs*2), rep(c("ML", "PV"), each = n_obs))

names(gg_v) <- c("SE", "PV_W", "Type")
head(gg_v)

ggplot(gg_v, aes(x = PV_W, y = SE, group = Type)) + geom_point(size = 3, aes(shape = Type, color = Type)) + ylim(c(0,.3)) +
    scale_color_manual(values = c("grey50", "black")) +
    scale_shape_discrete(solid=F) +
    theme_set(theme_grey(base_size = 15))

summary(var_prop)
summary(var_increase)
