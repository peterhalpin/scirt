# ------------------------------------------------------------
# Data simulation for NCME 2017 presentation: IRT-based models for online tasks that involve student collaboration By Peter Halpin and Yoav Bergner.

# Compare bias and MSE of w, theta, with MAP and ML, with small (20 items) and large (100 items) samples. n items on pretest is
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")

library("ggplot2")
library("dplyr")
source("~/github/cirt/R/cIRF_functions.R")
source("~/github/cirt/R/IRF_functions.R")
source("~/github/cirt/old_code/functions_WAI2.R")

NYU <- rgb(87, 6, 140, maxColorValue = 255)
# ------------------------------------------------------------
# Data simulation for weighted addtive model
# ------------------------------------------------------------
set.seed(101) # Use to replicate NCME results

# Data generating parameters
n_obs <- 500
n_items <- 100
alpha <- runif(n_items, .65, 2.5)
beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
parms <- data.frame(alpha, beta)
theta <- rnorm(n_obs*2)
a <- rnorm(n_obs)
w <- exp(a)/(1+exp(a))

# Other stuff
item_names <- paste0("item", 1:n_items)
row.names(parms) <- item_names
odd <- seq(1, n_obs*2, by = 2)

# Simulate data
theta1 <- theta[odd]
theta2 <- theta[odd + 1]

temp_col <- data_gen(1, w, parms, theta1, theta2)
names(temp_col)[-c(1:5)] <- paste0(names(temp_col)[-c(1:5)], "_COL")
temp_col<- temp_col[rep(1:nrow(temp_col), each = 2), ]

temp_ind <- data_gen(1, rep(.5, length(theta)), parms, theta, theta)
names(temp_ind) <- paste0(names(temp_ind), "_IND")


data <- cbind(temp_col, temp_ind[, -c(1:5)])
head(data)

# Information function for theta

# full data
col_items <- grep("COL", names(data))
ind_items <- grep("IND", names(data))
items1 <- c(col_items, ind_items)
ml1 <- est_WA(data[items], parms, method = "ml", SE = "exp")
map1 <- est_WA(data[items], parms, method = "map", SE = "exp")

plot(ml1$w, map1$w)
abline(a = 0, b = 1)

# short collab  data
items2 <- c(col_items[sample.int(100, 20)], ind_items)
temp2 <- data
temp2[-items2] <- NA
ml2 <- est_WA(temp2, parms, method = "ml", SE = "exp")
map2 <- est_WA(temp2, parms, method = "map", SE = "exp")
plot(data$w[odd], ml2$w)
mean(data$w[odd] - ml2$w)

plot(data$w[odd], map2$w)
mean(data$w[odd] - map2$w)

plot(ml2$w, map2$w)
abline(a = 0, b = 1)

# short ind  data
items3 <- c(col_items, ind_items[sample.int(100, 20)])
temp3 <- data
temp3[-items3] <- NA
ml3 <- est_WA(temp3, parms, method = "ml", SE = "exp")
map3 <- est_WA(temp3, parms, method = "map", SE = "exp")

plot(data$w[odd], ml3$w)
mean(data$w[odd] - ml3$w)

plot(data$w[odd], map3$w)
mean(data$w[odd] - map3$w)

plot(ml3$w, map3$w)
abline(a = 0, b = 1)

plot(ml3$w, map3$w)

gg <- c(map1$w)


# short both
items4 <- c(col_items[sample.int(100, 20)], ind_items[sample.int(100, 20)])
temp4 <- data
temp4[-items4] <- NA
ml4 <- est_WA(temp4, parms, method = "ml", SE = "exp")
map4 <- est_WA(temp4, parms, method = "map", SE = "exp")

plot(data$w[odd], ml4$w)
mean(data$w[odd] - ml4$w)

plot(data$w[odd], map4$w)
mean(data$w[odd] - map4$w)

plot(ml4$w, map4$w)
abline(a = 0, b = 1)


head(ml_obs)
plot(data$w[odd], ml_obs$w)
plot(data$w[odd], map_obs$w)
plot(ml_obs$w, map_obs$w)
cor(ml_obs$w, map_obs$w)

plot(map_obs$w, ml_obs$w_se)
points(map_obs$w, ml_exp$w_se, col = 5)
plot(ml_obs$w_se, ml_exp$w_se)
mean(ml_obs$w_se)
mean(ml_exp$w_se)

plot(map_obs$w, map_obs$w_se)
points(map_exp$w, map_exp$w_se, col = 5)
points(data$w[odd], map_exp$w_se, col = 5)

plot(ml_exp$w, ml_exp$w_se)
points(map_exp$w, map_exp$w_se, col = 5)
plot(ml_exp$w_se, map_exp$w_se)
abline(a = 0, b = 1)


# bias
mean(ml_obs$w - data$w)
mean(map_obs$w - data$w)

# variance
sqrt(mean(ml_obs$w_se^2))
sqrt(mean(ml_exp$w_se^2))
sqrt(mean(map_obs$w_se^2))
sqrt(mean(map_exp$w_se^2))





# Estimate weights for generated data
ml <- mle_WA(resp, parms, theta1, theta2, SE = "exp", starts = w)
map <- map_WA(resp, parms, theta1, theta2, SE = "exp", starts = w)
ml
plot(ml$w, w)
plot(map$w, w)
plot(ml$w, map$w, xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1)
cor(ml$w, map$w)
plot(ml$w, ml$se)
plot(map$w, map$psd)
plot(ml$se, map$psd)
 xlim = min(, .18), ylim = c(0,.2))



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
