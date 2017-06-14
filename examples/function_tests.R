# ------------------------------------------------------------
# Tests for functions
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")
rm(list = ls())
library("ggplot2")
library("dplyr")
source("~/github/cirt/R/IRF_functions.R")
source("~/github/cirt/R/cIRF_functions.R")

# ------------------------------------------------------------
# Parameter recovery with est_2PL
# ------------------------------------------------------------

# Data generating parameters
n_items <- 100
n_obs <- 1000
set.seed(101)

# Item parms
alpha <- runif(n_items, .65, 2.5)
beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
parms <- data.frame(alpha, beta)
names <- paste0("item", 1:n_items)
row.names(parms) <- names
names(parms) <- c("alpha", "beta")

# Respondents
theta <- rnorm(n_obs)

# Generate group data
resp <- sim_2PL(parms, theta)

# Estimation
ml <- est_2PL(resp, parms, method = "ML", parallel = F)
wml <- est_2PL(resp, parms, method = "WML", parallel = F)
map <- est_2PL(resp, parms, method = "MAP", parallel = F)

# Recovery is OK
cor(cbind(theta, ml$theta, wml$theta, map$theta))
plot(wml$theta, map$theta)

# Parallel doesnt make much of a difference ...
system.time(est_2PL(resp, parms, method = "MAP", parallel = F))
system.time(est_2PL(resp, parms, method = "MAP", parallel = T))

# ------------------------------------------------------------
# Parameter recovery with est_RSC
# Note that use of I, J
# ------------------------------------------------------------

# Data generating parameters
n_obs <- 100 # n respondents
K <- n_obs/2 # n groups
n_items <- 20 # n items on each assessment
e <- .05 # beta prior parm
set.seed(101)

# Individual test
ind_alpha <- runif(n_items, .65, 2.5)
ind_beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
ind_parms <- data.frame(ind_alpha, ind_beta)
ind_names <- paste0("item", 1:n_items, "_IND")
row.names(ind_parms) <- ind_names
names(ind_parms) <- c("alpha", "beta")

# Group test
col_alpha <- runif(n_items, .65, 2.5)
col_beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
col_parms <- data.frame(col_alpha, col_beta)
col_names <- paste0("item", 1:n_items, "_COL")
row.names(col_parms) <- col_names
names(col_parms) <- c("alpha", "beta")
parms <- rbind(col_parms, ind_parms)

# Respondents
odd <- seq(1, n_obs, by = 2)
theta <- rnorm(n_obs)
theta1 <- theta[odd] # odds
theta2 <- theta[odd + 1] # evens

# RSC parameter
w <- rbeta(K, 1 + e, 1 + e)

# Generate data
col_data <- data_gen(1, w, col_parms, theta1, theta2)
ind_data <- data_gen(1, rep(.5, n_obs), ind_parms, theta, theta)
data <- cbind(col_data[rep(1:K, each = 2),], ind_data[, -c(1:5)])
head(data)

# Estimate
ml <- est_RSC(data, parms, method = "ML", parallel = F)
map <- est_RSC(data, parms, method = "MAP", parallel = F)
plot(ml$w_se, map$w_se)
abline(a = 0, b = 1)
# ------------------------------------------------------------
# Parameter recovery with est_RSC2
# ------------------------------------------------------------

# Subset data
data2 <- data[grep("COL", names(data))][odd, ]
# Estimate
ml2 <- est_RSC2(data2, col_parms, theta1, theta2, method = "ML", parallel = F)
map2 <- est_RSC2(data2, col_parms, theta1, theta2, method = "MAP", parallel = F)
