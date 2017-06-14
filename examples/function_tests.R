# ------------------------------------------------------------
# Tests for functions
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")

library("ggplot2")
library("dplyr")
source("~/github/cirt/R/IRF_functions.R")

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
