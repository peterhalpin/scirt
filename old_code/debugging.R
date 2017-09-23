# ------------------------------------------------------------
# Debugging
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")
# devtools::use_data(sim_parms)
# source("~/github/sc_irt/R/cIRF_functions.R")
# source("~/github/sc_irt/R/IRF_functions.R")

# ------------------------------------------------------------
# Data simulation for weighted addtive model
# ------------------------------------------------------------

# Data generating parameters
n_obs <- 100 # n respondents
n_items <- 10 # n items
K <- n_obs/2 # n groups
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

# Generate group data
col_data <- data_gen(1, w, col_parms, theta1, theta2)
ind_data <- data_gen(1, rep(.5, n_obs), ind_parms, theta, theta)
data <- cbind(col_data[rep(1:K, each = 2),], ind_data[, -c(1:5)])
head(data)

# Only include item responses
ind_resp <- ind_data[,-(1:5)]
col_resp <- col_data[,-(1:5)]
resp <- data[,-(1:5)]

# ------------------------------------------------------------
# checking functions
# ------------------------------------------------------------

ind_resp <- ind_data[,-(1:5)]
m1 <- est_2PL(ind_resp, ind_parms, "MAP")
est_2PL(ind_resp, ind_parms, "MAP", parallel = F)

est_RSC2(col_resp, col_parms, m1$theta[odd], m1$theta[odd+1], method = "MAP", obs = F, epsilon = .05, parallel = T)

resp <- col_resp
parms <- col_parms
theta1 <- m1$theta[odd]
theta2 <-m1$theta[odd+1]

m2 <- est_RSC2(col_resp, col_parms, m1$theta[odd], m1$theta[odd+1], method = "ML", obs = T, epsilon = .05, parallel = T)

est_RSC(resp, parms, starts = NULL, method = "MAP", obs = F, epsilon = .05, parallel = F)
est_RSC(resp, parms, starts = NULL, method = "ML", obs = F, epsilon = .05, parallel = F)
est_RSC(resp, parms, starts = NULL, method = "ML", obs = F, epsilon = .05, parallel = T)
m3 <- est_RSC(resp, parms, starts = NULL, method = "ML", obs = T, epsilon = .05, parallel = F)

plot(m2$w, m3$w)

