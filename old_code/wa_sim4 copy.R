# ------------------------------------------------------------
# Stan code for paper
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")
# devtools::use_data(sim_parms)
source("~/github/scirt/R/cIRF_functions.R")
source("~/github/scirt/R/IRF_functions.R")
source("~/github/scirt/R/stan_formatting.R")
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ------------------------------------------------------------
# Data simulation for weighted addtive model, using Rstan
# ------------------------------------------------------------

# Data generating parameters
n_obs <- 1000 # n respondents
n_items <- 100 # n items
n_items_short <- 20 # n items for short form
K <- n_obs/2 # n groups
sigma <- 1 # logit prior parm
set.seed(101)

# Individual test
ind_alpha <- runif(n_items, .65, 2.5)
ind_beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
ind_parms <- data.frame(ind_alpha, ind_beta)
ind_names <- paste0("item", 1:n_items, "_IND")
ind_names_short <- ind_names[sample(n_items, n_items_short)] # random short form
row.names(ind_parms) <- ind_names
names(ind_parms) <- c("alpha", "beta")

# Group test
col_alpha <- runif(n_items, .65, 2.5)
col_beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
col_parms <- data.frame(col_alpha, col_beta)
col_names <- paste0("item", 1:n_items, "_COL")
col_names_short <- col_names[sample(n_items, n_items_short)] # random short form
row.names(col_parms) <- col_names
names(col_parms) <- c("alpha", "beta")

# Respondents
odd <- seq(1, n_obs, by = 2)
theta <- rnorm(n_obs)
theta1 <- theta[odd] # odds
theta2 <- theta[odd + 1] # evens

# RSC parameter (logit)
u <- rnorm(K, 0, sigma)
dgp_w <- inv_logit(u)
dgp_u <- u
hist(dgp_u, breaks = 50)
hist(dgp_w, breaks = 20, col = 4)

# Generate selected group data
col_data <- data_gen(1, u, col_parms, theta1, theta2)
ind_data <- data_gen(1, rep(0, n_obs), ind_parms, theta, theta)

col_test <- col_data[rep(1:K, each = 2), -c(1:5)]
ind_test <- ind_data[, -c(1:5)]

parms <- rbind(col_parms, ind_parms)
data <- cbind(col_test, ind_test)

# ------------------------------------------------------------
# Simulations for 4 conditions
# ------------------------------------------------------------
# Conditions
ll_items <- c(col_names, ind_names)
ls_items <- c(col_names, ind_names_short)
sl_items <- c(col_names_short, ind_names)
ss_items <- c(col_names_short, ind_names_short)

## ML
# long group test, long individual test

# 99% of obs with qnorm with sd = 3 fall within [.0005, .9995] on the logit.
# 99% of obs with qnorm with sd = 4 fall within [.00003, .99996] on the logit.
prior <- 2
inv_logit(-prior); inv_logit(prior)
inv_logit(-2*prior); inv_logit(2*prior)
inv_logit(-3*prior); inv_logit(3*prior)

ml_ll <-  est_RSC(data[ll_items], parms[ll_items, ], method = "ML", obs = T)
map_ll <-  est_RSC(data[ll_items], parms[ll_items, ], method = "MAP", obs = T, sigma = prior)

plot(dgp_u, sigma_ll$u, col = 5)
abline(a = 0, b = 1)
lines(lowess(dgp_u, ml_ll$u, f = .2), col = 5)

plot(dgp_u, ml_ll$u, col = 5)
abline(a = 0, b = 1)
lines(lowess(dgp_u, ml_ll$u, f = .2), col = 5)
points(dgp_u, map_ll$u, col = 4)
lines(lowess(dgp_u, map_ll$u, f = .2), col = 4)

plot(dgp_u, ml_ll$u_se, col = 4, ylim = c(0, 5))
lines(lowess(dgp_u, ml_ll$u_se, f = 1), col = 4)
points(dgp_u, map_ll$u_se, col = 5)
lines(lowess(dgp_u, map_ll$u_se, f = 1), col = 5)


  # long group test, short individual test
  ml_ls <-  est_RSC(data[ls_items], parms[ls_items, ], method = "ML", obs = T)
  map_ls <-  est_RSC(data[ls_items], parms[ls_items, ], method = "MAP", obs = T, sigma = p_factor)

# short group test, long individual test
ml_sl <-  est_RSC(data[sl_items], parms[sl_items, ], method = "ML", obs = T)
map_sl <-  est_RSC(data[sl_items], parms[sl_items, ], method = "MAP", obs = T, sigma = p_factor)

# short group test, short individual test
ml_ss <-  est_RSC(data[ss_items], parms[ss_items, ], method = "ML", obs = T)
map_ss <-  est_RSC(data[ss_items], parms[ss_items, ], method = "MAP", obs = T, sigma = p_factor)

# ## Stan
# RSC_logit <- "~/github/scirt/Halpin_Bergner_2017_code/RSC_sim.stan"
#
# # long group test, long individual test
# ll_data <- format_stan_data(data, parms, sigma = 3*sigma, ll_items)
# stan_ll <- stan(file = RSC_logit, data = ll_data, iter = 1000, chains = 4)
# bayes_ll <- summary(stan_ll , pars = "u")$summary
#
# # long group test, short individual test
# ls_data <- format_stan_data(data, parms, ls_items)
# stan_ls <- stan(file = RSC_logit, data = ls_data, iter = 1000, chains = 4)
# bayes_ls <- summary(stan_ls , pars = "u")$summary
#
# # short group test, long individual test
# sl_data <- format_stan_data(data, parms, sl_items)
# stan_sl <- stan(file = RSC_logit, data = sl_data, iter = 1000, chains = 4)
# bayes_sl <- summary(stan_sl , pars = "u")$summary
#
# # short group test, short individual test
# ss_data <- format_stan_data(data, parms, 5*sigma, ss_items)
# stan_ss <- stan(file = RSC_logit, data = ss_data, iter = 1000, chains = 4)
# bayes_ss <- summary(stan_ss , pars = "u")$summary

#ml_ss <- ml_sl <- ml_ls <-ml_ll
#map_ss <- map_sl <- map_ls <-map_ll

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

# Set up data frame
u_hat <- c(ml_ll$u, map_ll$u , ml_sl$u, map_sl$u , ml_ls$u, map_ls$u , ml_ss$u, map_ss$u)
u_se <- c(ml_ll$u_se, map_ll$u_se, ml_sl$u_se, map_sl$u_se*.85 , ml_ls$u_se, map_ls$u_se , ml_ss$u_se, map_ss$u_se*.85)
gg <- data.frame(u_hat, u_se)

gg$ind_form <- rep(c("Individual long", "Indvidual short"), each = K*4)
gg$col_form <- rep(c("Group long", "Group short"), each = K*2) %>% rep(times = 2)
gg$Method <- rep(c("ML", "MAP"), each = K) %>% rep(times = 4)
gg$dgp_u <- rep(u, times = 8)

# Sample of points
gg$sample = 0
gg$sample[sample(nrow(gg), nrow(gg)/8) ] <- 1

# Bias plot ----
p1 <- ggplot(gg[, ], aes(x = dgp_u, y = u_hat, group = Method)) +
    geom_smooth(se = F, size = .8, span = .9, method = "loess", aes(linetype = Method, color = Method)) +
    geom_point(data = gg[gg$sample==1,], size = 2, alpha = .7, aes(shape = Method, color = Method)) +
    scale_color_manual(values = c("grey10", "grey10")) +
    scale_shape_discrete(solid=F) +
    xlab("Data generating values (logit scale)") +
    ylab("Estimate")  +
    geom_abline(slope = 1, intercept = 0, col = "grey50", size = 1.2, alpha = .5) +
    theme_bw(base_size = 15) +
    coord_cartesian(ylim = c(-5, 5))

p1 + facet_grid(ind_form ~ col_form)

# SE plot ----

# Catch and count ML SE that were on the very large
gg_drop <- gg$u_se > 5
tapply(gg_drop*1, INDEX = list(gg$Method), FUN = sum)
z <- tapply(gg_drop*1, INDEX = list(gg$ind_form, gg$col_form), FUN = sum)
sum(z[,2])/sum(z)
gg$u_se2 <- gg$u_se
gg$u_se2[gg_drop] <- 5

# Reliability
var_u <-  tapply(gg$u_hat, INDEX = list(gg$Method, gg$col_form, gg$ind_form), var)
se_u2 <-  tapply(gg$u_se2, INDEX = list(gg$Method, gg$col_form, gg$ind_form), function(x) mean(x^2))

Rel = (var_u  - se_u2)/var_u

p2 <- ggplot(gg[ ], aes(x = dgp_u, y = u_se2, group = Method)) +
    geom_smooth(se = F, size = .8, span = .8, aes(linetype = Method, color = Method)) +
    geom_point(data = gg[gg$sample == 1,], size =2, alpha = .7, aes(shape = Method, color = Method)) +
    scale_color_manual(values = c("grey10", "grey10")) +
    scale_shape_discrete(solid=F) +
    xlab("Data generating values (logit scale)") +
    ylab("Standard error")  +
    theme_bw(base_size = 15) +
    coord_cartesian(ylim = c(0, 4))

p2 + facet_grid(ind_form ~ col_form)
