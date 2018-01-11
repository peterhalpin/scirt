

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
sigma <- 2 # logit prior parm
set.seed(111)

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
hist(dgp_w, breaks = 50)

# Generate selected group data
col_data <- data_gen(1, u, col_parms, theta1, theta2)
ind_data <- data_gen(1, rep(0, n_obs), ind_parms, theta, theta)

col_test <- col_data[rep(1:K, each = 2), -c(1:5)]
ind_test <- ind_data[, -c(1:5)]

parms <- rbind(col_parms, ind_parms)
data <- cbind(col_test, ind_test)


z <- function(u, theta1, theta2, alpha, stat = "synergy") {
  Theta <- cbind(theta1, theta2)
  theta_min <- apply(Theta, 1, min)
  theta_max <- apply(Theta, 1, max)
  u - mean(alpha) * (theta_max - theta_min)
}

sigma = 3
plot(u, z_hat)
hist(z_hat)
resp <- data
diag(var)
var <- var_RSC(resp, u, parms, theta1, theta2, method, obs, sigma)

z1 <- z_test(data, map_ll$u, parms, map_ll$theta1, map_ll$theta2, method = "MAP", obs = T, sigma = 3)
z1[z1$z_hat > 0,]
head(map_ll)
z1[z_hat > 0, ]

plot(map_ll$u_se, z1$z_se)
abline(a= 0, b = 1)

z_test <- function(resp, u, parms, theta1, theta2, method = "ML", obs = T, sigma = 3, stat = "synergy") {

   alpha_bar <- mean(parms[grep("COL", row.names(parms)), "alpha"])
   Theta <- cbind(theta1, theta2)
   theta_min <- apply(Theta, 1, which.min, arr.ind = T)
   theta_max <- apply(Theta, 1, which.max)
   z_hat <-  u - alpha_bar * (Theta[,which.max] - theta_min)
   z_se <- z_hat*0
   var <- var_RSC(resp, u, parms, theta_max, theta_min, method, obs, sigma)

   if (stat == "antagonism") {
     r <- c(1, -alpha_bar, alpha_bar)
   } else {
     r <- c(1, alpha_bar, -alpha_bar)
   }
   R <- matrix(r, nrow = 1, ncol = 3)
   ind <- seq(1, nrow(var), by = 3)

   for (i in 1:length(ind)) {
     j <- ind[i]
     k <- j + 2
     temp <- var[j:k, j:k]
     z_se[i] <-  R %*% temp %*% t(R) %>% sqrt %>% as.numeric
  }
  data.frame(z_hat, z_se)
}

warnings()
A <- matrix(1:64, nrow = 8, ncol = 8)
d <-
abs(row(A) - col(A))

 > 2
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
p_factor <- 5
qnorm(.005, 0, p_factor) %>% inv_logit
qnorm(.995, 0, p_factor) %>% inv_logit


ml_ll <-  est_RSC(data[ll_items], parms[ll_items, ], method = "ML", obs = T)
map_ll <-  est_RSC(data[ll_items], parms[ll_items, ], method = "MAP", obs = T, sigma = sigma)

# map_ll2 <-  est_RSC(data[ll_items], parms[ll_items, ], method = "MAP", obs = T, sigma = sigma)
#
#
# plot(theta1, map_ll$theta1)
# abline(a = 0, b = 1)
# lines(lowess(theta1, map_ll$theta1), col = 3)
# plot(dgp_u, map_ll$u)
# abline(a = 0, b = 1)
# lines(lowess(dgp_u, map_ll$u), col = 3)

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

# ------------------------------------------------------------
# Plots



#ml_ss <- ml_sl <- ml_ls <-ml_ll
#map_ss <- map_sl <- map_ls <-map_ll
u_hat <- c(ml_ll$u, map_ll$u , ml_sl$u, map_sl$u , ml_ls$u, map_ls$u , ml_ss$u, map_ss$u)
u_se <- c(ml_ll$u_se, map_ll$u_se , ml_sl$u_se, map_sl$u_se , ml_ls$u_se, map_ls$u_se , ml_ss$u_se, map_ss$u_se)


ind_form <- rep(c("Individual long", "Indvidual short"), each = K*4)
col_form <- rep(c("Group long", "Group short"), each = K*2) %>% rep(times = 2)
Method <- rep(c("ML", "MAP"), each = K) %>% rep(times = 4)
dgp_u <- rep(u, times = 8)

gg <- data.frame(u_hat, u_se, ind_form, col_form, Method, dgp_u)

#gg$sample[sample(nrow(gg), nrow(gg)/10)] <- 1

gg$sample <- 0
breaks <- seq(from = min(dgp_u), to = max(dgp_u), length.out = 10)
for (i in 1:length(breaks)) {
  ind <- which(gg$dgp_u > breaks[i] & gg$dgp_u < breaks[i+1])
  gg$sample[sample(ind, 30, replace = T)] <- 1
}

#p1 <- ggplot(gg, aes(x = inv_logit(dgp_u), y = inv_logit(u_hat), group = Method)) +
p1 <- ggplot(gg, aes(x = dgp_u, y = u_hat, group = Method)) +
    geom_smooth(se = F, size = .8, method = "loess", span = 1, aes(linetype = Method, color = Method)) +
    geom_point(data = gg[gg$sample == 1,], size = 3, aes(shape = Method, color = Method)) +
    scale_color_manual(values = c("grey10", "grey10")) +
    scale_shape_discrete(solid=F) +
    xlab("Data generating values") +
    ylab("Estimate")  +
    geom_abline(slope = 1, intercept = 0, col = "grey50", size = 1.2, alpha = .5) +
    theme_bw(base_size = 15)

p1 + facet_grid(ind_form ~ col_form)

p2 <- ggplot(gg, aes(x = dgp_u, y = u_se, group = Method)) +
    geom_point(data = gg[gg$sample == 1,], size = 3, aes(shape = Method, color = Method)) +
    geom_smooth(se = F, size = .8, aes(linetype = Method, color = Method)) +
    scale_color_manual(values = c("grey10", "grey10")) +
    scale_shape_discrete(solid=F) +
    xlab("Estimate") +
    ylab("Standard error")  +
    theme_bw(base_size = 15) +
    coord_cartesian(ylim=c(0,5), xlim=c(-5,5))

p2 + facet_grid(ind_form ~ col_form)
