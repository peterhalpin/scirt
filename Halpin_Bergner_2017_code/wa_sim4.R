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
dgp_w <- p(u)
dgp_u <- u

# Generate selected group data
col_data <- data_gen(1, u, col_parms, theta1, theta2)
ind_data <- data_gen(1, rep(.5, n_obs), ind_parms, theta, theta)

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
ml_ll <-  est_RSC(data[ll_items], parms[ll_items, ], method = "ML", obs = T)

# long group test, short individual test
ml_ls <-  est_RSC(data[ls_items], parms[ls_items, ], method = "ML", obs = T)

# short group test, long individual test
ml_sl <-  est_RSC(data[sl_items], parms[sl_items, ], method = "ML", obs = T)

# short group test, short individual test
ml_ss <-  est_RSC(data[ss_items], parms[ss_items, ], method = "ML", obs = T)

## Stan
RSC_logit <- "~/github/scirt/Halpin_Bergner_2017_code/RSC_logit.stan"

# long group test, long individual test
ll_data <- format_stan_data(data, parms, ll_items)
stan_ll <- stan(file = RSC_logit, data = ll_data, iter = 1000, chains = 4)
bayes_ll <- summary(stan_ll , pars = "u")$summary

# long group test, short individual test
ls_data <- format_stan_data(data, parms, ls_items)
stan_ls <- stan(file = RSC_logit, data = ls_data, iter = 1000, chains = 4)
bayes_ls <- summary(stan_ls , pars = "u")$summary

# short group test, long individual test
sl_data <- format_stan_data(data, parms, sl_items)
stan_sl <- stan(file = RSC_logit, data = sl_data, iter = 1000, chains = 4)
bayes_sl <- summary(stan_sl , pars = "u")$summary

# short group test, short individual test
ss_data <- format_stan_data(data, parms, ss_items)
stan_ss <- stan(file = RSC_logit, data = ss_data, iter = 1000, chains = 4)
bayes_ss <- summary(stan_ss , pars = "u")$summary

save(ml_ll, stan_ll, ml_ls, stan_ls, ml_sl, stan_sl, ml_ss, stan_ss, file = "sim_results")

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------
u_hat <- c(ml_ll$w, bayes_ll[,"mean"], ml_sl$w, bayes_sl[,"mean"], ml_ls$w, bayes_ls[,"mean"], ml_ss$w, bayes_ss[,"mean"])

u_se <- c(ml_ll$w_se, bayes_ll[,"sd"], ml_sl$w_se, bayes_sl[,"sd"], ml_ls$w_se, bayes_ls[,"sd"], ml_ss$w_se, bayes_ss[,"sd"])

ind_form <- rep(c("Individual long", "Indvidual short"), each = K*4)
col_form <- rep(c("Group long", "Group short"), each = K*2) %>% rep(times = 2)
Method <- rep(c("ML", "Bayes"), each = K) %>% rep(times = 4)
dgp_u <- rep(dgp_u, times = 8)

gg <- data.frame(u_hat, u_se, ind_form, col_form, Method, dgp_u)

gg$sample <- 0
gg$sample[sample(nrow(gg), nrow(gg)/10)] <- 1
head(gg)

p1 <- ggplot(gg, aes(x = dgp_u, y = u_hat, group = Method)) +
    geom_point(data = gg[gg$sample == 1,], size = 3, aes(shape = Method, color = Method)) +
    geom_smooth(se = F, size = .8, aes(linetype = Method, color = Method)) +
    scale_color_manual(values = c("grey10", "grey10")) +
    scale_shape_discrete(solid=F) +
    xlab("Data generating values") +
    ylab("Estimate")  +
    geom_abline(slope = 1, intercept = 0, col = "grey50", size = 1.2) +
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
    ylim(0, 5)

p2 + facet_grid(ind_form ~ col_form)
