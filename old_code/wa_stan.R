# ------------------------------------------------------------
# Stan code for paper
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")
# devtools::use_data(sim_parms)
# source("~/github/scirt/R/cIRF_functions.R")
# source("~/github/scirt/R/IRF_functions.R")

# ------------------------------------------------------------
# Data simulation for weighted addtive model, using Rstan
# ------------------------------------------------------------

# Data generating parameters
n_obs <- 200 # n respondents
n_items <- 20 # n items
i <- 20 # n items for short form
K <- n_obs/2 # n groups
sigma <- 2 # beta prior parm
set.seed(101)

# Individual test
ind_alpha <- runif(n_items, .65, 2.5)
ind_beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
ind_parms <- data.frame(ind_alpha, ind_beta)
ind_names <- paste0("item", 1:n_items, "_IND")
ind_names_short <- ind_names[sample(n_items, i)] # random short form
row.names(ind_parms) <- ind_names
names(ind_parms) <- c("alpha", "beta")

# Group test
col_alpha <- runif(n_items, .65, 2.5)
col_beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
col_parms <- data.frame(col_alpha, col_beta)
col_names <- paste0("item", 1:n_items, "_COL")
col_names_short <- col_names[sample(n_items, i)] # random short form
row.names(col_parms) <- col_names
names(col_parms) <- c("alpha", "beta")

# Respondents
odd <- seq(1, n_obs, by = 2)
theta <- rnorm(n_obs)
theta1 <- theta[odd] # odds
theta2 <- theta[odd + 1] # evens

# RSC parameter (logit)
u <- rnorm(K, 0, sigma)
dgw <- p(u)
dgu <- u

# Generate selected group data
col_data <- data_gen(1, u, col_parms, theta1, theta2)
ind_data <- data_gen(1, rep(.5, n_obs), ind_parms, theta, theta)

col_test <- col_data[rep(1:K, each = 2), -c(1:5)]
ind_test <- ind_data[, -c(1:5)]

# ------------------------------------------------------------
# Data formatting for stan input
# ------------------------------------------------------------
library("rstan")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

K <- n_obs/2 # n groups
J_ind <- n_items # n items on ind form
J_col <- n_items # n items on col form

ind_test1 <- ind_test[odd, ]
temp <- which(!is.na(ind_test1), arr.ind = T)
y_ind1 <- ind_test1[temp]
kk_ind1 <- temp[,1]
jj_ind1 <- temp[,2]
N_ind1 <- length(y_ind1)

ind_test2 <- ind_test[odd + 1, ]
temp <- which(!is.na(ind_test2), arr.ind = T)
y_ind2 <- ind_test2[temp]
kk_ind2 <- temp[,1]
jj_ind2 <- temp[,2]
N_ind2 <- length(y_ind2)

col_test1 <- col_test[odd, ]
temp <- which(!is.na(col_test1), arr.ind = T)
y_col <- col_test1[temp]
kk_col <- temp[,1]
jj_col <- temp[,2]
N_col <- length(y_col)

alpha_ind <- ind_parms$alpha
beta_ind <- ind_parms$beta
alpha_col <- col_parms$alpha
beta_col <- col_parms$beta

stan_data <- list(K = K,
                  J_ind = J_ind,
                  J_col = J_col,
                  N_ind1 = N_ind1,
                  jj_ind1 = jj_ind1,
                  kk_ind1 = kk_ind1,
                  y_ind1 = y_ind1,
                  N_ind2 = N_ind2,
                  jj_ind2 = jj_ind2,
                  kk_ind2 = kk_ind2,
                  y_ind2 = y_ind2,
                  N_col = N_col,
                  jj_col = jj_col,
                  kk_col = kk_col,
                  y_col = y_col,
                  alpha_ind = alpha_ind,
                  beta_ind = beta_ind,
                  alpha_col = alpha_col,
                  beta_col = beta_col)

stan_file <- "~/github/scirt/Halpin_Bergner_2017_code/RSC.stan"
fit1 <- stan(file = stan_file, data = stan_data, iter = 2000, chains = 2)





# Bias
theta1_hat <- summary(fit1, pars = "theta1")$summary
plot(theta1_hat[,"mean"], theta1)
abline(a = 0, b = 1)

theta2_hat <- summary(fit1, pars = "theta2")$summary
plot(theta2_hat[,"mean"], theta2)
abline(a = 0, b = 1)

w_hat <- summary(fit1, pars = "w")$summary
plot(w_hat[,"mean"], w)
abline(a = 0, b = 1)
(var(w_hat[,"mean"]) - mean(w_hat[,"sd"]^2)) / var(w_hat[,"mean"])

u_hat <- summary(fit1, pars = "u")$summary
plot(u_hat[,"mean"], u)
abline(a = 0, b = 1)
(var(u_hat[,"mean"]) - mean(u_hat[,"sd"]^2)) / var(u_hat[,"mean"])
hist(u_hat[,"mean"])

plot(u_hat[, "mean"], map_ll[,"w"])
plot(u_hat[, "50%"], map_ll[,"w"])
abline(a = 0, b = 1)
plot(u_hat[, "sd"], map_ll[,"w_se"])
abline(a = 0, b = 1)


# SE for u and W
u_index <- order(u_hat[,"mean"])
plot(u_hat[u_index,"mean"], ylim = c(-5, 5))
points(u_hat[u_index,"2.5%"], col = 3)
points(u_hat[u_index,"97.5%"], col = 2)
abline(a = 0, b = 0)

w_index <- order(w_hat[,"mean"])
plot(w_hat[w_index,"mean"], ylim = c(-0, 1))
points(w_hat[w_index,"2.5%"], col = 3)
points(w_hat[w_index,"97.5%"], col = 2)
abline(a = .5, b = 0)

# what is this n_eff?

plot(u_hat[,"mean"], u_hat[,"n_eff"])
plot(w_hat[,"mean"], w_hat[,"n_eff"])

plot(w_hat[,"n_eff"], u_hat[,"n_eff"])
plot(abs(theta1_hat[,"mean"] - theta2_hat[,"mean"]), w_hat[,"n_eff"])


# ------------------------------------------------------------
# AMT example data
# ------------------------------------------------------------

setwd("~/Dropbox/Academic/Projects/CA/Data/response_matrices")
temp_parms <- read.csv("calibration_parms.csv", row.names = 1)

# Drop DIF items
#dif_items <- "045|065"
#temp_parms <- temp_parms[!grepl(dif_items, row.names(temp_parms)),]

# Item names without versions suffix for easy extraction
items <- paste0(row.names(temp_parms), collapse = "|")

# Individual and collaborative versions
ind_parms <- col_parms <- temp_parms
row.names(ind_parms) <- paste0(row.names(ind_parms), "_IND")
row.names(col_parms) <- paste0(row.names(col_parms), "_COL")

# Final parameter set
parms <- rbind(ind_parms, col_parms)


collab <- read.csv("collaboration_2016_0.csv", check.names = F)
col_form <- format_resp(collab, row.names(col_parms), "COL")
ind_form <- format_resp(collab, row.names(ind_parms), "IND")

# Apply conjunctive scoring rule to collaborative form
odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]

# Drop 13 unsuable response patterns (all 1 or all 0)
# (But are they unusable with MAP??)
drop_groups <- c(
  collab$group_id[apply(col_form, 1, mean, na.rm = T) %in% c(1,0)],
  collab$group_id[apply(ind_form, 1, mean, na.rm = T) %in% c(1,0)])

col_form <-col_form[collab$group_id%in%drop_groups,]
ind_form <-ind_form[!collab$group_id%in%drop_groups,]

Reset odd for dropped items
odd <- seq(1, nrow(col_form), by = 2)

drop_groups

# ------------------------------------------------------------
# Stan step up
# ------------------------------------------------------------

library("rstan")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_obs <- dim(ind_form)[1]
K <- n_obs/2 # n groups
J_ind <- dim(ind_parms)[1] # n items on ind form
J_col <- dim(col_parms)[1] # n items on col form

ind_test1 <- ind_form[odd, ]
temp <- which(!is.na(ind_test1), arr.ind = T)
y_ind1 <- ind_test1[temp]
kk_ind1 <- temp[,1]
jj_ind1 <- temp[,2]
N_ind1 <- length(y_ind1)

ind_test2 <- ind_form[odd + 1, ]
temp <- which(!is.na(ind_test2), arr.ind = T)
y_ind2 <- ind_test2[temp]
kk_ind2 <- temp[,1]
jj_ind2 <- temp[,2]
N_ind2 <- length(y_ind2)

col_test1 <- col_form[odd, ]
temp <- which(!is.na(col_test1), arr.ind = T)
y_col <- col_test1[temp]
kk_col <- temp[,1]
jj_col <- temp[,2]
N_col <- length(y_col)

alpha_ind <- ind_parms$alpha
beta_ind <- ind_parms$beta
alpha_col <- col_parms$alpha
beta_col <- col_parms$beta

amt_data <- list(K = K,
                  J_ind = J_ind,
                  J_col = J_col,
                  jj_ind1 = jj_ind1,
                  kk_ind1 = kk_ind1,
                  y_ind1 = y_ind1,
                  N_ind1 = N_ind1,
                  jj_ind2 = jj_ind2,
                  kk_ind2 = kk_ind2,
                  y_ind2 = y_ind2,
                  N_ind2 = N_ind2,
                  jj_col = jj_col,
                  kk_col = kk_col,
                  y_col = y_col,
                  N_col = N_col,
                  alpha_ind = alpha_ind,
                  beta_ind = beta_ind,
                  alpha_col = alpha_col,
                  beta_col = beta_col)

RSC_logit <- "~/github/scirt/Halpin_Bergner_2017_code/RSC_logit.stan"
RSC_beta <- "~/github/scirt/Halpin_Bergner_2017_code/RSC_beta.stan"

fit_logit <- stan(file = RSC_logit, data = amt_data, iter = 2000, chains = 4)
#fit_beta <- stan(file = RSC_beta, data = amt_data, iter = 2000, chains = 4)

w_pars <- paste0("w[", 1:25, "]")
u_pars <- paste0("u[", 1:25, "]")

# Bias
fit_beta <- amt_fit_beta
fit_logit <- amt_fit_logit
probs <- c(.005 .025, .05, .1, .9, .95, .975, .995)

theta1_beta <- summary(fit_beta, pars = "theta1")$summary
theta1_logit <- summary(fit_logit, pars = "theta1")$summary
plot(theta1_beta[,"mean"], theta1_logit[,"mean"])
abline(a = 0, b = 1)

theta2_beta <- summary(fit_beta, pars = "theta2")$summary
theta2_logit <- summary(fit_logit, pars = "theta2")$summary
plot(theta2_beta[,"mean"], theta2_logit[,"mean"])
abline(a = 0, b = 1)

w_beta <- summary(fit_beta, pars = "w", probs = probs)$summary
w_logit <- summary(fit_logit, pars = "w", probs = probs)$summary
plot(w_beta[,"mean"], w_logit[,"mean"])
abline(a = 0, b = 1)

(var(w_beta[,"mean"]) - mean(w_beta[,"sd"]^2)) / var(w_beta[,"mean"])
(var(w_logit[,"mean"]) - mean(w_logit[,"sd"]^2)) / var(w_logit[,"mean"])
hist(w_beta[,"mean"])
hist(w_logit[,"mean"])

u_logit <- summary(fit_logit, pars = "u", probs = probs)$summary
plot(u_logit[,"mean"], w_logit[,"mean"])
abline(a = 0, b = 1)
(var(u_logit[,"mean"]) - mean(u_logit[,"sd"]^2)) / var(u_logit[,"mean"])
hist(u_logit[,"mean"])

# SE for u and W
u_index <- order(u_logit[,"mean"])
plot(u_logit[u_index,"mean"], ylim = c(-7, 7))
points(u_logit[u_index,"5%"], col = 3)
points(u_logit[u_index,"95%"], col = 2)
abline(a = 0, b = 0)

level <- .95
const <- qnorm(level)
u_index <- order(u_logit[,"mean"])
plot(u_logit[u_index,"mean"], ylim = c(-7, 7))
points(u_logit[u_index,"mean"] - const*u_logit[u_index,"sd"], col = 3)
points(u_logit[u_index,"mean"] + const*u_logit[u_index,"sd"], col = 2)
abline(a = 0, b = 0)

w_index <- order(w_logit[,"mean"])
plot(w_logit[w_index,"mean"], ylim = c(-0, 1))
points(w_logit[w_index,"10%"], col = 3)
points(w_logit[w_index,"90%"], col = 2)
abline(a = .5, b = 0)

w_index2 <- order(w_beta[,"mean"])
plot(w_beta[w_index,"mean"], ylim = c(-0, 1))
points(w_beta[w_index,"10%"], col = 3)
points(w_beta[w_index,"90%"], col = 2)
abline(a = .5, b = 0)

# ------------------------------------------------------------
# GOF  Plot
# ------------------------------------------------------------

# Quantiles to extract from stan fit
probs <- c(.005, .025, .05, .1, .9, .95, .975, .995)
int_names <- paste0(rep(c("lower", "upper"), each = 4), c("_99", "_95", "_90", "_80",  "_80",  "_90",  "_95", "_99"))

# Set up df
log_lik <- summary(fit_logit, pars = "log_lik", probs = probs)$summary
l_obs <- -2*l_RSC(col_form[odd,], u_logit[,"mean"], col_parms, theta1_logit[,"mean"], theta1_logit[,"mean"])

gg <- data.frame(log_lik)
gg$obs <- l_obs

# Rename intervals
names(gg)[grep("X", names(gg))] <- int_names
head(gg)

# Visual indicator for fit
gg$fit <- rep("<.95", times = K)
gg$fit[gg$obs > gg$upper_95] <- ">.95"
gg$fit <- ordered(gg$fit, c(">.95", "<.95"))
poor_fit <- which(gg$fit == ">.95")

# Order
gg <- gg[order(gg$mean), ]
gg$group <- 1:nrow(gg)

# Plot
head(gg)
y_max <- max(gg$upper_99)

ggplot(gg[, ], aes(x = group, y = mean, group = fit)) +
    geom_errorbar(aes(ymin = lower_90, ymax = upper_90, color = fit), size = 2, width = 0) +
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95, color = fit), size = .5, width = 0, alpha = 1) +
    geom_point(aes(x = group, y = obs, pch = fit)) +
    scale_shape_manual(values = c(4, 20)) +
    scale_size_manual(values = c(4, 1)) +
    scale_color_manual(values = c("#132B43", "#56B1F7")) +
    #scale_color_manual(values = c("grey10", "grey70")) +
    xlab("Groups") +
    ylab("-2 * log-likelihood") +
    ylim(c(0, y_max)) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
          text = element_text(size=15),
          legend.title=element_blank())


# who doesnt fit?
collab$group_id[odd[poor_fit]]
apply(col_form[odd[poor_fit], ], 1, sum, na.rm = T)
apply(ind_form[odd[poor_fit], ], 1, sum, na.rm = T)
apply(ind_form[odd[poor_fit]+1, ], 1, sum, na.rm = T)

abs(theta1_logit[poor_fit, "mean"] - theta2_logit[poor_fit, "mean"]) %>% mean
abs(theta1_logit[-poor_fit, "mean"] - theta2_logit[-poor_fit, "mean"]) %>% mean
p(u_logit[poor_fit, "mean"])



hist(abs(theta1_logit[, "mean"] - theta2_logit[, "mean"] ))

# ------------------------------------------------------------
# Results plots
# ------------------------------------------------------------

# Quantiles to extract from stan fit
probs <- c(.005, .025, .05, .1, .9, .95, .975, .995)
int_names <- paste0(rep(c("lower", "upper"), each = 4), c("_99", "_95", "_90", "_80",  "_80",  "_90",  "_95", "_99"))

# Set up df
u_logit <- summary(fit_logit, pars = "u", probs = probs)$summary
gg <- data.frame(u_logit)


# Indicator for match on ability
gg$Ability <- "Within 1/2 SD"
gg$Ability[abs(theta1_logit[,"mean"]  - theta2_logit[,"mean"] ) > .5] <- "Not Within 1/2 SD"
gg$Ability <- ordered(gg$Ability, c("Within 1/2 SD", "Not Within 1/2 SD"))

# Rename intervals
names(gg)[grep("X", names(gg))] <- int_names
head(gg)

# Drop poor fitting groups ?
gg <- gg[-poor_fit,]

# Order
gg <- gg[order(gg$mean), ]
gg$group <- 1:nrow(gg)
# Plot
ggplot(gg[, ], aes(x = group, y = mean, group = Ability)) +
    geom_errorbar(aes(ymin = lower_80, ymax = upper_80, color = Ability), size = 2, width = 0) +
    geom_errorbar(aes(ymin = lower_90, ymax = upper_90, color = Ability), size = .5, width = 0, alpha = 1) +
    geom_point(aes(size = Ability, color = Ability)) +
    scale_color_manual(values = c("#132B43", "#56B1F7")) +
    #scale_color_manual(values = c("grey10", "grey70")) +
    scale_shape_discrete(solid = F) +
    ylab("Estimated RSC weight (logit)") +
    xlab("Group")  +
    geom_abline(slope = 0, intercept = 0, col = "black", size = 1) +
    scale_size_manual(values = c(2, 2)) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), text = element_text(size=15))

# prop strong synergy
table(gg$Ability)
(gg$lower_80[gg$Ability == "Within 1/2 SD"] > 0) %>% sum
15/46

(gg$lower_90[gg$Ability == "Within 1/2 SD"] > 0) %>% sum
10/46

# reliability
(var(gg[,"mean"]) - mean(gg[,"sd"]^2))/var(gg[,"mean"])
hist(gg[,"mean"])

### Results zoom in -------
ggplot(gg[-(1:100),], aes(x = group, y = mean, group = Ability)) +
    geom_errorbar(aes(ymin = lower_80, ymax = upper_80, color = Ability), size = 6, width = 0) +
    geom_errorbar(aes(ymin = lower_90, ymax = upper_90, color = Ability), size = .5, width = 0, alpha = 1) +
    geom_point(aes(size = Ability, color = Ability)) +
    scale_color_manual(values = c("#132B43", "#56B1F7")) +
    #scale_color_manual(values = c("grey10", "grey70")) +
    scale_shape_discrete(solid = F) +
    ylab("Estimated RSC weight (logit)") +
    xlab("Group")  +
    geom_abline(slope = 0, intercept = 0, col = "black", size = 1) +
    scale_size_manual(values = c(2, 2)) +
    #theme(base_size = 15) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), text = element_text(size=15))

### correlations with ability -------

ability <- c(apply(cbind(theta1_logit[-poor_fit,"mean"], theta2_logit[-poor_fit,"mean"]), 1, max),
apply(cbind(theta1_logit[-poor_fit,"mean"], theta2_logit[-poor_fit,"mean"]), 1, min),
theta1_logit[-poor_fit,"mean"] + theta2_logit[-poor_fit,"mean"],
abs(theta1_logit[-poor_fit,"mean"] - theta2_logit[-poor_fit,"mean"]))
RSC <- rep(u_logit[-poor_fit,"mean"], times = 4)
length(RSC)/4

temp <- c("max", "min", "mean", "|dif|")
Summary <- rep(temp, each = (length(RSC)/4))
Summary <- ordered(Summary, c("min", "max", "mean", "|dif|"))
gg <- data.frame(ability, RSC, Summary)
head(gg)

p <- ggplot(gg, aes(x = RSC, y = ability)) +
    geom_point(size = 2.5, col = "#56B1F7") +
    scale_shape_discrete(solid=F) +
    xlab("Estimated RSC weight (logit)") +
    ylab("Individual Ability")  +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), text = element_text(size=15))

p + facet_wrap( ~Summary, ncol = 2)


## sampling distributions of weights --------------
sample.int(K, 25)

u_poor_fit <- paste0("u[", poor_fit, "]")
u_random <- paste0("u[", sample.int(K, 25), "]")

stan_hist(fit_logit, pars = u_poor_fit, include = TRUE, col = "white", fill = "grey40") + xlim(-10, 10)

stan_hist(fit_logit, pars = u_random, include = TRUE, col = "white", fill = "grey40") + xlim(-10, 10)
