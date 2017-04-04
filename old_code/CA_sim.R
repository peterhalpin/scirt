# devtools::install_github("peterhalpin/cirt")
# library("cirt")

library("ggplot2")
library("dplyr")
source("~/github/cirt/R/cIRF_functions.R")
source("~/github/cirt/R/IRF_functions.R")
source("~/github/cirt/R/bootstrapping.R")

# ------------------------------------------------------------
# Data simulation
# ------------------------------------------------------------

# Constants

set.seed(101)
models <- c("Ind", "Min", "Max", "AI")
Models <- ordered(models, models)
n_models <- 4
n_obs <- 500
n_items <- 100
n_reps <- 50
test_length <- 25
odd <- seq(1, n_obs*2, by = 2)
rep_order <- order(rep(1:n_obs, times = n_reps))

# Data generating parameters
theta <- rnorm(n_obs*2)
beta <- sort(rnorm(n_items, mean = .35, sd = 1.3))
alpha <- runif(n_items, .65, 2.5)
parms <- data.frame(alpha, beta)
row.names(parms) <- paste0("item", 1:n_items)

# Determined classes
mix_prop_ones <- matrix(diag(rep(1, 4)), nrow = n_models, ncol = n_obs) %>% t

# load("gg")
# gg_a$sample <- sample.int(10, nrow(gg_a), replace = T)
# gg_b$sample <- sample.int(10, nrow(gg_b), replace = T)
# gg_c$sample <- sample.int(10, nrow(gg_c), replace = T)


# ------------------------------------------------------------
# Test non-redudancy
# ------------------------------------------------------------
theta1 <- theta2 <- seq(-5, 5, by = .1)
z <- outer(theta1, theta2, PL_Max, parms = parms)
z[upper.tri(z, diag = F)] <- NA
persp(theta2, theta1, z, theta = 65, phi = 15, expand = .5, ticktype = "detailed", nticks = 3, zlab = "", ylab = "", xlab = "", col = "grey")


# ------------------------------------------------------------
# Simulated example A: Randomly selected partners
# ------------------------------------------------------------

theta_1a <- theta[odd]
theta_2a <- theta[odd+1]
summary(theta_2a - theta_1a)

data_a <- data_gen(n_reps, mix_prop_ones, parms, theta_1a, theta_2a, fixed_class = T)

table(data_a$model[data_a$sample == 1])/ n_obs
head(data_a)
# Boostrapped EM estimates and SEs
fun <- function(i) {
  ind <- data_a$samples == i
  temp <- data_a[ind, grep("item", names(data_a))]
  EM(models, temp, parms, data_a$theta1[ind], data_a$theta2[ind], sorted = T)
}

temp_em <- parallel::mclapply(1:n_reps, fun)

# Save priors and se
out_a <- lapply(temp_em, function(x) c(x$prior, x$se^2)) %>% unlist %>% matrix(nrow = n_reps, ncol = 8, byrow = T) %>% data.frame

prior_a <- apply(out_a[,1:4], 2, mean) %>% round(3)
prior_se_a <- apply(out_a[,1:4], 2, sd) %>% round(3)

# Save posteriors
temp <- lapply(temp_em, function(x) x$posterior) %>% {do.call(rbind, .)}
data_a[models] <- temp[rep_order,]

# Expected log likelihoods
p_max_a <- matrix(c(0,0,1,0), ncol = 4, nrow = nrow(data_a), byrow = T)
e_max_a <- e_likelihood(p_max_a, parms, data_a$theta1, data_a$theta2, sorted = T)
e_mix_a <- e_likelihood(data_a[models], parms, data_a$theta1, data_a$theta2, sorted = T)


# Process loss
pl_max_a <- 1 - NPL(e_max_a, parms, data_a$theta1, data_a$theta2)
pl_mix_a <- 1 - NPL(e_mix_a, parms, data_a$theta1, data_a$theta2)
summary(pl_mix_a)



# Classification matrix and SEs
cp_a <- boot_cp(data_a, mix_prop[1,], parms)
post_a <- round(cp_a$mean, 3)
post_se_a <- round(cp_a$se, 3)

# Relationship bewteen item detal and model classification, for tests of fixed length
gg_a <- item_delta_gg(test_length, data_a, mix_prop[1,], parms, parallel = T)

pa <-
ggplot(gg_a, aes(x = delta, y = prob, group = model)) +
  geom_smooth(se = F, color = "black", method = "loess", aes(linetype = model)) +
  geom_point(data = gg_a[gg_a$sample == 1,], aes(pch = model), color = "black", alpha = .02) +
  scale_linetype_manual(values=c(2,4,6,1)) +
  scale_shape_manual(values=c(0,1,2,3), name = "") +
  ylim(c(0, 1)) +
  #xlab("Average scaled item delta") +
  #ylab("Probability of correct classification") +
  xlab("") +
  ylab("") +
  ggtitle("Random Ability") +
  #guides(line = guide_legend(order = 2)) +
  #guides(shape = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2))) +
  #theme(legend.spacing = unit(-.8, "cm")) +
  #theme(legend.key.size = unit(.9, "cm")) +
  #theme(legend.box = "horizontal") +
  theme(legend.position = "none")  +
  theme(plot.title = element_text(hjust = 0.5))

# ------------------------------------------------------------
# Simulated example B: high and low ability partners
# ------------------------------------------------------------

ind <- order(theta)
theta_1b <- theta[ind[1:n_obs]]
theta_2b <- theta[ind[(n_obs+1):(n_obs*2)]]
summary(theta_2a - theta_2b)

data_b <- data_gen(n_reps, mix_prop_ones, parms, theta_1b, theta_2b, fixed_class = T)
table(data_b$model[data_b$sample == 1])/ n_obs

# Boostrapped EM estimates and SEs
head(data_b)
# Boostrapped EM estimates and SEs
fun <- function(i) {
  ind <- data_b$samples == i
  temp <- data_b[ind, grep("item", names(data_b))]
  EM(models, temp, parms, data_b$theta1[ind], data_b$theta2[ind], sorted = T)
}

temp_em <- parallel::mclapply(1:n_reps, fun)

# Save priors and se
out_b <- lapply(temp_em, function(x) c(x$prior, x$se^2)) %>% unlist %>% matrix(nrow = n_reps, ncol = 8, byrow = T) %>% data.frame

prior_b <- apply(out_b[,1:4], 2, mean) %>% round(3)
prior_se_b <- apply(out_b[,1:4], 2, sd) %>% round(3)

# Save posteriors
temp <- lapply(temp_em, function(x) x$posterior) %>% {do.call(rbind, .)}
data_b[models] <- temp[rep_order,]


# Expected log likelihoods
p_max_b <- matrix(c(0,0,1,0), ncol = 4, nrow = nrow(data_b), byrow = T)
e_max_b <- e_likelihood(p_max_b, parms, data_b$theta1, data_b$theta2, sorted = T)
e_mix_b <- e_likelihood(data_b[models], parms, data_b$theta1, data_b$theta2, sorted = T)

# Process loss
pl_max_b <- 1 - NPL(e_max_b, parms, data_b$theta1, data_b$theta2)
pl_mix_b <- 1 - NPL(e_mix_b, parms, data_b$theta1, data_b$theta2)
max(pl_mix_b)
# Classification matrix and SEs
cp_b <- boot_cp(data_b, mix_prop[1,], parms)
post_b <- round(cp_b$mean, 3)
post_se_b <- round(cp_b$se, 3)




# Relationship bewteen item detal and model classification, for tests of fixed length
load("gg")
gg_b <- item_delta_gg(test_length, data_b, mix_prop[1,], parms, parallel = T)

gg_b$sample <- sample.int(10, nrow(gg_b), replace = T)


pb <-
ggplot(gg_b, aes(x = delta, y = prob, group = model)) +
  geom_smooth(se = F, color = "black", method = "loess", aes(linetype = model)) +
  geom_point(data = gg_b[gg_b$sample == 1,], aes(pch = model), color = "black", alpha = .015) +
  scale_linetype_manual(values=c(2,4,6,1)) +
  scale_shape_manual(values=c(0,1,2,3), name = "") +
  ylim(c(0, 1)) +
  xlab("Average scaled item delta") +
  #ylab("Probability of correct classification") +
  ylab("") +
  ggtitle("Disparate Ability") +
  guides(line = guide_legend(order = 2)) +
  guides(shape = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2))) +
  theme(legend.spacing = unit(-.8, "cm")) +
  theme(legend.key.size = unit(.9, "cm")) +
  theme(legend.box = "horizontal") +
  theme(legend.position = c(.5, .15)) +
  theme(plot.title = element_text(hjust = 0.5))




# ------------------------------------------------------------
# Simulated example C: matching partners close in ability
# ------------------------------------------------------------
ind <- order(theta)
theta_1c <- theta[ind[odd]]
theta_2c <- theta[ind[odd+1]]
summary(theta_2c - theta_1c)

data_c <- data_gen(n_reps, mix_prop_ones, parms, theta_1c, theta_2c, fixed_class = T)
table(data_c$model[data_c$sample == 1])/ n_obs

# Boostrapped EM estimates and SEs
fun <- function(i) {
  ind <- data_c$samples == i
  temp <- data_c[ind, grep("item", names(data_c))]
  EM(models, temp, parms, data_c$theta1[ind], data_c$theta2[ind], sorted = T)
}

temp_em <- parallel::mclapply(1:n_reps, fun)

# Save priors and se
out_b <- lapply(temp_em, function(x) c(x$prior, x$se^2)) %>% unlist %>% matrix(nrow = n_reps, ncol = 8, byrow = T) %>% data.frame

prior_c <- apply(out_b[,1:4], 2, mean) %>% round(3)
prior_se_c <- apply(out_b[,1:4], 2, sd) %>% round(3)

# Save posteriors
temp <- lapply(temp_em, function(x) x$posterior) %>% {do.call(rbind, .)}
data_c[models] <- temp[rep_order,]



# ------------------------------------------------------------
# Process loss
# ------------------------------------------------------------

# Expected log likelihoods
p_max_c <- matrix(c(0,0,1,0), ncol = 4, nrow = nrow(data_c), byrow = T)
e_max_c <- e_likelihood(p_max_c, parms, data_c$theta1, data_c$theta2, sorted = T)
e_mix_c <- e_likelihood(data_c[models], parms, data_c$theta1, data_c$theta2, sorted = T)

# Process loss
pl_max_c <- 1 - NPL(e_max_c, parms, data_c$theta1, data_c$theta2)
pl_mix_c <- 1 - NPL(e_mix_c, parms, data_c$theta1, data_c$theta2)
max(pl_mix_c)


# Classification matrix and SEs
cp_c <- boot_cp(data_c, mix_prop[1,], parms)
post_c <- round(cp_c$mean, 3)
post_se_c <- round(cp_c$se, 3)


# Relationship bewteen item detal and model classification, for tests of fixed length
gg_c <- item_delta_gg(test_length, data_c, mix_prop[1,], parms, parallel = T)

pc <-
ggplot(gg_c, aes(x = delta, y = prob, group = model)) +
  geom_smooth(se = F, color = "black", method = "loess", aes(linetype = model)) +
  geom_point(data = gg_c[gg_c$sample == 1,], aes(pch = model), color = "black", alpha = .02) +
  scale_linetype_manual(values=c(2,4,6,1)) +
  scale_shape_manual(values=c(0,1,2,3), name = "") +
  ylim(c(0, 1)) +
  xlab("") +
  #xlab("Average scaled item delta") +
  ylab("Probability of correct classification") +
  ggtitle("Similar Ability") +
  #guides(line = guide_legend(order = 2)) +
  #guides(shape = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2))) +
  #theme(legend.spacing = unit(-.8, "cm")) +
  #theme(legend.key.size = unit(.9, "cm")) +
  #theme(legend.box = "horizontal") +
  theme(legend.position = "none")  +
  theme(plot.title = element_text(hjust = 0.5))



# ------------------------------------------------------------
# Tables and Figures
# ------------------------------------------------------------

# Priors
xtable::xtable(rbind(
    paste0(prior_c, " (", prior_se_c, ") "),
    paste0(prior_b, " (", prior_se_b, ") "),
    paste0(prior_a, " (", prior_se_a, ") ")))

# Classification probabilities
post_a <- round(cp_a$mean, 3)
post_se_a <- round(cp_a$se, 3)
xtable::xtable(rbind(
    paste0(post_c, " (", post_se_c, ") ") %>% array(, dim = c(4,4)),
    paste0(post_b, " (", post_se_b, ") ") %>% array(, dim = c(4,4)),
    paste0(post_a, " (", post_se_a, ") ") %>% array(, dim = c(4,4))))

gridExtra::grid.arrange(pc, pb, pa, nrow = 1, ncol = 3)


# process loss

gg <- data.frame(c(tapply(pl_mix_a, data_a$pairs, mean), tapply(pl_mix_b, data_b$pairs, mean), tapply(pl_mix_c, data_c$pairs, mean)))
names(gg) <- "pl"
gg$model <- rep(Models[data_a$model[!duplicated(data_a$pairs)]], times = 3)
gg$condition <- rep(c("Random", "Disparate", "Similar"), each = n_obs)
gg$condition <- ordered(gg$condition, c("Similar", "Disparate", "Random"))
gg$int <- interaction(gg$condition, gg$model)
head(gg)

ggplot(gg, aes(x = int, y = pl)) +
  geom_boxplot(outlier.size = 0, size = .2, aes(fill = condition)) +
  scale_fill_manual(values = c("grey90", "grey60", "grey30")) +
  ylab("1 - Process Loss") +
  xlab("Condition by Model")

gg <- data.frame(c(tapply(e_mix_a, data_a$pairs, mean), tapply(e_mix_b, data_b$pairs, mean), tapply(e_mix_c, data_c$pairs, mean)))
names(gg) <- "pl"
gg$model <- rep(Models[data_a$model[!duplicated(data_a$pairs)]], times = 3)
gg$condition <- rep(c("Random", "Disparate", "Similar"), each = n_obs)
gg$condition <- ordered(gg$condition, c("Similar", "Disparate", "Random"))
gg$int <- interaction(gg$condition, gg$model)
head(gg)

ggplot(gg, aes(x = int, y = pl)) +
  geom_boxplot(outlier.size = 0, size = .2, aes(fill = condition)) +
  scale_fill_manual(values = c("grey90", "grey60", "grey30")) +
  ylab("1 - Process Loss") +
  xlab("Condition by Model")
