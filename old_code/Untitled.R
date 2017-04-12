# devtools::install_github("peterhalpin/cirt")
# library("cirt")

library("ggplot2")
library("dplyr")
source("~/github/cirt/R/cIRF_functions.R")
source("~/github/cirt/R/IRF_functions.R")
source("~/github/cirt/R/bootstrapping.R")

# ------------------------------------------------------------
# dataata simulation
# ------------------------------------------------------------

# Constants
set.seed(101)
models <- c("Ind", "Min", "Max", "AI")
Models <- ordered(models, models)
conditions <- c("Similar", "Disparate", "Random")
Conditions <- ordered(conditions, conditions)

# Data generating parameters
n_models <- length(models)
n_obs <- 500
n_items <- 100
n_reps <- 500
alpha <- runif(n_items, .65, 2.5)
beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
parms <- data.frame(alpha, beta)
theta <- rnorm(n_obs*2)
theta_se <- SE(parms[sample.int(n_items, 25), ], theta)
mix_prop <- matrix(diag(rep(1, 4)), nrow = n_models, ncol = n_obs) %>% t

# Other stuff
item_names <- paste0("item", 1:n_items)
row.names(parms) <- item_names
odd <- seq(1, n_obs*2, by = 2)
rep_order <- order(rep(1:n_obs, times = n_reps))
test_length <- 25


# ------------------------------------------------------------
# Simulated example A: Similar abililty
# ------------------------------------------------------------
ind <- order(theta)
a1 <- ind[odd]
a2 <- ind[odd + 1]

A <- data_gen(1, mix_prop, parms, theta[a1], theta[a2])
em_a <- EM(models, A[item_names], parms, theta[a1], theta[a2])

pv_a <- pv_gen(n_reps, A[item_names], parms, theta[a1], theta[a2], theta_se[a1], theta_se[a2], model = A$model)
em_pv_a <- parallel::mclapply(1:n_reps, em_parallel, sim_data = pv_a, parms = parms)

# ------------------------------------------------------------
# Simulated example B: Disparate ability
# ------------------------------------------------------------
b1 <- ind[1:n_obs]
b2 <- ind[(n_obs + 1):(n_obs * 2)]
B <- data_gen(1, mix_prop, parms, theta[b1], theta[b2])
em_b <- EM(models, B[item_names], parms, theta[b1], theta[b2])

pv_b <- pv_gen(n_reps, B[item_names], parms, theta[b1], theta[b2], theta_se[b1], theta_se[b2], model = B$model)
em_pv_b <- parallel::mclapply(1:n_reps, em_parallel, sim_data = pv_b, parms = parms)

# ------------------------------------------------------------
# Simulated example C: Randomly selected partners
# ------------------------------------------------------------
c1 <- odd
c2 <- odd + 1
C <- data_gen(1, mix_prop, parms, theta[c1], theta[c2])
em_c <- EM(models, C[item_names], parms, theta[c1], theta[c2])

pv_c <- pv_gen(n_reps, C[item_names], parms, theta[c1], theta[c2], theta_se[c1], theta_se[c2], model = C$model)
em_pv_c <- parallel::mclapply(1:n_reps, em_parallel, sim_data = pv_c, parms = parms)

# ------------------------------------------------------------
# Parameter recovery
# ------------------------------------------------------------
em_priors_a <-  rbind(em_a$prior, em_a$se) %>% round(3)
em_priors_b <-  rbind(em_b$prior, em_b$se) %>% round(3)
em_priors_c <-  rbind(em_c$prior, em_c$se) %>% round(3)

pv_priors_a <- pv_priors(em_pv_a) %>% round(3)
pv_priors_b <- pv_priors(em_pv_b) %>% round(3)
pv_priors_c <- pv_priors(em_pv_c) %>% round(3)

xtable::xtable(rbind(
  paste0(em_priors_a[1,], " (", em_priors_a[2,], ")"),
  paste0(pv_priors_a[1,], " (", pv_priors_a[4,], ")"),
  paste0(pv_priors_a[5,]),
  paste0(em_priors_b[1,], " (", em_priors_b[2,], ")"),
  paste0(pv_priors_b[1,], " (", pv_priors_b[4,], ")"),
  paste0(pv_priors_b[5,]),
  paste0(em_priors_c[1,], " (", em_priors_c[2,], ")"),
  paste0(pv_priors_c[1,], " (", pv_priors_c[4,], ")"),
  paste0(pv_priors_c[5,])))

# ------------------------------------------------------------
# Process loss with pv
# ------------------------------------------------------------
pv_a[models] <- pv_posteriors(em_pv_a, rep_order)
pv_b[models] <- pv_posteriors(em_pv_b, rep_order)
pv_c[models] <- pv_posteriors(em_pv_c, rep_order)

pv_a$e_mix <- e_mix(pv_a[models], parms, pv_a$theta1, pv_a$theta2, sorted = T)
pv_b$e_mix <- e_mix(pv_b[models], parms, pv_b$theta1, pv_b$theta2, sorted = T)
pv_c$e_mix <- e_mix(pv_c[models], parms, pv_c$theta1, pv_c$theta2, sorted = T)

pv <- rbind(pv_a, pv_b, pv_c)
pv$pl <-  1 - PL(pv$e_mix, parms, pv$theta1, pv$theta2, normed = T)
pv$pl2 <-  PL(pv$e_mix, parms, pv$theta1, pv$theta2, normed = F)
pv$Condition <- rep(Conditions, each = n_obs*n_reps)
pv$Model <- Models[pv$model]
pv$int <- interaction(pv$Condition, pv$Model)

ggplot(pv, aes(x = int, y = pl)) +
  geom_boxplot(outlier.shape = NA, outlier.color = "white", size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  ylab("1 - Process Loss") +
  xlab("Condition by Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(pv, aes(x = int, y = pl2)) +
  geom_boxplot(outlier.shape = NA, size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  ylab("Process Loss") +
  xlab("Condition by Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(c(0,120))

# ------------------------------------------------------------
# Process loss with original data
# ------------------------------------------------------------

A$e_mix <- e_mix(em_a$posterior, parms, theta[a1], theta[a2])
B$e_mix <- e_mix(em_b$posterior, parms, theta[b1], theta[b2])
C$e_mix <- e_mix(em_c$posterior, parms, theta[c1], theta[c2])

D <- rbind(A, B, C)
D$pl <-  1 - PL(D$e_mix, parms, D$theta1, D$theta2, normed = T)
D$pl2 <-  PL(D$e_mix, parms, D$theta1, D$theta2, normed = F)
D$Condition <- rep(Conditions, each = n_obs)
D$Model <- Models[D$model]
D$int <- interaction(D$Condition, D$Model)

ggplot(D, aes(x = int, y = pl)) +
  geom_boxplot(outlier.shape = NA, outlier.color = "white", size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  #  scale_fill_manual(values = c("grey90", "grey60", "grey30")) +
  ylab("1 - Process Loss") +
  xlab("Condition by Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(D, aes(x = int, y = pl2)) +
  geom_boxplot(outlier.shape = NA, size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  ylab("Process Loss") +
  xlab("Condition by Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(c(0,120))

# ------------------------------------------------------------
# Test redundancy
# Use pv or orginal data?
  # data <- D
    data <- pv
# ------------------------------------------------------------

deltas <- item_delta(parms, data$theta1, data$theta2)
delta_order <- apply(deltas, 1, order) %>% t
n_tests <- n_items - test_length

fun <- function(j, normed) {
  screen <- deltas*NA
  i <- j - test_length + 1
  ind <- cbind(rep(1:nrow(data), times = test_length), c(delta_order[,i:j]))
  screen[ind] <- 1
  redundancy(parms, data$theta1, data$theta2, normed = normed, NA_pattern = screen)
}

temp1 <- parallel::mclapply(test_length:(n_items-1), fun, normed = F)
temp2 <- parallel::mclapply(test_length:(n_items-1), fun, normed = T)
gg <- unlist(temp1) %>% data.frame
names(gg) <- "redun"
gg$pl <- 1- unlist(temp2)
gg$test <- as.ordered(rep(1:n_tests, each = nrow(data)))
gg$test2 <- as.ordered(rep(1:(n_tests/5)*5, each = nrow(data)*5))
gg$Condition <- rep(data$Condition, times = n_tests)
gg$Model <- data$Model

ggplot(gg, aes(x = test2, y = redun)) +
  geom_boxplot(outlier.shape = NA, size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(breaks=seq(0,75,5)) +
  ylab("Process Loss") +
  xlab("Test")

ggplot(gg, aes(x = test2, y = pl)) +
  geom_boxplot(outlier.shape = NA, size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(breaks=seq(0,75,5)) +
  ylab("1 - Process Loss") +
  xlab("Test")
