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
conditions <- c("Similar", "Disparate", "Random")
Conditions <- ordered(conditions, conditions)
n_models <- 4
n_obs <- 1000
n_items <- 100
n_reps <- 1
test_length <- 25
odd <- seq(1, n_obs*2, by = 2)
rep_order <- order(rep(1:n_obs, times = n_reps))

# Data generating parameters
alpha <- runif(n_items, .65, 2.5)
beta <- sort(rnorm(n_items, mean = .35, sd = 1.3))
parms <- data.frame(alpha, beta)
item_names <- paste0("item", 1:n_items)
row.names(parms) <- item_names
theta <- rnorm(n_obs*2)
theta_se <- SE(parms, theta)
mix_prop <- matrix(diag(rep(1, 4)), nrow = n_models, ncol = n_obs) %>% t
#mix_prop <- matrix(c(.02, .03, .05, .9), nrow = n_models, ncol = n_obs) %>% t


# ------------------------------------------------------------
# Simulated example A: Similar abililty
# ------------------------------------------------------------
ind <- order(theta)
a1 <- ind[odd]
a2 <- ind[odd + 1]
A <- data_gen(n_reps, mix_prop, parms, theta[a1], theta[a2])

em_a <- EM(models, A[item_names], parms, theta[a1], theta[a2])
em_a$prior

# ------------------------------------------------------------
# Simulated example B: Disparate ability
# ------------------------------------------------------------
b1 <- ind[1:n_obs]
b2 <- ind[(n_obs + 1):(n_obs * 2)]
B <- data_gen(n_reps, mix_prop, parms, theta[b1], theta[b2])
em_b <- EM(models, B[item_names], parms, theta[b1], theta[b2])
em_b$prior

# ------------------------------------------------------------
# Simulated example C: Randomly selected partners
# ------------------------------------------------------------
c1 <- odd
c2 <- odd + 1
C <- data_gen(n_reps, mix_prop, parms, theta[c1], theta[c2])
em_c <- EM(models, C[item_names], parms, theta[c1], theta[c2])
em_c$prior

# ------------------------------------------------------------
# Process loss
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
# ------------------------------------------------------------
deltas <- item_delta(parms, D$theta1, D$theta2)
delta_order <- apply(deltas, 1, order) %>% t

fun <- function(j, normed) {
  screen <- deltas*NA
  i <- j - test_length + 1
  ind <- cbind(rep(1:nrow(D), times = test_length), c(delta_order[,i:j]))
  screen[ind] <- 1
  redundancy(parms, D$theta1, D$theta2, normed = normed, NA_pattern = screen)
}

temp1 <- parallel::mclapply(test_length:(n_items-1), fun, normed = F)
temp2 <- parallel::mclapply(test_length:(n_items-1), fun, normed = T)
gg <- unlist(temp1) %>% data.frame
names(gg) <- "redun"
gg$pl <- 1- unlist(temp2)
gg$test <- as.ordered(rep(1:(n_items - test_length), each = nrow(D)))
gg$test2 <- as.ordered(rep(1:((n_items - test_length)/5)*5, each = nrow(D)*5))

gg$Condition <- rep(D$Condition, times = (n_items - test_length))
gg$Model <- D$Model

ggplot(gg, aes(x = test2, y = redun)) +
  geom_boxplot(outlier.shape = NA, size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(breaks=seq(0,75,5))

ggplot(gg, aes(x = test2, y = pl)) +
  geom_boxplot(outlier.shape = NA, size = .2, aes(fill = Condition)) +
  scale_fill_manual(values = c("grey30", "grey60", "white")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(breaks=seq(0,75,5))
