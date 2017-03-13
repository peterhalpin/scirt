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
n_reps <- 500
test_length <- 25
odd <- seq(1, n_obs*2, by = 2)

# Data generating parameters
theta <- rnorm(n_obs*2)
beta <- sort(rnorm(n_items, mean = .35, sd = 1.3))
alpha <- runif(n_items, 1, 2.5)
parms <- data.frame(alpha, beta)
row.names(parms) <- paste0("item", 1:n_items)

# Determined classes
mix_prop_ones <- matrix(diag(rep(1, 4)), nrow = n_models, ncol = n_obs) %>% t

# load("gg")
# gg_a$sample <- sample.int(10, nrow(gg_a), replace = T)
# gg_b$sample <- sample.int(10, nrow(gg_b), replace = T)
# gg_c$sample <- sample.int(10, nrow(gg_c), replace = T)

# ------------------------------------------------------------
# Simulated example A: Randomly selected partners
# ------------------------------------------------------------

theta_1a <- theta[odd]
theta_2a <- theta[odd+1]
dif_a <- theta_2a - theta_1a

data_a <- data_gen(n_reps, mix_prop_ones, parms, theta_1a, theta_2a, fixed_class = T)
table(data_a$model[data_a$sample == 1])/ n_obs

# Boostrapped EM estimates and SEs
em_a <- boot_em(data_a, parms, parallel = T)
prior_a <- apply(em_a[,1:4], 2, mean) %>% round(3)
prior_se_a <- apply(em_a[,1:4], 2, sd) %>% round(3)

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
em_b <- boot_em(data_b, parms, parallel = T)
prior_b <- apply(em_b[,1:4], 2, mean) %>% round(3)
prior_se_b <- apply(em_b[,1:4], 2, sd) %>% round(3)

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
em_c <- boot_em(data_c, parms, parallel = T)
prior_c <- apply(em_c[,1:4], 2, mean) %>% round(3)
prior_se_c <- apply(em_c[,1:4], 2, sd) %>% round(3)

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


# ------------------------------------------------------------
# Real data example: demographics
#------------------------------------------------------------

setwd("~/Dropbox/Academic/Projects/CA/Data")
demo <- read.csv("Long_AMT_IDs_Match(1.30.17).csv")
names(demo)
vars <- c("Age", "Gender", "Ethnicity", "English", "leveledu", "liveusa")

temp <- demo[vars]
temp$Ethnicity[temp$Ethnicity != 10] <- 0
temp$Ethnicity[temp$Ethnicity != 0] <- 1
temp$Gender <- temp$Gender - 1
temp$English <- abs(temp$English - 2)
temp$leveledu[temp$leveledu < 3] <- 0
temp$leveledu[temp$leveledu != 0] <- 1
temp$liveusa <- abs(temp$liveusa - 2)

head(temp)

apply(temp, 2, mean)
summary(temp$Age)



#------------------------------------------------------------
# DIF: Items 45 and 65 identified as problematic in scalar model
# After dropping, scalar model fits OK

# Metric            36.971        37       0.4704
# Scalar            101.071       74       0.0200
# Scalar (w/ drop)  84.358        70       0.1161
# ------------------------------------------------------------

 # Load calibrated item parms
setwd("~/Dropbox/Academic/Projects/CA/Data/response_matrices")
parms <- read.csv("calibration_parms.csv", row.names = 1)
grep_items <- paste0(row.names(parms), collapse = "|")
scalar.out <- "./DIF/collaboration_DIF_scalar.out"

temp <- readLines(scalar.out)
Begin <- grep("Item Difficulties", temp)
End <- grep("Variances", temp)

calib_temp <- temp[Begin[1]:End[End > Begin[1] & End < Begin[2]]]
collab_temp <- temp[Begin[2]:End[End > Begin[2] & End < Begin[3]][1]]

calib_beta <- substr(calib_temp[grepl(grep_items, calib_temp)], 23, 28) %>% as.numeric
collab_beta <- substr(collab_temp[grepl(grep_items, collab_temp)], 23, 28) %>% as.numeric
item_names <- substr(collab_temp[grepl(grep_items, collab_temp)], 6, 10)

gg <- data.frame(item_names, collab_beta, calib_beta)
gg$dif <- resid(lm(collab_beta ~ calib_beta, gg)) > 1
dif_items <- paste0(gg$item_names[gg$dif == T], collapse = "|")

# DIF Figure
ggplot(gg, aes(x = calib_beta, y = collab_beta)) +
  geom_point(aes(pch = dif, color = dif, size = dif)) +
  stat_smooth(method = "lm", col = "black", se = F) +
  xlab("Calibration sample") +
  ylab("Collaborative testing condition") +
  theme(legend.position = "none") +
  scale_shape_manual(values = c(20, 4)) +
  scale_color_manual(values = c("black", "red")) +
  scale_size_manual(values = c(2, 4))


# ------------------------------------------------------------
# Real data example: load data and estimate individual thetas
# ------------------------------------------------------------

# Load calibrated item parms
setwd("~/Dropbox/Academic/Projects/CA/Data/response_matrices")
parms <- read.csv("calibration_parms.csv", row.names = 1)
summary(parms)

# Drop DIF items
dif_items <- "045|065"
parms <- parms[!grepl(dif_items, row.names(parms)),]

# Load collaboration data and split into forms
collab <- read.csv("collaboration_2016.csv", check.names = F)
col_form <- format_resp(collab, row.names(parms), "COL")
ind_form <- format_resp(collab, row.names(parms), "IND")

 # Apply conjunctive scoring rule
odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]

# Drop 13 unsuable response patterns (all 1 or all 0)
drop_groups <- c(
  collab$group_id[apply(col_form, 1, mean, na.rm = T) %in% c(1,0)],
  collab$group_id[apply(ind_form, 1, mean, na.rm = T) %in% c(1,0)])

col_form <-col_form[!collab$group_id%in%drop_groups,]
ind_form <-ind_form[!collab$group_id%in%drop_groups,]

# Reset odd for dropped items
odd <- odd[1:(nrow(col_form)/2)]

# Estimate theta for ind forms
ind <- MLE(ind_form, parms, WMLE = T)
plot(ind$theta, ind$se)

theta1 <- ind$theta[odd]
theta2 <- ind$theta[odd+1]
theta1_se <- ind$se[odd]
theta2_se <- ind$se[odd+1]

resp <- col_form[odd, ]

# ------------------------------------------------------------
# real data: EM
# ------------------------------------------------------------
em <- EM(models, resp, parms, theta1, theta2)

# mixing proportions
round(em$prior, 3)
round(em$se, 3)

# classification probabilities
classify <- class_probs(em$posterior)
round(classify, 3)

# ------------------------------------------------------------
# real data: person fit using contribution to likelihood
# ------------------------------------------------------------

# Get observed values
components <- likelihood(models, resp, parms, theta1, theta2, Log = F)
logL <- incomplete_data(components, em$posterior, Sum = F)

# Simulate null distribution
n_reps <- 20
mix_prop <- em$posterior

temp <- data_gen(n_reps, mix_prop, parms, theta1, theta2, theta1_se, theta2_se, NA_pattern = resp)

temp_components <- likelihood(models, temp[,grep(grep_items, names(temp))], parms, temp$theta1, temp$theta2, sorted = T, Log = F)

temp_logL <- incomplete_data(temp_components, temp[models], Sum = F)

# Get quantiles
quant <- tapply(temp_logL, temp$pair, function(x) quantile(x, p = c(.5, .025))) %>% unlist %>% matrix(, nrow = 2, ncol = length(theta1)) %>% t()

# Visual key for misfit
fit <- rep("<.95", times = length(theta1))
fit[logL < quant[,2]] <- ">.95"
fit <- ordered(fit, c(">.95", "<.95"))

# Set up and plot
gg <- data.frame(-2*temp_logL, -2*rep(logL, each = n_reps), rep(fit, each = n_reps), rep(-2*quant[,1], each = n_reps))

names(gg) <- c("l_dist", "l_obs", "fit", "median")
gg <- gg[order(gg$median), ]
gg$pair <- rep(1:length(theta1), each = n_reps)

ggplot(gg, aes(x = pair, y = l_dist, group = pair)) +
  geom_boxplot(outlier.size = 0, outlier.color = "white", aes(fill = fit)) +
  geom_point(aes(x = pair, y = l_obs, pch = fit, size = fit)) +
  scale_shape_manual(values = c(4, 20)) +
  scale_fill_grey(start = 0.1, end = 0.8) +
  xlab("Groups") +
  ylab("-2 * loglikelihood") +
  scale_size_manual(values = c(4, 1)) +
  theme(legend.title=element_blank())


# ------------------------------------------------------------
# real data:  Plausible Values
# ------------------------------------------------------------
n_reps <- 500
set.seed(101)
pv_data <- pv_gen(n_reps, resp, parms, theta1, theta2, theta1_se, theta2_se)
pv_data[models] = 0
head(pv_data)
out <- data.frame(matrix(0, nrow = n_reps, ncol <- n_models*2))
names(out) <- paste0(rep(c("prior", "se"), each = n_models), 1:n_models)

for(i in 1:n_reps) {
  ind <- pv_data$samples == i
  temp <- pv_data[ind, grep(paste0(row.names(parms), collapse = "|"), names(pv_data))]
  temp_em <- EM(models, temp, parms, pv_data$theta1[ind], pv_data$theta2[ind], sorted = T)
  out[i, ] <- c(temp_em$prior, temp_em$se^2)
  pv_data[ind, models] <- temp_em$posterior
}


# ------------------------------------------------------------
# real data:  priors
# ------------------------------------------------------------

mean_out <- apply(out, 2, mean, na.rm = T)
var_out <- apply(out, 2, var, na.rm = T)
pv_prior <- round(mean_out[1:4], 3)
pv_se <- sqrt(mean_out[5:8] + (1 + 1/n_reps) * var_out[1:4])
pvl <- (1 + 1/n_reps) * var_out[1:4]/mean_out[5:8]

temp <- rbind(em$prior, em$se, pv_prior, pv_se, pvl)
row.names(temp) <- c("em_prior", "em_se", "pv_prior", "pv_se", "pvl")
colnames(temp) <- models
xtable::xtable(temp, digits = 3)
round(temp, 4)



# ------------------------------------------------------------
# real data: posteriors
# ------------------------------------------------------------

# Get the "true model" from PVs

temp <- lapply(pv_data[models], function(x) tapply(x, pv_data$pairs, mean)) %>% data.frame
pv_data$model <- rep(apply(temp, 1, which.max), each = n_reps)

cp  <- boot_cp(pv_data, pv_prior, parms)
post <- round(cp$mean, 3)
post_se <- round(cp$se, 3)

xtable::xtable(paste0(post, " (", post_se, ") ") %>% array(, dim = c(4,4)))


# ------------------------------------------------------------
# real data: classification probabilities
# ------------------------------------------------------------

gg <- data.frame(unlist(pv_data[models]))
names(gg) <- "prob"
gg$q <- rep(as.matrix(pv_data[models]) %*% 1:n_models, times = n_models)
gg$model <- rep(Models, each = nrow(pv_data))
gg$sample <- sample.int(500, nrow(gg), replace = T)
sum(gg$sample == 1)

p1 <-
ggplot(gg, aes(x = q, y = prob, group = model)) +
  stat_smooth(aes(linetype = model), se = F, lwd = 1, color = "black") +
  geom_point(data = gg[gg$sample == 1,], aes(pch = model), color = "black", alpha = .5) +
  xlab("Expectation of posterior distribtuion") +
  ylab("Posterior probability of each model") +
  scale_linetype_manual(values=c(2,4,6,1)) +
  scale_shape_manual(values=c(0,1,2,3), name = "") +
  guides(line = guide_legend(order = 2)) +
  guides(shape = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2))) +
  guides(fill = FALSE) +
  theme(legend.spacing = unit(-.8, "cm")) +
  theme(legend.key.size = unit(.9, "cm")) +
  theme(legend.box = "horizontal")+
  ylim(c(0, 1))



# ------------------------------------------------------------
# real data: standard error of posterior
# ------------------------------------------------------------

temp <- as.matrix(pv_data[models]) %*% 1:4
pv_q <- tapply(temp, pv_data$pairs, mean)
pv_se <- tapply(temp, pv_data$pairs, sd)
model <- lapply(pv_data[models], function(x) tapply(x, pv_data$pairs, mean)) %>% data.frame() %>% apply(1, which.max)
gg$model <- Models[gg$model]
gg <- data.frame(pv_q, pv_se, model)
names(gg) <- c("pv", "se", "model")


p2 <-
ggplot(gg, aes(x = pv, y = se, group = model)) +
  stat_smooth(aes(group = 1), se = F, lwd = 1, color = "black") +
  geom_point(aes(pch = model), color = "black", alpha = .5) +
  scale_linetype_manual(values=c(2,4,6,1)) +
  scale_shape_manual(values=c(0,1,2,3), name = "") +
  guides(line = guide_legend(order = 2)) +
  guides(shape = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2))) +
  guides(fill = FALSE) +
  theme(legend.spacing = unit(-.8, "cm")) +
  theme(legend.key.size = unit(.9, "cm")) +
  theme(legend.box = "horizontal")+
  xlab("Expectation of posterior") +
  ylab("Standard error")
  # ggtitle("Standard error of expectation of posterior distribution")


tapply(gg$pv, gg$model, mean)
tapply(gg$pv, gg$model, sd)
gridExtra::grid.arrange(p1, p2, nrow = 1)

# ------------------------------------------------------------
# AI / Max by theta deltas
# ------------------------------------------------------------



gg <- lapply(pv_data[models], function(x) tapply(x, pv_data$pairs, mean)) %>% data.frame()

temp <- item_delta(parms, pv_data$theta1, pv_data$theta2, sorted = T, NA_pattern = resp)/.25
temp[temp > 1] <- 1
temp <- apply(temp, 1, mean, na.rm = T)
gg$delta <-  tapply(temp, pv_data$pairs, mean)
head(gg)

gg$model <- Models[apply(gg[models], 1, which.max)]
ggplot(gg, aes(x = delta, fill = model)) + geom_histogram(color = "black") + scale_fill_grey(start = 0.9, end = 0.2, na.value = "grey50")

plot(gg$delta, gg$AI)
hist(gg$delta)

temp <- item_delta(parms, pv_data$theta1, pv_data$theta2, sorted = T, NA_pattern = resp)

head(temp)

gg$delta <-  apply(temp, 1, mean, na.rm = T) %>% rep(times = n_models)

gg$model <- rep(Models, each = nrow(pv_data))
gg$sample <- sample.int(500, nrow(gg), replace = T)
head(gg)

ggplot(gg[model == "AI"], aes(x = delta, y = prob, group = model)) +
  stat_smooth(aes(linetype = model), se = F, lwd = 1, color = "black") +
  geom_point(aes(pch = model), color = "black", alpha = .01) +
  scale_linetype_manual(values=c(2,4,6,1)) +
  scale_shape_manual(values=c(0,1,2,3), name = "") +
  guides(line = guide_legend(order = 2)) +
  guides(shape = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2))) +
  guides(fill = FALSE) +
  theme(legend.spacing = unit(-.8, "cm")) +
  theme(legend.key.size = unit(.9, "cm")) +
  theme(legend.box = "horizontal")

  +
  xlab("Expectation of posterior") +
  ylab("Standard error")
  # ggtitle("Standard error of expectation of posterior distribution")


tapply(gg$pv, gg$model, mean)
tapply(gg$pv, gg$model, sd)


gg <- data.frame("model" = apply(em$posterior, 1, which.max))
gg$t_delta <- abs(theta1 - theta2)
gg$i_delta <- apply(deltas(parms, theta1, theta2)/.25, 1 mean)
gg[models] <- em$posterior
head(gg)

plot(gg$delta, gg$AI)


# ------------------------------------------------------------
# CP versus item deltas
# ------------------------------------------------------------
deltas <- item_delta(parms, theta1, theta2)/.25
delta_order <- t(apply(deltas, 1, order))
round(deltas[1,delta_order[1,]], 4)
screen <- deltas*NA
n_short <- 50
ind <- cbind(1:4, c(2,3,2,3))
ind <- cbind(1:4, 1:4)
out <- data.frame(rbind((classify)[ind]), sum(deltas), n_items)
names(out) <- c(models, "deltas", "n_items")
i = 50
for(i in n_short:n_items) {
  screen <- deltas*NA
  screen[delta_order < (i+1) & delta_order > (i - n_short)] <- 1
  apply(screen, 1, sum, na.rm = T)
  components <- likelihood(models, resp*screen, parms, theta1, theta2, Log = F)
  temp_post <- posterior(components, c(sample_mix_prop))
  #temp_post <- posterior(components, em$prior)
  out[(i+1), 1:4] <- class_probs(temp_post, data$model)[ind]
  out$deltas[i+1] <- sum(deltas*screen, na.rm = T)
  #out$n_items[i+1] <- n_items - i
  cat(i)
}

# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# Old simulation stuff
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------


# Compute standard errors on theta using a half of the items
ind_beta <- sort(rnorm(n_ind_items, mean = .35, sd = 1.3))
ind_alpha <- runif(n_ind_items, .7, 2.5)
ind_parms <- data.frame(ind_alpha, ind_beta)
names(ind_parms) <- c("alpha", "beta")
theta_se <- SE(ind_parms, theta)

# Pair odd and even elements of theta / theta_se
odd <- seq(1, n_obs*2, by = 2)
theta1 <- theta[odd]
theta2 <- theta[odd+1]
theta1_se <- theta_se[odd]
theta2_se <- theta_se[odd+1]

# Pair high and low
# ind <- order(theta)
# theta1 <- theta[ind[1:n_obs]]
# theta2 <- theta[ind[(n_obs+1):(n_obs*2)]]
# theta1_se <- theta_se[ind[1:n_obs]]
# theta2_se <- theta_se[ind[(n_obs+1):(n_obs*2)]]
#
# # Pair neighbors
# ind <- order(theta)
# theta1 <- theta[ind[odd]]
# theta2 <- theta[ind[odd+1]]
# theta1_se <- theta_se[ind[odd]]
# theta2_se <- theta_se[ind[odd+1]]

# Generate data
data <- data_gen(1, mix_prop, parms, theta1, theta2, expected = F)
sample_mix_prop <- table(data$model) / n_obs
head(data)


# ------------------------------------------------------------
# EM
# ------------------------------------------------------------
resp <- data[grep("item", names(data))]
em <- EM(models, resp, parms, theta1, theta2)

# mixing proportions
round(em$prior, 3)
round(em$se, 3)

# classification probabilities
classify <- class_probs(em$posterior, data$model)
#classify <- class_probs(em$posterior)

round(classify, 3)



# ------------------------------------------------------------
# CP versus item deltas
# ------------------------------------------------------------
deltas <- item_delta(parms, theta1, theta2)/.25
delta_order <- t(apply(deltas, 1, order))
round(deltas[1,delta_order[1,]], 4)
screen <- deltas*NA
n_short <- 50
ind <- cbind(1:4, c(2,3,2,3))
ind <- cbind(1:4, 1:4)
out <- data.frame(rbind((classify)[ind]), sum(deltas), n_items)
names(out) <- c(models, "deltas", "n_items")
i = 50
for(i in n_short:n_items) {
  screen <- deltas*NA
  screen[delta_order < (i+1) & delta_order > (i - n_short)] <- 1
  apply(screen, 1, sum, na.rm = T)
  components <- likelihood(models, resp*screen, parms, theta1, theta2, Log = F)
  temp_post <- posterior(components, c(sample_mix_prop))
  #temp_post <- posterior(components, em$prior)
  out[(i+1), 1:4] <- class_probs(temp_post, data$model)[ind]
  out$deltas[i+1] <- sum(deltas*screen, na.rm = T)
  #out$n_items[i+1] <- n_items - i
  cat(i)
}

out <- out[n_short:n_items,]
gg <- data.frame(unlist(out[models]),
  rep(models, each = nrow(out)),
  rep(out$deltas, times = n_models),
  rep(out$n_items, times = n_models)
  )

head(gg)
names(gg) <- c("prob", "model", "delta", "items")

ggplot(gg, aes(x = delta, y = prob, group = model)) + geom_line(se = F, aes(color = model))
ggplot(gg, aes(x = items, y = prob, group = model)) + geom_line(aes(color = model))

# ------------------------------------------------------------
# Plausible Values
# ------------------------------------------------------------
n_reps <- 20

pv_data <- pv_gen(n_reps, resp, parms, theta1, theta2, theta1_se, theta2_se, true_model = data$model)
pv_data[models] = 0

out <- data.frame(matrix(0, nrow = n_reps, ncol <- n_models*2))
names(out) <- paste0(rep(c("prior", "se"), each = n_models), 1:n_models)

for(i in 1:n_reps) {
  ind <- pv_data$samples == i
  temp <- pv_data[ind, grep("item", names(pv_data))]
  temp_em <- EM(models, temp, parms, pv_data$theta1[ind], pv_data$theta2[ind], sorted = T)

  out[i, ] <- c(temp_em$prior, temp_em$se^2)
  pv_data[ind, models] <- temp_em$posterior
}

# ------------------------------------------------------------
# Table 1
# ------------------------------------------------------------

mean_out <- apply(out, 2, mean, na.rm = T)
var_out <- apply(out, 2, var, na.rm = T)
pv_prior <- round(mean_out[1:4], 3)
pv_se <- sqrt(mean_out[5:8] + (1 + 1/n_reps) * var_out[1:4])
pvl <- (1 + 1/n_reps) * var_out[1:4]/mean_out[5:8]

temp <- rbind(mix_prop[1,], sample_mix_prop, em$prior, em$se, pv_prior, pv_se, pvl)
row.names(temp) <- c("mix_prop", "sample_mix_prop", "em_prior", "em_se", "pv_prior", "pv_se", "pvl")
colnames(temp) <- models
xtable::xtable(temp, digits = 3)
round(temp, 4)



# ------------------------------------------------------------
# Figure 1
# ------------------------------------------------------------

pv_posterior <- lapply(pv_data[models], function(x) tapply(x, pv_data$pairs, mean)) %>% data.frame()

#raster_plot(pv_posterior, sort = F, grey_scale = F)

gg <- data.frame(unlist(pv_posterior))
names(gg) <- "prob"
gg$q <- rep(as.matrix(pv_posterior) %*% 1:n_models, times = n_models)
gg$model <- rep(Models, each = n_obs)

ggplot(gg[], aes(x = q, y = prob, group = model)) +
  stat_smooth(aes(linetype = model), lwd = 1, se = F, color = "black") +
  geom_point(aes(pch = model), color = "black", alpha = .5) +
  xlab("Expectation of posterior distribtuion") +
  ylab("Posterior probability of each model") +
  guides(pch= guide_legend(override.aes = list(alpha = 1, size = 1))) +
  ylim(c(0, 1)) +
  scale_linetype_manual(values=c(1,2,3,6)) +
  scale_shape_manual(values=c(0,1,2,3))


# ------------------------------------------------------------
# Figure 2 (??)
# ------------------------------------------------------------

temp <- as.matrix(pv_data[models]) %*% 1:4
pv_q <- tapply(temp, pv_data$pairs, mean)
pv_se <- tapply(temp, pv_data$pairs, sd)

gg <- data.frame(pv_q, pv_se, models[data$model])
names(gg) <- c("pv", "se", "model")
gg$model <- ordered(gg$model, models)

ggplot(gg, aes(x = pv, y = se, group = model)) +
  geom_point(aes(pch = model)) +
  stat_smooth(se = F, lwd = 1, color = "black", aes(group = 1)) +
  xlab("Expectation of posterior") +
  ylab("Standard error") +
  scale_shape_manual(values=c(0,1,2,3))


tapply(gg$pv, gg$model, mean)
tapply(gg$pv, gg$model, sd)

# Worth showing?
delta <- abs(data$theta1 - pv_data$theta2)
pv_delta <- tapply(delta, pv_data$pairs, mean)
plot(pv_delta, pv_se)

# ------------------------------------------------------------
# Item deltas
# ------------------------------------------------------------
pv_temp_delta <- item_delta(parms, pv_data$theta1, pv_data$theta2, sorted = T)/.25
mean_delta <- tapply(apply(pv_temp_delta, 1, mean), pv_data$pair, mean)
temp_posterior <- pv_posterior
hist(mean_delta)


temp_posterior[,2] <- temp_posterior[,2] + temp_posterior[,3]
temp_model <- data$model
temp_model[temp_model == 3] <- 2


gg <- data.frame(temp_model, mean_delta, temp_posterior[cbind(1:nrow(temp_posterior), temp_model)])

names(gg) <- c("model", "delta", "prob")
gg$model[gg$model == 1] <- "Ind"
gg$model[gg$model == 2] <- "Min + Max"
gg$model[gg$model == 4] <- "AI"
gg$model <- ordered(gg$model, c("Ind", "Min + Max", "AI"))
# Figure X

ggplot(gg, aes(x = delta, y = prob, group = model)) +
  geom_point(aes(pch = model), color = "black", alpha = .5) +
  stat_smooth(aes(linetype = model), lwd = 1, se = F, color = "black") +
  xlab("Mean Item Delta") +
  ylab("Probability of Correct Classification") +
  scale_shape_manual(values=c(0,1,3))
