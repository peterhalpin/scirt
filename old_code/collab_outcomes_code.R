# devtools::install_github("peterhalpin/cirt")
# library("cirt")

library("ggplot2")
library("dplyr")
source("~/github/cirt/R/cIRF_functions.R")
source("~/github/cirt/R/IRF_functions.R")

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
n_ind_items <- 25

# Data generating parameters
theta <- rnorm(n_obs*2)

beta <- sort(rnorm(n_items, mean = .35, sd = 1.3))
alpha <- runif(n_items, 1, 2.5)
parms <- data.frame(alpha, beta)
row.names(parms) <- paste0("item", 1:n_items)
mix_prop <- matrix(.25, nrow= n_models, ncol = n_obs) %>% t

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
ind <- order(theta)
theta1 <- theta[ind[1:n_obs]]
theta2 <- theta[ind[(n_obs+1):(n_obs*2)]]
theta1_se <- theta_se[ind[1:n_obs]]
theta2_se <- theta_se[ind[(n_obs+1):(n_obs*2)]]

# Pair neighbors
ind <- order(theta)
theta1 <- theta[ind[odd]]
theta2 <- theta[ind[odd+1]]
theta1_se <- theta_se[ind[odd]]
theta2_se <- theta_se[ind[odd+1]]

# Generate data
data <- data_gen(1, mix_prop, parms, theta1, theta2, expected = F)
sample_mix_prop <- table(data$model) / n_obs
head(data)

data <- data_gen(1, mix_prop, parms, theta1, theta2, expected = T)
resp <- data[grep("item", names(data))]

temp_components1 <- likelihood(models, resp, parms, theta1, theta2, Log = F)
temp_components2 <- likelihood(models, resp, parms, theta1, theta2, Log = T)

table(data$model, apply(temp_components1, 1, which.max))
table(data$model, apply(temp_components2, 1, which.max))

temp_post1 <- posterior((temp_components1), rep(.25, times = 4))
temp_post2 <- posterior(exp(temp_components2), rep(.25, times = 4))

round(class_probs(temp_post1, data$model), 3)
round(class_probs(temp_post2, data$model), 3)


# ------------------------------------------------------------
# 3 examples
# ------------------------------------------------------------

# Pair odd and even elements of theta / theta_se


odd <- seq(1, n_obs*2, by = 2)
theta1 <- theta[odd]
theta2 <- theta[odd+1]
data <- data_gen(200, mix_prop, parms, theta1, theta2, fixed_class = T)
head(data)
resp <- data[grep("item", names(data))]

data[paste0("component_", models)] <- likelihood(models, resp, parms, theta1, theta2, Log = F)

data[paste0("posterior_", models)]  <- posterior(data[paste0("component_", models)], rep(.25, times = 4))
head(data)


random
pv_post <- lapply(pv_data[models], function(x) tapply(x, pv_data$pairs, mean)) %>% data.frame()
cp <- round(class_probs(temp_post, data$model), 3)
em <- EM(models, resp, parms, theta1, theta2)
em$prior
em$se

# Pair high and low

ind <- order(theta)
far1 <- theta[ind[1:n_obs]]
far2 <- theta[ind[(n_obs+1):(n_obs*2)]]
far_data <- data_gen(1, mix_prop, parms, far1, far2, expected = F)
far_resp <- far_data[grep("item", names(far_data))]
temp_components <- exp(likelihood(models, far_resp, parms, far1, far2, Log = T))
temp_post <- posterior(temp_components, rep(.25, times = 4))
far_cp <- round(class_probs(temp_post, far_data$model), 3)
far_em <- EM(models, far_resp, parms, far1, far2)
far_em$prior
far_em$se


# Pair neighbors
ind <- order(theta)
near1 <- theta[ind[odd]]
near2 <- theta[ind[odd+1]]
near_data <- data_gen(1, mix_prop, parms, near1, near2, expected = F)
near_resp <- near_data[grep("item", names(near_data))]
temp_components <- exp(likelihood(models, near_resp, parms, near1, near2, Log = T))
temp_post <- posterior(temp_components, rep(.25, times = 4))
near_cp <- round(class_probs(temp_post, near_data$model), 3)
near_em <- EM(models, near_resp, parms, near1, near2)
near_em$prior
near_em$se

# ------------------------------------------------------------
# CP versus item deltas
# ------------------------------------------------------------
t1 <- random1
t2 <- random2
data <- random_data
resp <- data[grep("item", names(data))]
classify <- random_cp


t1 <- far1
t2 <- far2
data <- far_data
resp <- data[grep("item", names(data))]
classify <- far_cp

t1 <- near1
t2 <- near2
data <- near_data
resp <- data[grep("item", names(data))]


deltas <- item_delta(parms, t1, t2)
delta_order <- apply(deltas, 1, order) %>% t

n_short <- 50
j <- n_short
cp_ind <- cbind(1:4, 1:4)
cp_ind <- cbind(1:4, c(2,3,2,3))

out <- data.frame(rbind((classify)[cp_ind]), sum(deltas))
names(out) <- c(models, "deltas")


for(j in n_short:n_items) {
  screen <- deltas*NA
  i <- j - n_short + 1
  ind <- cbind(rep(1:n_obs, times = n_short), c(delta_order[,i:j]))
  screen[ind] <- 1
  components <- likelihood(models, resp*screen, parms, t1, t2, Log = F)
  temp_post <- posterior(components, rep(.25, times = 4))
  out[i, 1:4] <- class_probs(temp_post, data$model)[cp_ind]
  out$deltas[i] <- sum(deltas*screen, na.rm = T)
  cat(i)
}
out
# q_out <- out
gg <- data.frame(unlist(out[models]),
  rep(models, each = nrow(out)),
  rep(out$deltas, times = n_models)
  )

names(gg) <- c("prob", "model", "delta")
head(gg)

#gg1 <- gg

ggplot(gg, aes(x = delta, y = prob, group = model)) +
  geom_line(aes(color = model)) +
  #geom_smooth(se = F, aes(color = model)) +
  ylim(c(0, 1))

ggplot(gg, aes(x = items, y = prob, group = model)) + geom_line(aes(color = model))


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


# ------------------------------------------------------------
# real data: DIF
# Items 45 and 65 identified as problematic in scalar model
# After dropping, scalar model fits OK

# Metric            36.971        37       0.4704
# Scalar            101.071       74       0.0200
# Scalar (w/ drop)  84.358        70       0.1161
# ------------------------------------------------------------

 # Load calibrated item parms
setwd("~/Dropbox/Academic/Projects/CA/Data/response_matrices")
parms <- read.csv("calibration_parms.csv", row.names = 1)

scalar <- "~/Dropbox/Academic/Projects/CA/Data/response_matrices/DIF/collaboration_DIF_scalar.out"

grep_items <- paste0(row.names(parms), collapse = "|")
q <- readLines(scalar)
Begin <- grep("Item Difficulties", q)
End <- grep("Variances", q)

calib_temp <- q[Begin[1]:End[End > Begin[1] & End < Begin[2]]]
collab_temp <- q[Begin[2]:End[End > Begin[2] & End < Begin[3]]]

calib_beta <- substr(calib_temp[grepl(grep_items, calib_temp)], 23, 28) %>% as.numeric
collab_beta <- substr(collab_temp[grepl(grep_items, collab_temp)], 23, 28) %>% as.numeric
item_names <-substr(collab_temp[grepl(grep_items, collab_temp)], 6, 10)

gg <- data.frame(item_names, collab_beta, calib_beta)
head(gg)
gg$dif <- resid(lm(collab_beta ~ calib_beta, gg)) > 1


ggplot(gg, aes(x = calib_beta, y = collab_beta)) +
  geom_point(aes(pch = dif, color = dif, size = dif)) +
  stat_smooth(method = "lm", col = "black", se = F) +
  # ggtitle("Differential item functioning: difficulty estimates") +
  xlab("Calibration sample") +
  ylab("Collaborative testing condition") +
  # theme_bw() +
  theme(legend.position = "none") +
  scale_shape_manual(values = c(20, 4)) +
  scale_color_manual(values = c("black", "red")) +
  scale_size_manual(values = c(2, 4))

dim(gg)
dif_items <- paste0(gg$item_names[gg$dif == T], collapse = "|")

# ------------------------------------------------------------
# real data example
# ------------------------------------------------------------

 # Load calibrated item parms
setwd("~/Dropbox/Academic/Projects/CA/Data/response_matrices")
parms <- read.csv("calibration_parms.csv", row.names = 1)
summary(parms)
# drop DIF items
#dif_items <- "045|065"
parms <- parms[!grepl(dif_items, row.names(parms)),]

#  collaboration data and split into forms
collab <- read.csv("collaboration_2016.csv", check.names = F)
col_form <- format_resp(collab, row.names(parms), "COL")
ind_form <- format_resp(collab, row.names(parms), "IND")
head(ind_form)

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
head(resp)

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
hist(-2*logL)

# Simulate null distribution
n_reps <- 250
mix_prop <- em$posterior

temp <- data_gen(n_reps, mix_prop, parms, theta1, theta2, theta1_se, theta2_se, NA_pattern = resp)

temp_components <- likelihood(models, temp[,grep(grep_items, names(temp))], parms, temp$theta1, temp$theta2, sorted = T, Log = F)


temp_logL <- incomplete_data(temp_components, temp[models], Sum = F)
q <- tapply(temp_logL, temp$pair, function(x) quantile(x, p = c(.5, .025))) %>% unlist %>% matrix(, nrow = 2, ncol = length(theta1)) %>% t()

fit <- rep("<.95", times = length(theta1))
fit[logL < q[,2]] <- ">.95"
fit <- ordered(fit, c(">.95", "<.95"))

gg <- data.frame(-2*temp_logL, -2*rep(logL, each = n_reps), rep(fit, each = n_reps), rep(-2*q[,1], each = n_reps))

names(gg) <- c("l_dist", "l_obs", "fit", "median")
gg <- gg[order(gg$median), ]
gg$pair <- rep(1:length(theta1), each = n_reps)

ggplot(gg, aes(x = pair, y = l_dist, group = pair)) +
  geom_boxplot(outlier.size = 0, outlier.color = "white", aes(fill = fit)) +
  geom_point(aes(x = pair, y = l_obs, pch = fit, size = fit)) +
  scale_shape_manual(values = c(4, 20)) +
  # theme_bw() +
  scale_fill_grey(start = 0.1, end = 0.8) +
  xlab("Groups") +
  ylab("-2 * loglikelihood") +
  scale_size_manual(values = c(4, 1)) +
  theme(legend.title=element_blank())



# ------------------------------------------------------------
# real data:  Plausible Values
# ------------------------------------------------------------
n_reps <- 200
theta1_se
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
# real data:  Table 1
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
# real data: raster plot
# ------------------------------------------------------------
# ------------------------------------------------------------
# Figure 1
# ------------------------------------------------------------

pv_posterior <- lapply(pv_data[models], function(x) tapply(x, pv_data$pairs, mean)) %>% data.frame()

#raster_plot(pv_posterior, sort = F, grey_scale = F)

gg <- data.frame(unlist(pv_posterior))
names(gg) <- "prob"
gg$q <- rep(as.matrix(pv_posterior) %*% 1:n_models, times = n_models)
gg$model <- rep(Models, each = nrow(resp))

ggplot(gg[], aes(x = q, y = prob, group = model)) +
  stat_smooth(aes(linetype = model), lwd = 1, se = T, color = "black") +
  geom_point(aes(pch = model), color = "black", alpha = .5) +
  xlab("Expectation of posterior distribtuion") +
  ylab("Posterior probability of each model") +
  guides(pch= guide_legend(override.aes = list(alpha = 1, size = 1))) +
  ylim(c(-.05, 1)) +
  scale_linetype_manual(values=c(1,2,3,6)) +
  scale_shape_manual(values=c(0,1,2,3)) +
   guides(fill=guide_legend(override.aes=list(fill=NA)))



# ------------------------------------------------------------
# real data: standard error of posterior
# ------------------------------------------------------------

temp <- as.matrix(pv_data[models]) %*% 1:4
pv_q <- tapply(temp, pv_data$pairs, mean)
pv_se <- tapply(temp, pv_data$pairs, sd)
model <- apply(em$posterior, 1, function(x) models[which.max(x)])
gg <- data.frame(pv_q, pv_se, model)
names(gg) <- c("pv", "se", "model")

ggplot(gg, aes(x = pv, y = se, group = model)) +
  geom_point(aes(pch = model)) +
  stat_smooth(se = F, lwd = 1, color = "black", aes(group = 1)) +
  xlab("Expectation of posterior") +
  ylab("Standard error") +
  ggtitle("Standard error of expectation of posterior distribution") +
  theme_bw()

tapply(gg$pv, gg$model, mean)
tapply(gg$pv, gg$model, sd)

# ------------------------------------------------------------
# real data: Item deltas
# ------------------------------------------------------------
pv_temp_delta <- item_delta(parms, pv_data$theta1, pv_data$theta2, sorted = T)/.25
mean_delta <- tapply(apply(pv_temp_delta, 1, mean), pv_data$pair, mean)
temp_posterior <- pv_posterior
hist(mean_delta)

temp_posterior[,2] <- temp_posterior[,2] + temp_posterior[,3]
temp_model <- apply(temp_posterior, 1, which.max)
temp_model[temp_model == 3] <- 2


gg <- data.frame(temp_model, mean_delta, temp_posterior[cbind(1:nrow(temp_posterior), temp_model)])

names(gg) <- c("model", "delta", "prob")
gg$model[gg$model == 1] <- "Ind"
gg$model[gg$model == 2] <- "Min + Max"
gg$model[gg$model == 4] <- "AI"
head(gg)
# Figure X

ggplot(gg, aes(x = delta, y = prob, group = model)) +
  geom_point(color = "grey60", aes(pch = model)) +
  stat_smooth(aes(linetype = model), color = "black", se = F ) +
  xlab("Mean Item Delta") +
  ylab("Probability of Correct Classification") +
  ggtitle("Posterior Probability of Correct Model Classification") +
  theme_bw() + scale_color_grey() + scale_fill_grey()


# ------------------------------------------------------------
# real data: person fit
# ------------------------------------------------------------

mean_l <- function
var_l

scalar <- "~/Dropbox/Academic/Projects/CA/Data/response_matrices/DIF/collaboration_DIF_scalar.out"

grep_items <- paste0(row.names(parms), collapse = "|")
q <- readLines(scalar)
Begin <- grep("Item Difficulties", q)
End <- grep("Variances", q)

calib_temp <- q[Begin[1]:End[End > Begin[1] & End < Begin[2]]]
collab_temp <- q[Begin[2]:End[End > Begin[2] & End < Begin[3]]]

calib_beta <- substr(calib_temp[grepl(grep_items, calib_temp)], 23, 28) %>% as.numeric
collab_beta <- substr(collab_temp[grepl(grep_items, collab_temp)], 23, 28) %>% as.numeric
item_names <-substr(collab_temp[grepl(grep_items, collab_temp)], 6, 10)

gg <- data.frame(item_names, collab_beta, calib_beta)
gg$dif <- resid(lm(collab_beta ~ calib_beta, gg)) > 1


ggplot(gg, aes(x = calib_beta, y = collab_beta)) +
  geom_point(size = 2, aes(pch = dif, color = dif)) +
  stat_smooth(method = "lm", col = "black", se = F) +
  ggtitle("Differential item functioning: difficulty estimates") +
  xlab("Calibration sample") +
  ylab("Collaborative testing condition") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_shape_manual(values = c(20, 4)) +
  scale_color_manual(values = c("black", "red"))


dif_items <- gg$item_names[gg$dif == T]


  # ------------------------------------------------------------
  # real data: DIF
  # Items 45 and 65 identified as problematic in scalar model
  # After dropping, scalar model fits
  #Models Compared
  # Metric            36.971        37       0.4704
  # Scalar            101.071       74       0.0200
  # Scalar (w/ drop)  84.358        70       0.1161

  # ------------------------------------------------------------




# Plot IRFS for numerical example ---------------------
theta <- seq(-4, 4, by = .25)
nIRF <- 4
a <- 2
alpha <- c(a, a, a, a/3)
beta <- c(-1, 1, 0, 0)
temp <- cbind(rep(1:nIRF, each = length(theta)), rep(theta, times = length(nIRF)))
q <- c()
for(i in 1:nIRF) q <- c(q, twoPL(alpha[i], beta[i], theta))

temp <- data.frame(cbind(temp, q))
names(temp) <- c("item", "theta", "p")
temp$item <- as.factor(temp$item)

ggplot(temp, aes(x = theta, y = p, group = item)) +
  geom_line(size = 1.25, aes(col = item)) +
  geom_vline(size = 1.25, xintercept = 1, linetype = 2 ,) +
  geom_vline(size = 1.25, xintercept = -1, linetype = 3) +
  annotate("text", x = .25, y = 1.02, size = 5, label = "theta_2") +
  annotate("text", x = -1.8, y = 1.02, size = 5, label = "theta_1"))

# Weights items ------------------

weights <- delta(parms$alpha, parms$beta, theta1, theta2)/.25
temp <- data.frame(apply(weights*col_form, 1, sum, na.rm = T) / apply(!is.na(col_form), 1, sum))
names(temp) <- "weights"

ggplot(n_item, aes(x = weights)) +
  geom_histogram(color = "white", fill =  "#132B43", binwidth = .04)  +
  xlab("proportion") +
  labs(title = "Effective proportion of items for each dyad")

# Screen items
s <- screen(.06, parms$alpha, parms$beta, theta1, theta2)
col2 <- EM(models, col_form[odd,]*s, theta1, theta2, parms)
col2$prior
class_accuracy(col2)
raster_plot(col2)

# Table for slides
xtable::xtable(round(rbind(col$prior,
  class_accuracy(col),
  col2$prior,
  class_accuracy(col2)),3))

# Screening dyads ?
ind <- temp[odd,]
temp_form <- col_form[odd,]
col3 <- EM(models, temp_form[ind>.05,], theta1[ind>.05], theta2[ind>.05], parms)
col3$prior
raster_plot(col3)
class_accuracy(col3)
plot(col3$posterior%*%1:4, col$posterior[ind > .05,]%*%1:4)
