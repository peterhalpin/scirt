# devtools::install_github("peterhalpin/cirt")
# library("cirt")
rm(list = ls())
library("ggplot2")
library("dplyr")
source("~/github/cirt/R/cIRF_functions.R")
source("~/github/cirt/R/IRF_functions.R")
source("~/github/cirt/R/bootstrapping.R")

# ------------------------------------------------------------
# Load data and estimate individual thetas
# ------------------------------------------------------------

# Load calibrated item parms
setwd("~/Dropbox/Academic/Projects/CA/Data/response_matrices")
parms <- read.csv("calibration_parms.csv", row.names = 1)
summary(parms)

# Drop DIF items
dif_items <- "045|065"
parms <- parms[!grepl(dif_items, row.names(parms)),]
items <- paste(row.names(parms), collapse = "|")

# Load collaboration data and split into forms
collab <- read.csv("collaboration_2016_0.csv", check.names = F)
col_form <- format_resp(collab, row.names(parms), "COL")
ind_form <- format_resp(collab, row.names(parms), "IND")

 # Apply conjunctive scoring rule
odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]

# Drop 12 unsuable response patterns (all 1 or all 0)?
drop_group <- c(
  collab$group_id[apply(col_form, 1, mean, na.rm = T) %in% c(0,1)],
  collab$group_id[apply(ind_form, 1, mean, na.rm = T) %in% c(0,1)]) %>% unique

col_form <-col_form[!collab$group_id%in%drop_group,]
ind_form <-ind_form[!collab$group_id%in%drop_group,]


# Estimate theta for ind forms
ind <- MLE(ind_form, parms, WMLE = T)
plot(ind$theta, ind$se)

# Reset odd for dropped items
odd <- seq(1, nrow(col_form), by = 2)

theta1 <- ind$theta[odd]
theta2 <- ind$theta[odd+1]
theta1_se <- ind$se[odd]
theta2_se <- ind$se[odd+1]

resp <- col_form[odd, ]



# ------------------------------------------------------------
#  EM + PV
# ------------------------------------------------------------

# Constants
n_reps <- 250
n_obs <- length(theta1)
rep_order <- order(rep(1:length(theta1), times = n_reps))
models <- c("Ind", "Min", "Max", "AI")
Models <- ordered(models, models)
n_models <- 3

# Run EM on observed data
em <- EM(models, resp, parms, theta1, theta2)
em$prior
# Generate PV
set.seed(101)
pv <- pv_gen(n_reps, resp, parms, theta1, theta2, theta1_se, theta2_se)

# Run EM on PV data
fun <- function(i){
  ind <- pv$samples == i
  temp <- pv[ind, grep(items, names(pv))]
  EM(models, temp, parms, pv$theta1[ind], pv$theta2[ind], sorted = T)
}

temp_em <- parallel::mclapply(1:n_reps, fun)

# Save priors and se
out <- lapply(temp_em, function(x) c(x$prior, x$se^2)) %>% unlist %>% matrix(nrow = n_reps, ncol = n_models*2, byrow = T) %>% data.frame
names(out) <- paste0(rep(c("prior", "se"), each = n_models), 1:n_models)
out <- out[!apply(is.na(out), 1, sum) > 0,]

# Save posteriors
temp <- lapply(temp_em, function(x) x$posterior) %>% {do.call(rbind, .)}
pv[models] <- temp[rep_order,]

# Compute PV on priors
pv_prior <- apply(out[1:n_models], 2, mean)
var_w <- apply(out[(n_models+1):(n_models*2)], 2, mean)
var_b <- apply(out[1:n_models], 2, var)
pv_se <- sqrt(var_w + (1 + 1/n_reps) * var_b)
var_increase <- (1 + 1/n_reps) * var_b/var_w
var_prop <- (1 + 1/n_reps) * var_b/pv_se^2
summary(var_prop)
summary(var_increase)

# Table of parameters
temp <- rbind(em$prior, em$se, pv_prior, pv_se, var_increase)
row.names(temp) <- c("em_prior", "em_se", "pv_prior", "pv_se", "var_increase")
colnames(temp) <- models
xtable::xtable(temp, digits = 3)
round(temp, 4)



# ------------------------------------------------------------
# Person fit using contribution to likelihood
# Using prior or posterior? Posterior is more consistent with person fit...
# ------------------------------------------------------------

# Get observed values
components <- likelihood(models, resp, parms, theta1, theta2, Log = F)

# Using PV posterior
#logl <- incomplete_data(components, em$posterior, Sum = F)
#mix_prop <- em$posterior

# Using PV prior
logl <- incomplete_data(components, pv_prior, Sum = F)
mix_prop <- matrix(pv_prior, nrow = n_obs, ncol = n_models, byrow = T)

# Simulate null distribution
set.seed(101)
temp <- data_gen(n_reps, mix_prop, parms, theta1, theta2, theta1_se, theta2_se, NA_pattern = resp)

temp_components <- likelihood(models, temp[,grep(items, names(temp))], parms, temp$theta1, temp$theta2, sorted = T, Log = F)

temp_logl <- incomplete_data(temp_components, temp[models], Sum = F)

# Get quantiles
quant <- tapply(temp_logl, temp$pair, function(x) quantile(x, p = c(.5, .025))) %>% unlist %>% matrix(, nrow = 2, ncol = length(theta1)) %>% t()

# Visual key for misfit
fit <- rep(1, times = length(theta1))
fit[logl < quant[,2]] <- 0
fit <- as.factor(fit)

# Plot
gg <- data.frame(-2*temp_logl, -2*rep(logl, each = n_reps), rep(fit, each = n_reps), rep(-2*quant[,1], each = n_reps))

names(gg) <- c("l_dist", "l_obs", "fit", "median")
gg <- gg[order(gg$median), ]
gg$pair <- rep(1:length(theta1), each = n_reps)
e1 <- expression(phantom(0) %in% ".95 CI")
e2 <- expression(phantom(0) %notin% ".95 CI")

ggplot(gg, aes(x = pair, y = l_dist, group = pair)) +
  geom_boxplot(outlier.size = 0, outlier.color = "grey90", size = .1, aes(fill = fit)) +
  geom_point(aes(x = pair, y = l_obs, pch = fit, size = fit)) +
  scale_shape_manual(values = c(4, 20), guide = F) +
  xlab("Groups") +
  ylab("-2 * log-likelihood") +
  scale_size_manual(values = c(4, 1), guide = F) +
  ylim(c(0, 45)) +
  theme_set(theme_grey(base_size = 15)) +
  scale_fill_grey(start = 0.1, end = 0.6, breaks=c(1, 0), name = "", labels=c(e1, e2))


# ------------------------------------------------------------
# Process loss
# ------------------------------------------------------------

# Expected log likelihoods
p_ai <- format_NA(cIRF("AI", parms, pv$theta1, pv$theta2), NA_pattern = resp)
#p_ai <-  cIRF("AI", parms, pv$theta1, pv$theta2)
e_max <- likelihood("Max", p_ai, parms, pv$theta1, pv$theta2, sorted = T)
e_ai <- likelihood("AI", p_ai, parms, pv$theta1, pv$theta2)
e_ind <- likelihood("Ind", p_ai, parms, pv$theta1, pv$theta2)

out_prior <- lapply(pv_prior, function(x) rep(x, times = n_obs*n_reps)) %>% data.frame
components <- likelihood(models, p_ai, parms, pv$theta1, pv$theta2, sorted = T)

e_mod <-  apply(components * pv[models], 1, sum)
# + apply(log(out_prior) * pv[models], 1, sum, na.rm = T) -
#    apply(log(pv[models])* pv[models], 1, sum, na.rm = T)

pl_mod <- 1 - (e_ai - e_mod) / (e_ai - e_ind)
pl_max <- 1 - (e_ai - e_max) / (e_ai - e_ind)
min(pl_mod)

# quantiles for mix model
quant_mod <- tapply(pl_mod, pv$pairs, function(x) quantile(x, p = c(.5, .025))) %>% unlist %>% matrix(, nrow = 2, ncol = length(theta1)) %>% t()

# quantiles for max model
quant_max <- tapply(pl_max, pv$pairs, function(x) quantile(x, p = c(.5, .975))) %>% unlist %>% matrix(, nrow = 2, ncol = length(theta1)) %>% t()

# order groups my wa medians
ind <- order(rep(quant_mod[,1], each = n_reps))
sum(quant_mod[,2] > quant_max[,1])/n_obs


# overlapping 95% confidence intervals for wa and max?
overlap <- rep(as.factor(quant_mod[,2] > quant_max[,1]), each = n_reps)
levels(overlap) <- c("< Max", "> Max")


gg <- data.frame(pl_max[ind], pl_mod[ind], overlap[ind], rep(1:n_obs, each = n_reps))
names(gg) <- c("max", "obs", "overlap", "pair")
head(gg)

ggplot(gg, aes(x = pair, y = obs, fill = overlap)) +
  geom_boxplot(aes(group = pair), outlier.size = 0, outlier.color = "grey90", size = .2) +
  scale_fill_manual(values = c("grey50", 5)) +
  ylab("1 - process loss") +
  xlab("Pair") +
  theme(legend.title=element_blank()) +
  theme_set(theme_grey(base_size = 15)) +
  geom_rect(aes(xmin = 125, xmax = 162, ymin = .9, ymax = 1.02),
                 fill = "transparent", color = "transparent", size = 1.5)


ggplot(gg, aes(x = pair, y = obs, fill = overlap)) +
  geom_boxplot(aes(group = pair), outlier.size = 0, outlier.color = "grey90", size = .2) +
  scale_fill_manual(values = c("grey50", 5)) +
  ylab("1 - process loss") +
  xlab("Pair") +
  theme(legend.title=element_blank()) +
  geom_rect(aes(xmin = 80, xmax = 162, ymin = .8, ymax = 1.02),
               fill = "transparent", color = 5, size = 1.5)


ggplot(gg[gg$pair > 76, ], aes(x = pair, y = obs, fill = overlap)) +
  geom_boxplot(aes(group = pair), outlier.size = 0, outlier.color = "grey90", size = .3) +
  scale_fill_manual(values = c("grey50", 5)) +
  ylab("1 - process loss") +
  xlab("Pair") +
  ylim(c(.75, 1)) +
  theme(legend.title=element_blank()) +
  theme(panel.border = element_rect(colour = 5, fill = NA, size = 1.5))
















# ------------------------------------------------------------
# Process loss
# ------------------------------------------------------------

# Expected log likelihoods
p_max <- matrix(c(0,0,1,0), ncol = 4, nrow = nrow(pv_data), byrow = T)
e_max <- e_likelihood(p_max, parms, pv_data$theta1, pv_data$theta2, sorted = T)
e_mix <- e_likelihood(pv_data[models], parms, pv_data$theta1, pv_data$theta2, sorted = T)

# Process loss
pl_max <- 1 - PL(e_max, parms, pv_data$theta1, pv_data$theta2)
pl_mix <- 1 - PL(e_mix, parms, pv_data$theta1, pv_data$theta2)

quant <- tapply(pl_mix, pv_data$pair, function(x) quantile(x, p = c(.5, .05))) %>% unlist %>% matrix(, nrow = 2, ncol = length(theta1)) %>% t()

ind <- order(rep(quant[,1], each = n_reps))
ref <- rep(as.factor(quant[,2] > tapply(pl_max, pv_data$pair, median)), each = n_reps)
levels(ref) <- c("", "> Max")
ref2 <- rep(tapply(pl_max, pv_data$pair, mean), each = n_reps)
gg <- data.frame(pl_max[ind], pl_mix[ind], ref[ind], ref2[ind], rep(1:n_obs, each = n_reps))
names(gg) <- c("max", "obs", "ref", "ref2", "pair")
head(gg)

ggplot(gg, aes(x = pair, y = obs, fill = ref)) +
  geom_boxplot(aes(group = pair), outlier.size = 0, outlier.color = "grey90", size = .05) +
  scale_fill_manual(values = c("grey70", "grey20")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("1 - process loss") +
  xlab("Pair") +
  theme(legend.title=element_blank())


# And thats all we need???


# ------------------------------------------------------------
# PV posterior
# ------------------------------------------------------------

# Get the "true model" from PVs
temp <- lapply(pv_data[models], function(x) tapply(x, pv_data$pairs, mean)) %>% data.frame
pv_data$model <- rep(apply(temp, 1, which.max), each = n_reps)
pv_data$model[pv_data$model == 3] <- 2
pv_data$model[pv_data$model == 4] <- 3

cp  <- boot_cp(pv_data, pv_prior, parms)
post <- round(cp$mean, 3)
post_se <- round(cp$se, 3)
xtable::xtable(paste0(post, " (", post_se, ") ") %>% array(, dim = c(3,3)))

boot_cp <- function(data, mix_prop, parms) {
  n_reps <- max(data$samples)
  models <- c("Ind", "Min", "Max", "AI")
  item_names <- paste(row.names(parms) , collapse = "|")
  resp <- data[grep(item_names, names(data))]
  components <- likelihood(models, resp, parms, data$theta1, data$theta2, sorted = T, Log = F)
  post  <- posterior(components, mix_prop)
  post <- cbind(post[,1], post[,2] + post[,3], post[,4])
  temp_cp <- lapply(1:n_reps, function(x) class_probs(post[data$samples == x,], data$model[data$samples == x]))
  mean <- Reduce(`+`, temp_cp) / n_reps
  se <- Reduce(`+`, lapply(temp_cp, function(x) (x - mean)^2 / (n_reps-1))) %>% sqrt
  out <- list(mean, se)
  names(out) <- c("mean", "se")
  out
}

# ------------------------------------------------------------
# Posteriors by deltas
# ------------------------------------------------------------

# Get the "true model" from PVs

gg <- lapply(pv_data[models], function(x) tapply(x, pv_data$pairs, mean)) %>% unlist %>% data.frame

head(gg)
names(gg) <- "prob"
gg$model <- rep(Models, each = length(theta1))

temp <- item_delta(parms, pv_data$theta1, pv_data$theta2, sorted = F, NA_pattern = resp)/.25
temp[temp > 1] <- 1
temp <- apply(temp, 1, mean, na.rm = T)
hist(temp)
gg$delta <-  tapply(temp, pv_data$pairs, mean) %>% rep(4)
head(gg, 200)
gg$delta2 <- as.ordered(round(gg$delta,10))

ggplot(gg, aes(x = delta2, y = prob, fill = model)) +
    geom_bar(position = "fill", stat = "identity")  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_x_discrete(breaks = levels(gg$delta2)[seq(1, 151, 5)], label = abbreviate) +
    #scale_fill_grey(start = 0.8, end = 0.2, na.value = "grey50")
     scale_fill_manual(values = c("grey80", "grey50", "grey50", "grey20"))

t1 <- t2 <- seq(-3, 3, .1)

d <- function(t) {
  -apply(item_delta(parms, t[1], t[2])/.25, 1, mean, na.rm = T)
}

optim(c(0,0), d)

d <- function(t1, t2) {
  apply(item_delta(parms, t1, t2)/.25, 1, mean, na.rm = T)
}

q <- outer(t1,t2, d)

persp(t1, t2, q, theta = 15, phi = 25, main = "SE", expand = .5, ticktype = "detailed", nticks = 5)


# ------------------------------------------------------------
# real data: classification probabilities
# ------------------------------------------------------------

gg <- data.frame(unlist(pv_data[models]))
names(gg) <- "prob"
gg$q <- rep(as.matrix(pv_data[models]) %*% 1:n_models, times = n_models)
gg$model <- rep(Models, each = nrow(pv_data))
gg$sample <- sample.int(500, nrow(gg), replace = T)
sum(gg$sample == 1)
head(gg)
p1 <-
ggplot(gg, aes(x = q, y = prob, group = model)) +
  stat_smooth(aes(linetype = model), se = T, lwd = 1, color = "black") +
  geom_point(data = gg[gg$sample > 1,], aes(pch = model), color = "black", alpha = .5) +
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
gg$model <-  Models[gg$model]
gg <- data.frame(pv_q, pv_se, model)
names(gg) <- c("pv", "se", "model")
gg$model

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
# Demographics
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
  scale_size_manual(values = c(3, 4))





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
