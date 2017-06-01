# ------------------------------------------------------------
# Real data example for paper
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")

library("ggplot2")
library("dplyr")
source("~/github/cirt/R/cIRF_functions.R")
source("~/github/cirt/R/IRF_functions.R")
source("~/github/cirt/R/wa_functions.R")

# ------------------------------------------------------------
# Load data and estimate individual thetas
# ------------------------------------------------------------

# Load calibrated item parms
setwd("~/Dropbox/Academic/Projects/CA/Data/response_matrices")
parms <- read.csv("calibration_parms.csv", row.names = 1)
summary(parms)

# Drop DIF items
#dif_items <- "045|065"
temp_parms <- parms[!grepl(dif_items, row.names(parms)),]
items <- paste0(row.names(temp_parms), collapse = "|")
ind_parms <- col_parms <- temp_parms
row.names(ind_parms) <- paste0(row.names(ind_parms), "_IND")
row.names(col_parms) <- paste0(row.names(col_parms), "_COL")
parms <- rbind(col_parms, ind_parms)

# Load collaboration data and split into forms
collab <- read.csv("collaboration_2016_0.csv", check.names = F)
head(collab)
col_form <- format_resp(collab, row.names(temp_parms), "COL")
ind_form <- format_resp(collab, row.names(temp_parms), "IND")

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
odd <- seq(1, nrow(col_form), by = 2)

resp <- cbind(ind_form, col_form)
head(resp)

# Estimate RSC
est <- est_WA(resp, parms, SE = "exp", method = "map", parallel = F)

est$lower <- est$w - 1.96*est$w_se
est$upper <- est$w + 1.96*est$w_se
est$lower[est$lower < 0] <- 0
est$upper[est$upper < 0] <- 0
est$lower[est$lower > 1] <- 1
est$upper[est$upper > 1] <- 1
head(est)
est <- est[order(est$w), ]
est$group <- 1:nrow(est)
est$Ability <- "Within 1 SD"
est$Ability[abs(est$theta1 - est$theta2) > .5] <- "Not Within 1/2 SD"
est$Ability <- ordered(est$Ability, c("Within 1/2 SD", "Not Within 1/2 SD"))
table(est$Ability)

ggplot(est, aes(x = group, y = w, group = Ability)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color = Ability), width = .1) +
    geom_point(aes(size = Ability, shape = Ability, color = Ability)) +
    scale_color_manual(values = c("grey10", "grey70")) +
    scale_shape_discrete(solid = T) +
    ylab("Estimated RSC weight") +
    xlab("Group")  +
    theme_bw(base_size = 15) +
    geom_abline(slope = 0, intercept = 0.5, col = "black", size = 1.2) +
    scale_size_manual(values = c(4, 2)) +
    theme_bw(base_size = 15) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

17/47

# ------------------------------------------------------------
# Goodness of fit
# ------------------------------------------------------------

n_reps <- 50
n_obs <- length(est$w)

sim <- data_gen(n_reps, est$w, col_parms, est$theta1, est$theta2, NA_pattern = col_form)

sim_est <- est_WA(resp, parms, SE = "exp", method = "map", parallel = T)

map_WA(sim[grep(items, names(sim))], col_parms, sim$theta1, sim$theta2, SE = "exp", starts = sim$w, parallel = T)

sim$l <- l_WA(sim[grep(items, names(sim))], sim_est$w, col_parms, sim$theta1, sim$theta2)

l <- l_WA(col_form, est$w, col_parms, est$theta1, est$theta2)
hist(-2*l)
hist(-2*sim$l)
sim_est$w[1:500]
head(sim)

quant <- tapply(sim$logL, sim$pairs, function(x) quantile(x, p = c(.5, .025))) %>% unlist %>% matrix(, nrow = 2, ncol = n_obs) %>% t()

# Visual key for misfit
fit <- rep("<.95", times = n_obs)
fit[logL < quant[,2]] <- ">.95"
fit <- ordered(fit, c(">.95", "<.95"))

# Set up and plot
gg <- data.frame(-2*sim$logL, -2*rep(logL, each = n_reps), rep(fit, each = n_reps), rep(-2*quant[,1], each = n_reps))

names(gg) <- c("l_dist", "l_obs", "fit", "median")
gg <- gg[order(gg$median), ]
gg$pair <- rep(1:n_obs, each = n_reps)

ggplot(gg, aes(x = pair, y = l_dist, group = pair)) +
  geom_boxplot(outlier.size = 0, outlier.color = "white", aes(fill = fit)) +
  geom_point(aes(x = pair, y = l_obs, pch = fit, size = fit)) +
  scale_shape_manual(values = c(4, 20)) +
  scale_fill_grey(start = 0.3, end = 0.7) +
  xlab("Groups") +
  ylab("-2 * loglikelihood") +
  scale_size_manual(values = c(4, 1)) +
  theme(legend.title=element_blank()) +
  theme_bw(base_size = 15)





# ------------------------------------------------------------
#  WA parameter estimation
# ------------------------------------------------------------
set.seed(101)
n_reps <- 250
n_obs <- length(theta1)
w <- runif(n_obs, 0, 1)

ell <- function(i, w){
  nw <- length(w)
  t1 <- rep(theta1[i], nw)
  t2 <- rep(theta2[i], nw)
  s1 <- rep(theta1_se[i], nw)
  s2 <- rep(theta2_se[i], nw)
  R <- matrix(as.numeric(resp[i,col]), nrow = nw, ncol = nrow(parms), byrow = T)
  m_WA(R, w, parms, t1, t2)
}

w <- seq(0.0001, .99999, by = .1)
i = 1

y <- ell(i, w)
plot(w, y, type = "l", main = round(c(w[which.max(y)], ml$w[i]), 3))
sum(resp[i,], na.rm = T)
sum(ind_form[odd[i],], na.rm = T)
sum(ind_form[odd[i]+1,], na.rm = T)
resp[i,]
theta1[i]
theta2[i]
w[which.max(y)]
i = i + 1
col <- grep("COL", names(resp))

map <- map_WA(resp[odd,col], parms, theta1, theta2, SE = "exp")
map
ml <- mle_WA(resp, parms, theta1, theta2, SE = "exp")
plot(map$w, ml$w)
abline(a = 0, b = 1)
plot(map$w, map$psd)
plot(map$psd, ml$se)
abline(a = 0, b = 1)


m <- function(par, resp) {
  -1*m_full(resp, par[1], parms, par[2], par[3])
}
m(c(.999, -5, -5), resp[c(odd[i], odd[i]+ 1),])

i = 7
  q <- optim(c(.5,0,0),
         m,
         resp = resp[c(odd[i], odd[i]+ 1),],
         #method = "Nelder-Mead",
         method = "L-BFGS-B",
         lower = c(.00001, -4, -4),
         upper = c(.99999, 4, 4),
         hessian = T
         )

q$par
sqrt(diag(solve(q$hessian)))
theta1[i]
theta2[i]
theta1_se[i]
theta2_se[i]
map[i,]


points(ml$w, ml$se, col = 2)


hist(map$w)

sqrt(1/diag(map$hessian))/10
#MLE(ind_form[odd[i + 1],], parms, WMLE = T)
# Estimate weights for empirical data
ml <- mle_WA(resp, parms, theta1, theta2, SE = "exp")
plot(ml$w, ml$se)

ml <- mle_WA2(resp, parms, theta1, theta2, theta1_se, theta2_se, SE = "obs")
sum(ml$w == 1)
plot(ml$w, ml$se)

plot(abs(theta1-theta2), ml$w)
q <- apply(item_delta(parms, theta1, theta2, NA_pattern = resp), 1, mean, na.rm = T)
plot(q, ml$w)
hist(ml$w)


# Generate plausible values for theta estimates
pv <- pv_gen(n_reps, resp, parms, theta1, theta2, theta1_se, theta2_se, weights = NA)

# Estimate weights for generated data
ml_pv <- mle_WA(pv[grep(items, names(pv))], parms, pv$theta1, pv$theta2, SE = "exp", starts = pv$w, parallel = T)
pv <- cbind(pv, ml_pv)

pv_w <- tapply(pv$w, pv$pairs, mean)
var_w <- tapply(pv$se^2, pv$pairs, mean)
var_b <- tapply(pv$w, pv$pairs, var)
pv_se <- sqrt(var_w + (1 + 1/n_reps) * var_b)

var_increase <- (1 + 1/n_reps) * var_b/var_w
var_prop <- (1 + 1/n_reps) * var_b/pv_se^2
summary(var_prop)
summary(var_increase)


# SEs Figure ---------------------

gg_v <- data.frame(c(ml$se, pv_se), rep(pv_w, times = n_obs*2), rep(c("ML", "PV"), each = n_obs))
names(gg_v) <- c("SE", "PV_W", "Type")

ggplot(gg_v, aes(x = PV_W, y = SE, group = Type)) + geom_point(size = 3, aes(shape = Type, color = Type)) + ylim(c(0,.8)) +
    scale_color_manual(values = c("grey50", "black")) +
    scale_shape_discrete(solid=F) +
    theme_set(theme_grey(base_size = 15))

# ------------------------------------------------------------
# Goodness of fit
# ------------------------------------------------------------

logL <- ml$logL
sim <- data_gen(n_reps, pv_w, parms, theta1, theta2, NA_pattern = resp)
temp <- mle_WA(sim[grep(items, names(sim))], parms, sim$theta1, sim$theta2, SE = "exp", starts = sim$w, parallel = T)
sim <- cbind(sim, temp)
head(sim)

quant <- tapply(sim$logL, sim$pairs, function(x) quantile(x, p = c(.5, .025))) %>% unlist %>% matrix(, nrow = 2, ncol = length(theta1)) %>% t()

# Visual key for misfit
fit <- rep("<.95", times = length(theta1))
fit[logL < quant[,2]] <- ">.95"
fit <- ordered(fit, c(">.95", "<.95"))

# Set up and plot
gg <- data.frame(-2*sim$logL, -2*rep(logL, each = n_reps), rep(fit, each = n_reps), rep(-2*quant[,1], each = n_reps))

names(gg) <- c("l_dist", "l_obs", "fit", "median")
gg <- gg[order(gg$median), ]
gg$pair <- rep(1:length(theta1), each = n_reps)

ggplot(gg, aes(x = pair, y = l_dist, group = pair)) +
  geom_boxplot(outlier.size = 0, outlier.color = "grey90", aes(fill = fit)) +
  geom_point(aes(x = pair, y = l_obs, pch = fit, size = fit)) +
  scale_shape_manual(values = c(4, 20)) +
  scale_fill_grey(start = 0.1, end = 0.8) +
  xlab("Groups") +
  ylab("-2 * loglikelihood") +
  scale_size_manual(values = c(4, 1)) +
  theme(legend.title=element_blank()) +
  theme_set(theme_grey(base_size = 15))


# ------------------------------------------------------------
# Process loss
# ------------------------------------------------------------
dim(pv)
dim(resp)
dim(parms)
head(pv)
# Expected log likelihoods
p_add <- format_NA(cIRF("Add", parms, pv$theta1, pv$theta2), NA_pattern = resp)
e_max <- likelihood("Max", p_add, parms, pv$theta1, pv$theta2)
e_ai <- likelihood("Add", p_add, parms, pv$theta1, pv$theta2)
e_ind <- likelihood("Ind", p_add, parms, pv$theta1, pv$theta2)
e_wa <- l_WA(p_add, rep(pv_w, each = n_reps), parms, pv$theta1, pv$theta2)
e_wa2 <- l_WA(p_add, pv$w, parms, pv$theta1, pv$theta2)


# Process loss
pl_wa <- 1 - (e_ai - e_wa)/(e_ai - e_ind)
pl_wa2 <- 1 - (e_ai - e_wa2)/(e_ai - e_ind)
plot(pl_wa, pl_wa2)
pl_max <- (e_ai - e_max)

/(e_ai - e_ind)

# quantiles for wa model
quant_wa <- tapply(pl_wa2, pv$pairs, function(x) quantile(x, p = c(.5, .025))) %>% unlist %>% matrix(, nrow = 2, ncol = length(theta1)) %>% t()

# quantiles for max model
quant_max <- tapply(pl_max, pv$pairs, function(x) quantile(x, p = c(.5, .975))) %>% unlist %>% matrix(, nrow = 2, ncol = length(theta1)) %>% t()

# order groups my wa medians
ind <- order(rep(quant_wa[,1], each = n_reps))

sum(quant_wa[,2] > quant_max[2])/151

# overlapping 95% confidence intervals for wa and max?
overlap <- rep(as.factor(quant_wa[,2] > quant_max[2]), each = n_reps)
levels(overlap) <- c("< Max", "> Max")

gg <- data.frame(pl_max[ind], pl_wa[ind], overlap[ind], rep(1:n_obs, each = n_reps))
names(gg) <- c("max", "obs", "overlap", "pair")
head(gg)

ggplot(gg, aes(x = pair, y = obs, fill = overlap)) +
  geom_boxplot(aes(group = pair), outlier.size = 0, outlier.color = "grey90", size = .2) +
  scale_fill_manual(values = c("grey50", 5)) +
  ylab("1 - process loss") +
  xlab("Pair") +
  theme(legend.title=element_blank()) +
  theme_set(theme_grey(base_size = 15))


ggplot(gg, aes(x = pair, y = obs, fill = overlap)) +
  geom_boxplot(aes(group = pair), outlier.size = 0, outlier.color = "grey90", size = .2) +
  scale_fill_manual(values = c("grey50", 5)) +
  ylab("1 - process loss") +
  xlab("Pair") +
  theme(legend.title=element_blank()) +
  geom_rect(aes(xmin = 70, xmax = 151, ymin = .8, ymax = 1.05),
               fill = "transparent", color = 5, size = 1.5)


ggplot(gg[gg$pair > 70, ], aes(x = pair, y = obs, fill = overlap)) +
  geom_boxplot(aes(group = pair), outlier.size = 0, outlier.color = "grey90", size = .3) +
  scale_fill_manual(values = c("grey50", 5)) +
  ylab("1 - process loss") +
  xlab("Pair") +
  theme(legend.title=element_blank()) +
  theme(panel.border = element_rect(colour = 5, fill = NA, size = 1.5))










# ------------------------------------------------------------
# Sample Demographics
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
# DIF with Mplus: Items 45 and 65 identified as problematic in
# scalar model. After dropping, scalar model fits OK

# Metric            36.971        37       0.4704
# Scalar           101.071        74       0.0200
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
gg$DIF <- resid(lm(collab_beta ~ calib_beta, gg)) > 1
dif_items <- paste0(gg$item_names[gg$DIF == T], collapse = "|")

# DIF Figure
ggplot(gg, aes(x = calib_beta, y = collab_beta)) +
  geom_point(solid = F, aes(pch = DIF, color = DIF, size = DIF)) +
  stat_smooth(method = "lm", col = "black", se = F) +
  xlab("Calibration sample") +
  ylab("Collaborative testing condition") +
  theme(legend.position = "none") +
  scale_shape_manual(values = c(1, 4)) +
  scale_color_manual(values = c("black", "red")) +
  scale_size_manual(values = c(3, 4)) +
  theme_bw(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
