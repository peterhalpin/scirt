# ---------------------------------------------------------------------
#  Analyses for "NAEd / Spencer presentation 2016"
# ---------------------------------------------------------------------

#devtools::install_github("peterhalpin/cirt")
#source("~/github/cirt/R/functions.R")
library("ltm")
library("ggplot2")
library("gridExtra")
library("dplyr")
library("cirt")

NYU <- rgb(87, 6, 140, maxColorValue = 255)

# Set up data ------------------------------------
  # drop  032, 075, 095, 097 based on DIF analysis
  # drop  092 because of administration problem

drop_items <- c("032", "075", "095", "097", "092")

# Calibration data
setwd(paste0("~/Dropbox/Academic/Projects/CA/Data/response_matrices"))
calib <- read.csv("calibration_2016.csv", check.names = F)
calib <- calib[-grep(paste0(drop_items, collapse = "|"), names(calib))]

# Summarize
summary(apply(!is.na(calib[,-1]), 2, sum))
#sort(apply(!is.na(calib[,-1]), 2, sum))

# Calibrate items
calib_ltm <- ltm(calib[,-1] ~ z1)
parms <- coef(calib_ltm) %>% data.frame
names(parms) <- c("beta", "alpha")

# Load collaboration data and split into forms
collab <- read.csv("collaboration_2016.csv", check.names = F)
col_form <- format_resp(collab, calib[,-1], "COL")
ind_form <- format_resp(collab, calib[,-1], "IND")

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

# Summarize
dim(col_form)
summary(apply(!is.na(col_form), 1, sum))
summary(apply(!is.na(ind_form), 1, sum))

# Estimate theta for both forms
ind_theta <- factor.scores(calib_ltm, ind_form, type = "EB", prior = F)$score.dat
$z1
col_theta <- factor.scores(calib_ltm, col_form, type = "EB", prior = F)$score.dat
$z1

# Take a look
par(mfrow  = c(1,2))
hist(col_theta)
hist(ind_theta)

# EM ----------------------------------------
  # two sanity checks and one analysis

models <- c("Ind", "Min", "Max", "AI")
theta1 <- ind_theta[odd]
theta2 <- ind_theta[odd+1]
plot(theta1, theta2)

# Sanity check #1
sanity <- EM(models, ind_form[odd,]*ind_form[odd+1,], theta1, theta2, parms)
sanity$prior
raster_plot(sanity)

# Sanity check #2
n_obs <- 5000
n_items <- 20
sim <- sim_mix(n_obs, n_items)
sim$prior
raster_plot(sim, sort = F)
table(apply(sim$posterior, 1, which.max))
class_accuracy(sim)

# Real deal: Collaborative responses
col <- EM(models, col_form[odd,], theta1, theta2, parms)
col$prior
group_score <- rep(col$posterior %*% 1:4, each = 2)
class_accuracy(col)
raster_plot(col)

# Random sample for barbell plot
temp <- sample(odd, 20)
ind <- sort(c(temp, temp + 1))
barbell_plot(ind_theta[ind], col_theta[ind], group_score[ind], legend = "right")

# Plots of group scores ------------------

# Histogram
temp <- data.frame(group_score)
names(temp) <- "u"
1 <- ggplot(temp, aes(x = u)) + geom_histogram(color = "white", fill =  "#132B43", binwidth = .25)  + xlab("\'group score\'") #+ labs(title = "Effective proportion of items for each dyad")

# Scatter plot
theta_min <- apply(cbind(theta1, theta2), 1, min)
theta_max <- apply(cbind(theta1, theta2), 1, max)
temp <- data.frame(
  c(theta_min, theta_max),
  group_score,
  rep(c("theta_min", "theta_max"), each = length(theta_min)))

temp[,1] <- rep(col_theta[odd], times = 2) - temp[,1]
names(temp) <- c("delta_theta", "U", "member")

p2 <- ggplot(temp, aes(x = U, y = delta_theta, group = member)) +
  geom_smooth(aes(color = member)) +
  geom_point(aes(color = member), size = 1) +
  scale_y_continuous(limits = c(-3,3)) +
  xlab("\'group score\'")+
  ylab("collaborative theta minus individual theta") +
  geom_abline(intercept = 0, slope = 0, col = "grey") +
  scale_color_manual(values = c("#132B43", "#56B1F7")) +
  theme(legend.position = c(.8, .2))

grid.arrange(p1, p2, ncol = 2)

#plot(col$posterior%*%1:4, col_theta[odd] -theta_min, pch = 20)
#points(col$posterior%*%1:4, col_theta[odd] -theta_max, col = 2, pch = 20)

# DIF analysis --------------------------------

library(MplusAutomation)
file <- "~/Dropbox/Academic/Projects/CA/Data/response_matrices/DIF/collaboration_DIF2a.out"
q <- extractModelParameters(file)$stdy.standardized
thresholds <- q[q$paramHeader == "Thresholds",]
ind <- grep("IN", thresholds$param)
col <- grep("CO", thresholds$param)

temp <- data.frame(thresholds$est[ind], thresholds$est[col])
names(temp) <- c("ti", "tc")
temp$dif <- (abs((thresholds$est[ind]-.3) - thresholds$est[col]) > .34)*1
temp$dif <- factor(temp$dif)

ggplot(temp, aes(x = ti, y = tc, group = dif)) + geom_point(size = 2, aes(color = dif)) +
  geom_abline(intercept = -.32, slope = 1, col = "grey") +
  ggtitle("Estimated Item Thresholds, Research Sample") +
  xlab("Individual condition") +
  ylab("Collaborative condition") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#132B43", "#56B1F7"))

# Chi square dif tests scalar invariance is OK without 39, 75, 96, 97
l_con  <- -4976.624
l_full <- -5048.337
l_metric  <- -5001.881
l_scalar  <- -5019.747

c_full <- (68*1.328 - 129*1.323)/(68 - 129)
c_metric <-  (104*1.328 - 129*1.323)/(104 - 129)
c_scalar <-  (76*1.316 - 129*1.323)/(76 - 129)

lr_full <- -2*(l_full - l_con)/c_full
lr_metric <- -2*(l_metric - l_con)/c_metric
lr_scalar <- -2*(l_scalar - l_con)/c_scalar

pchisq(lr_full, 129-68)
pchisq(lr_metric, 129-104)
pchisq(lr_scalar, 129-76)

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


#BARRRRF output for Noreen, in a rush ---------------------------------
#
# t1 <- collab2[temp_ind,]
# t1$member_id
# temp_ind <- temp_ind[-c(1:6)]
# t1 <- collab2[temp_ind,]
#
#
# t2 <- round(cbind(rep(1:(nrow(t1)/2), each = 2), t1$group_id, t1$member_id, ind_theta[temp_ind], col_theta[temp_ind], col_theta[temp_ind] - ind_theta[temp_ind]), 3)
# colnames(t2) <- c("pairs", "group_id", "member_id", "ind_theta", "collab_theta", "delta_theta")
#
# t2 <- data.frame(t2)
# groups <- unique(t2$group_id)
# table(temp$group_id)
#
# temp <- read.csv(paste0("/Users/", machine, "/Dropbox/Academic/Projects/CA/Data/serversql_0908_2016/chat_time_series.csv"))
#
# head(temp)
# gg <- table(temp$user_id)
# gg[gg%in%groups]
# temp2 <- temp[temp$user_id%in%t2$member_id, ]
# names(temp2)[1] <- "pairs"
#
#
# for(i in 1:nrow(t2)) { temp2$pairs[temp2$user_id == t2$member_id[i]] <- t2$pairs[i] }
# head(temp2)
# write.csv(temp2, "AI_chats.csv")
# table(temp2$pairs)
# temp2$pair
#
# table(temp2$group_id)
