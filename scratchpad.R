# Analayses for "Pscychometric Models for Small Group Collaborations"

#devtools::install_github("peterhalpin/BearShare")
# library("BearShare")
library("ltm")
library("ggplot2")
source("~/github/cirt/R/functions.R")
NYU <- rgb(87, 6, 140, maxColorValue = 255)

# Stuff from Collab_outcomes.Rmd---------------------------------------------------------
friends
calib_ltm <- ltm(calibration ~ z1)
beta <- coef(calib_ltm)[,1]
alpha <- coef(calib_ltm)[,2]

ind_form <- col_form <- friends[,-1]
ind_form[, grep("C", names(ind_form))] <- NA
col_form[, grep("I", names(col_form))] <- NA

odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]

ind_theta <- factor.scores(calib_ltm, ind_form, type = "EB", prior = F)$score.dat$z1
col_theta <- factor.scores(calib_ltm, col_, type = "EB", prior = F)$score.dat$z1
barbell_plot(ind_theta, col_theta, legend = "left")

parms <- coef(calib_ltm)
beta_C <- parms[grep("C", row.names(parms)), 1]
alpha_C <- parms[grep("C", row.names(parms)), 2]

# New data -----------------------------------------------------------------------
  # drop  032, 075, 095, 097 based on DIF analysis
  # drop  092 because of administration problem

drop_items <- c("032", "075", "095", "097", "092")

# Load calibration data and summarize
setwd(paste0("~/Dropbox/Academic/Projects/CA/Data/response_matrices"))
calib <- read.csv("calibration_2016.csv", check.names = F)
calib <- calib[-grep(paste0(drop_items, collapse = "|"), names(calib))]
dim(calib)
summary(apply(!is.na(calib[,-1]), 2, sum))

#sort(apply(!is.na(calib[,-1]), 2, sum))

# Calibrate items
calib_ltm <- ltm(calib[,-1] ~ z1)
parms <- coef(calib_ltm) %>% data.frame
names(parms) <- c("beta", "alpha")

# Load collaboration data, split forms,
collab <- read.csv("collaboration_2016.csv", check.names = F)
names(collab)

# split into forms
col_form <- format_resp(collab, calib[,-1], "COL")
ind_form <- format_resp(collab, calib[,-1], "IND")

# Apply conjunctive scoring rule
odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]

# Drop unsuable response patterns (all 1 or all 0)
drop_groups <- c(collab$group_id[apply(col_form, 1, mean, na.rm = T)%in%c(1,0)],
  collab$group_id[apply(ind_form, 1, mean, na.rm = T)%in%c(1,0)])

length(unique(drop_groups))

col_form <-col_form[!collab$group_id%in%drop_groups,]
ind_form <-ind_form[!collab$group_id%in%drop_groups,]
odd <- seq(1, nrow(col_form), by = 2)

# Summarize
summary(apply(!is.na(col_form), 1, sum))
dim(col_form)
summary(apply(!is.na(col_form), 2, sum))

#sort(apply(!is.na(col_form), 1, sum))

summary(apply(!is.na(ind_form), 1, sum))
#sort(apply(!is.na(ind_form), 1, sum))

# Estimate theta for both forms
ind_theta <- factor.scores(calib_ltm, ind_form, type = "EB", prior = F)$score.dat$z1
col_theta <- factor.scores(calib_ltm, col_form, type = "EB", prior = F)$score.dat$z1

# Take a look
hist(col_theta)
hist(ind_theta)
barbell_plot(ind_theta[1:60], col_theta[1:60])

# Set up EM -------------------------------------------------------------------
models <- c("Ind", "Min", "Max", "AI")
theta1 <- ind_theta[odd]
theta2 <- ind_theta[odd+1]
plot(theta1, theta2)

# Sanity check
sanity <- EM(models, ind_form[odd,]*ind_form[odd+1,], theta1, theta2, parms)
sanity$prior
raster_plot(sanity)

# Collaborative repsonses
col <- EM(models, col_form[odd,], theta1, theta2, parms)
round(col$prior,3)
raster_plot(col)
class_accuracy(col)


# Effective proportion of items?
weights <- delta(parms$alpha, parms$beta, theta1, theta2)/.25

n_item <- data.frame(apply(weights*col_form, 1, sum, na.rm = T) / apply(!is.na(col_form), 1, sum))
names(n_item) <- "weights"

ggplot(n_item, aes(x = weights)) + geom_histogram(color = "white", fill = NYU, binwidth = .04)  + xlab("proportion") #+ labs(title = "Effective proportion of items for each dyad")

# Screening items
s <- screen(.05, parms$alpha, parms$beta, theta1, theta2)
col <- EM(models, col_form[odd,]*s, theta1, theta2, parms)
col$prior
class_accuracy(col)

class_accuracy <- function(EM){
  ind <- apply(EM$posterior, 1, which.max)
  arr_ind <- cbind(1:nrow(EM$posterior), ind)
  cp <- EM$posterior[arr_ind]
  tapply(cp, ind, mean)
}

raster_plot(col)

# Screening dyads
ind <- n_item[odd,]
ind > .05
col <- EM(models, col_form[ind>.05,], theta1[ind>.05], theta2[ind>.05], parms)
col$prior
raster_plot(col)



# Some follow up

temp <- col$posterior
temp_ind <- which(temp%*% 1:4 > 3.8 & temp%*% 1:4 < 4.1)
temp_ind <- sort(c(temp_ind*2, (temp_ind*2)-1))
barbell_plot(ind_theta[temp_ind], col_theta[temp_ind], "right")
table(apply(temp, 1, which.max))
hist(temp[,4])
hist(temp[,3])

apply(temp)
hist(temp%*% 1:4, breaks = 20)



#BARRRRF output of Noreen, in a rush ---------------------------------
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

twoPL(alpha, beta, -1) * (1-twoPL(alpha, beta, 1))


# Mplus DIF stuff -----------
library(MplusAutomation)
file <- "~/Dropbox/Academic/Projects/CA/Data/response_matrices/DIF/collaboration_DIF2a.out"


head(q)
q <- extractModelParameters(file)$stdy.standardized
thresholds <- q[q$paramHeader == "Thresholds",]
ind <- grep("IN", thresholds$param)
col <- grep("CO", thresholds$param)

plot(thresholds$est[ind], thresholds$est[col], pch = 20)
abline(a = -.32, b = 1)
dif <- abs((thresholds$est[ind]-.3) - thresholds$est[col]) > .34
points(thresholds$est[ind[dif]], thresholds$est[col[dif]], col = 2, pch = 20)

thresholds[dif,]


# chi square dif tests scalar invariance is OK without 39, 75, 96, 97

lc  <- -4976.624
lm  <- -5001.881
ls  <- -5019.747

ls2 <- -5017.613

cm <-  (104*1.328 - 129*1.323)/(104 - 129)
cs <-  (76*1.316 - 129*1.323)/(76 - 129)
cs2 <- (78*1.326 - 129*1.323)/(78 - 129)
csm <- (76*1.328 - 104*1.316)/(76 - 104)

lrm <- -2*(lm - lc)/cm
lrs <- -2*(ls - lc)/cs
lrs2 <- -2*(ls2 - lc)/cs2
lrsm <- 2*(lm - ls)/csm

pchisq(lrm, 129-104)
1 - pchisq(lrs, 129-76)
pchisq(lrs2, 129-78)
pchisq(lrsm, 104-76)











# Functions for WAI -----------------------------------------------------------------

WAI <- function(w1, w2, alpha, beta, theta1, theta2){
  w1 * twoPL(alpha, beta, theta1) + w2 * twoPL(alpha, beta, theta2) + (1 - w1 - w2) * twoPL(alpha, beta, theta1) * twoPL(alpha, beta, theta2)
}


R2 <- outer(theta1, theta2, WAI2, alpha, beta, w1, w2)
persp(theta1, theta2, R2, theta = -55, phi = 20, main = " ", col = "grey", expand = .5, ticktype = "detailed", nticks = 2)

# plots IRFS
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
head(temp)
temp$item <- as.factor(temp$item)
ggplot(temp, aes(x = theta, y = p, group = item)) + geom_line(size = 1.25, aes(col = item)) + geom_vline(size = 1.25, xintercept = 1, linetype = 2 ,) + geom_vline(size = 1.25, xintercept = -1, linetype = 3) + annotate("text", x = .25, y = 1.02, size = 5, label = "theta_2") +  annotate("text", x = -1.8, y = 1.02, size = 5, label = "theta_1")
