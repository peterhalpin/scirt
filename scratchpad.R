 # Analayses for "Pscychometric Models for Small Group Collaborations"

#devtools::install_github("peterhalpin/BearShare")
library("BearShare")
library("ltm")
library("ggplot2")
library("stats4")

# Stuff from Collab_outcomes.Rmd-----------------------------------------------------------

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
drop_items <- c("032", "075", "095", "097", "092")

machine <- "peterfrancis"
machine <- "peterhalpin"
setwd(paste0("/Users/", machine, "/Dropbox/Academic/Projects/CA/Data/response_matrices"))
calib <- read.csv("calibration_2016.csv", check.names = F)
names(calib)
dim(calib)
summary(apply(!is.na(calib[,-1]), 2, sum))

calib <- calib[-grep(paste0(drop_items, collapse = "|"), names(calib))]

# The whoe thing?
calib_ltm <- ltm(calib[,-1] ~ z1)
parms <- coef(calib_ltm) %>% data.frame
names(parms) <- c("beta", "alpha")
hist(parms$alpha)

collab <- read.csv("collaboration_2016.csv", check.names = F)
names(collab)

version = "IND"

format_resp <- function(collab, calib, version){
  temp <- collab[grep(version, names(collab))]
  names(temp) <- substr(names(temp), 1, 5)
  temp <- temp[names(temp)%in%names(calib)]
  temp[names(calib)[!names(calib)%in%names(temp)]] <- NA
  temp <- temp[names(calib)]
  temp
}

col_form <- format_resp(collab, calib[,-1], "COL")
ind_form <- format_resp(collab, calib[,-1], "IND")
summary(apply(!is.na(ind_form), 2, sum))

odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]

ind_theta <- factor.scores(calib_ltm, ind_form, type = "EB", prior = F)$score.dat$z1
col_theta <- factor.scores(calib_ltm, col_form, type = "EB", prior = F)$score.dat$z1

hist(col_theta)
hist(ind_theta)

# temp fix for shite
ind_theta[abs(ind_theta) > 4] <- NA
col_theta[abs(col_theta) > 4] <- NA

drop_groups <- collab$group_id[is.na(ind_theta) | is.na(col_theta)]
col_form <- col_form[!collab$group_id%in%drop_groups, ]
ind_form <- ind_form[!collab$group_id%in%drop_groups, ]
collab2 <- collab[!collab$group_id%in%drop_groups, ]
ind_theta <- ind_theta[!collab$group_id%in%drop_groups]
col_theta <- col_theta[!collab$group_id%in%drop_groups]
odd <- seq(1, nrow(col_form), by = 2)

barbell_plot(ind_theta[1:60], col_theta[1:60])

dim(ind_form)
# EM -------------------------------------------------------------------
models <- c("Ind", "Min", "Max", "AI")
theta1 <- ind_theta[odd]
theta2 <- ind_theta[odd+1]

cor(theta1, theta2)
# effective number of items?
weights <- delta(parms$alpha, parms$beta, theta1, theta2)/.25
apply(weights, 1, sum)

s <- screen(.05, parms$alpha, parms$beta, theta1, theta2)


sanity <- EM(models, ind_form[odd,]*ind_form[odd+1,]*s, theta1, theta2, parms)
sanity$prior

col <- EM(models, col_form[odd,], theta1, theta2, parms)
col$prior
dim(temp)
temp <- col$posterior


temp_ind <- which(temp%*% 1:4 > 3.6 & temp%*% 1:4 < 4.1)
temp_ind <- sort(c(temp_ind*2, (temp_ind*2)-1))

barbell_plot(ind_theta[temp_ind], col_theta[temp_ind], "right")
table(apply(temp, 1, which.max))
hist(temp[,4])
hist(temp[,1])

apply(temp)
hist(temp%*% 1:4, breaks = 20)

temp <- data.frame(cbind(1:nrow(temp), temp))
names(temp) <- c("pair", "Ind", "Min", "Max", "AI")
head(temp)

q <- reshape(temp,
  varying = names(temp)[-1],
  v.names = "prob",
  timevar = "model",
  times = names(temp)[-1],
  direction = "long"
  )

  scale_fill_gradient2( high=muted('NYU'))
  NYU <- rgb(87, 6, 140, maxColorValue = 255)
  NYU <- rgb(87, 6, 140, maxColorValue = 255)

scale_colour_manual(NYU)
q$model <- ordered(q$model, c("Ind", "Min", "Max", "AI"))
ggplot(q, aes(pair, model, fill = prob)) + geom_raster()+   scale_fill_gradient2(high = NYU)
 +theme_bw()


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
