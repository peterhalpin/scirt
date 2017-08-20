# ------------------------------------------------------------
# Real data example for paper
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")

library("ggplot2")
library("dplyr")
install.packages("dplyr")
source("~/github/cirt/R/cIRF_functions.R")
source("~/github/cirt/R/IRF_functions.R")

# ------------------------------------------------------------
# Load data and estimate individual thetas
# ------------------------------------------------------------

# Load calibrated item parms
setwd("~/Dropbox/Academic/Projects/CA/Data/response_matrices")
parms <- read.csv("calibration_parms.csv", row.names = 1)
summary(parms)

# Drop DIF items
dif_items <- "045|065"
temp_parms <- parms[!grepl(dif_items, row.names(parms)),]
items <- paste0(row.names(temp_parms), collapse = "|")
ind_parms <- col_parms <- temp_parms
row.names(ind_parms) <- paste0(row.names(ind_parms), "_IND")
row.names(col_parms) <- paste0(row.names(col_parms), "_COL")
parms <- rbind(col_parms, ind_parms)

# Load collaboration data and split into forms
collab <- read.csv("collaboration_2016_0.csv", check.names = F)
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
est
head(est)
legend = "none"
plot(c(est$w, est$w), c(est$theta1, est$theta2))
cor(c(est$w, est$w), c(est$theta1, est$theta2))
barbell_plot( c(est$theta1, est$theta2), c(est$w, est$w))

ind_theta <- c(rbind(est$theta1, est$theta2))
col_theta <- rep(est$w, each = 2)
data <- data.frame(ind_theta, col_theta)
data$pairs <- factor(rep(1:(length(ind_theta)/2), each = 2))
if (is.null(group_score)) {
  data$group_score <- data$pairs
  legend_title <- "pairs"
} else {
  data$group_score <- group_score
  legend_title <- "group_score"
}
ggplot(data = data, aes(x = ind_theta, y = col_theta, group = pairs)) +
  geom_line(aes(color = group_score)) +
  geom_point(aes(color = group_score), size = 4) +
  #scale_x_continuous(limits = lim) +
  #scale_y_continuous(limits = lim) +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  theme(legend.position = legend) +
  labs(color = legend_title) +
  ggtitle("Collaborative vs Individual Performance") +
  xlab("Individual Theta")+
  ylab("Collaborative Theta")+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)
 )
data$diff <- rep(est$theta1 - est$theta2, each = 2)
data$se <- abs(data$diff) < 1
(data$diff[data$se == T]
  sum(data$se)/2
cor((est$theta1 + est$theta2) / 2, est$w)
)
col
 ggplot(data = data[data$se == T, ], aes(x = ind_theta, y = col_theta, group = pairs)) +
   geom_point(aes(color = group_score), size = 4) +
   geom_line(aes(color = group_score)) +
   theme(legend.position = "none") +
   #ggtitle("Collaborative vs Individual Performance") +
   xlab("Individual Proficiency")+
   ylab("Baseline Group Performance")+
   theme(axis.text.x = element_text(size = 13),
         axis.text.y = element_text(size = 13)
  )

plot(est$w, est$theta2)

# Plot confidence intervals
gg <- est
gg$lower <- gg$w - 1.96*gg$w_se
gg$upper <- gg$w + 1.96*gg$w_se
gg$lower[gg$lower < 0] <- 0
gg$upper[gg$upper < 0] <- 0
gg$lower[gg$lower > 1] <- 1
gg$upper[gg$upper > 1] <- 1

gg <- gg[order(gg$w), ]
gg$group <- 1:nrow(gg)
gg$Ability <- "Within 1/2 SD"
gg$Ability[abs(gg$theta1 - gg$theta2) > .5] <- "Not Within 1/2 SD"
gg$Ability <- ordered(gg$Ability, c("Within 1/2 SD", "Not Within 1/2 SD"))
table(gg$Ability)

ggplot(gg, aes(x = group, y = w, group = Ability)) +
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

(gg$lower[gg$Ability == "Within 1/2 SD"] > .5) %>% mean
(gg$upper[gg$Ability == "Within 1/2 SD"] < .5) %>% mean

# ------------------------------------------------------------
# Goodness of fit
# ------------------------------------------------------------

n_reps <- 500
n_obs <- length(est$w)
sim <- data_gen(n_reps, est$w, col_parms, est$theta1, est$theta2, NA_pattern = col_form[odd,])
sim_resp <- sim[grep(items, names(sim))]
sim_est <- est_RS2(sim_resp, col_parms, sim$theta1, sim$theta2, SE = "exp", parallel = T)

sim$l <- l_RSC(sim_resp, sim_est$w, col_parms, sim$theta1, sim$theta2)
l <- l_RSC(col_form[odd,], est$w, col_parms, est$theta1, est$theta2)

quant <- tapply(sim$l, sim$pairs, function(x) quantile(x, p = c(.5, .05))) %>% unlist %>% matrix(, nrow = 2, ncol = n_obs) %>% t()

# Visual key for misfit
fit <- rep("<.95", times = n_obs)
fit[l < quant[,2]] <- ">.95"
fit <- ordered(fit, c(">.95", "<.95"))

# Set up and plot
gg <- data.frame(-2*sim$l, -2*rep(l, each = n_reps), rep(fit, each = n_reps), rep(-2*quant[,1], each = n_reps))

names(gg) <- c("l_dist", "l_obs", "Fit", "median")
gg <- gg[order(gg$median), ]
gg$pair <- rep(1:n_obs, each = n_reps)

ggplot(gg, aes(x = pair, y = l_dist, group = pair)) +
  geom_boxplot(outlier.size = 0, outlier.color = "white", aes(color = Fit, fill = Fit)) +
  geom_point(aes(x = pair, y = l_obs, pch = Fit, size = Fit)) +
  scale_shape_manual(values = c(4, 20)) +
  scale_fill_grey(start = .1, end = 0.7) +
  scale_color_grey(start = .2, end = 0.6) +
  xlab("Groups") +
  ylab("-2 * log-likelihood") +
  scale_size_manual(values = c(4, 1)) +
  theme(legend.title=element_blank()) +
  ylim(c(0, 37)) +
  theme_bw(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))









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
