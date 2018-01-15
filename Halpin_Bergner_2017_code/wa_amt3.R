# ------------------------------------------------------------
# Real data example for paper: Psychometric models of small group collaborations
# Last updated: 8/23/2017
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
library("cirt")

# ------------------------------------------------------------
# Load item parms
# ------------------------------------------------------------
setwd("~/Dropbox/Academic/Projects/CA/Data/response_matrices")
temp_parms <- read.csv("calibration_parms.csv", row.names = 1)

# Drop DIF items
dif_items <- "045|065"
temp_parms <- temp_parms[!grepl(dif_items, row.names(temp_parms)),]

# Item names without versions suffix for easy extraction
items <- paste0(row.names(temp_parms), collapse = "|")

# Individual and collaborative versions
ind_parms <- col_parms <- temp_parms
row.names(ind_parms) <- paste0(row.names(ind_parms), "_IND")
row.names(col_parms) <- paste0(row.names(col_parms), "_COL")

# Final parameter set
parms <- rbind(ind_parms, col_parms)

# ------------------------------------------------------------
# Load response data
# ------------------------------------------------------------
collab <- read.csv("collaboration_2016_0.csv", check.names = F)
col_form <- format_resp(collab, row.names(col_parms), "COL")
ind_form <- format_resp(collab, row.names(ind_parms), "IND")

# Apply conjunctive scoring rule to collaborative form
odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]

# Drop 13 unsuable response patterns (all 1 or all 0)
# (But are they unusable with simultaneous MAP??)
drop_groups <- c(
  collab$group_id[apply(col_form, 1, mean, na.rm = T) %in% c(1,0)],
  collab$group_id[apply(ind_form, 1, mean, na.rm = T) %in% c(1,0)])

col_form <-col_form[!collab$group_id%in%drop_groups,]
ind_form <-ind_form[!collab$group_id%in%drop_groups,]

# Reset odd for dropped items
odd <- seq(1, nrow(col_form), by = 2)

# Final response data
resp <- cbind(ind_form, col_form)
nrow(resp)/2
# ------------------------------------------------------------
# Estimate RSC using logit parameterization
# ------------------------------------------------------------
est <- est_RSC(resp, parms, obs = T, sigma = 3)



# ------------------------------------------------------------
# Goodness of fit
# ------------------------------------------------------------
# Simulate data
n_reps <- 500
n_obs <- nrow(est)
sim <- data_gen(n_reps, est$u, col_parms, est$theta1, est$theta2, NA_pattern = col_form[odd,])

# Re-Estimate model on simulated data
sim_est <- est_RSC2(sim_resp, col_parms, sim$theta1, sim$theta2, sigma = 3, parallel = T)
head(sim_est, 100)

  # Compute quantiles of simulated reference distribtuion
  l_sim <- -2*l_RSC(sim_resp, sim_est$u, col_parms, sim$theta1, sim$theta2)
  l_sim[1:100]
  probs <- c(.025, .05, .10, .5, .80, .93, .975)
  quant <- tapply(l_sim, sim$pairs, function(x) quantile(x, p = probs)) %>% unlist %>% matrix(, nrow = length(probs), ncol = n_obs) %>% t()

  colnames(quant) <- c("p_025", "p_05", "p_10", "p_50", "p_90", "p_95", "p_975")

  # Compute likelihoods from original data
  l_obs <- -2*l_RSC(col_form[odd,], est$u, col_parms, est$theta1, est$theta2)

  # Format data
  gg <- data.frame(l_obs, quant)
  head(gg)

  # Visual indicator for fit
  gg$Fit <- rep("<.95", times = n_obs)
  gg$Fit[gg$l_obs > gg$p_95] <- ">.95"
  gg$Fit <- ordered(gg$Fit, c(">.95", "<.95"))
  poor_fit <- which(gg$Fit == ">.95")

  # Sort
  gg <- gg[order(gg$p_50), ]
  gg$group <- 1:nrow(gg)

  # Plot
  ggplot(gg[, ], aes(x = group, y = l_obs, group = Fit)) +
      geom_errorbar(aes(ymin = p_10, ymax = p_90, color = Fit), size = 2, width = 0) +
      geom_errorbar(aes(ymin = p_05, ymax = p_95, color = Fit), size = .5, width = 0, alpha = 1) +
      geom_point(aes(x = group, y = l_obs, pch = Fit, size = Fit)) +
      scale_shape_manual(values = c(4, 20)) +
      scale_size_manual(values = c(4, 1)) +
      #scale_color_manual(values = c("#132B43", "#56B1F7")) +
      scale_color_manual(values = c("grey10", "grey70")) +
      xlab("Group") +
      ylab("-2 * log-likelihood") +
      theme_bw() +
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
            text = element_text(size=20))   #, legend.title=element_blank()))



# ------------------------------------------------------------
# Estimates of w
# ------------------------------------------------------------

# 95% and 80 %Confidence intervals
gg <- est
q80 <- qnorm(.18)
q95 <- qnorm(.09)

gg$lower_80 <- gg$u + q80*gg$u_se
gg$upper_80 <- gg$u - q80*gg$u_se

gg$lower_95 <- gg$u + q95*gg$u_se
gg$upper_95 <- gg$u - q95*gg$u_se

# Indicator for match on ability
gg$Ability <- "Within 1/2 SD"
gg$Ability[abs(gg$theta1  - gg$theta2) > .5] <- "Not Within 1/2 SD"
gg$Ability <- ordered(gg$Ability, c("Within 1/2 SD", "Not Within 1/2 SD"))

# Drop poor fitting groups ?
gg <- gg[-poor_fit,]

# Sort
gg <- gg[order(gg$u), ]
gg$group <- 1:nrow(gg)
head(gg)
# Plot
ggplot(gg[, ], aes(x = group, y = u, group = Ability)) +
    geom_errorbar(aes(ymin = lower_80, ymax = upper_80, color = Ability), size = 2, width = 0) +
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95, color = Ability), size = .5, width = 0, alpha = 1) +
    geom_point(color = "white", aes(size = Ability)) +
    #scale_color_manual(values = c("#132B43", "#56B1F7")) +
    scale_color_manual(values = c("grey10", "grey70")) +
    scale_shape_discrete(solid = F) +
    ylab("Estimated RSC weight (logit)") +
    xlab("Group")  +
    geom_abline(slope = 0, intercept = 0, col = "black", size = 1) +
    scale_size_manual(values = c(1, 1)) +
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), text = element_text(size=20))



(gg$lower_80[gg$Ability == "Within 1/2 SD"] > 0) %>% mean
(gg$lower_95[gg$Ability == "Within 1/2 SD"] > 0) %>% mean

(gg$upper_80[gg$Ability == "Within 1/2 SD"] < 0) %>% mean
(gg$upper_95[gg$Ability == "Within 1/2 SD"] < 0) %>% mean

(var(est$u) - mean(est$u_se^2))/var(est$u)





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
  scale_shape_manual(values = c(20, 4)) +
  scale_color_manual(values = c("black", "black")) +
  scale_size_manual(values = c(3, 4)) +
  theme_bw(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
