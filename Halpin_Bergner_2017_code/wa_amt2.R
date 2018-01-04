# ------------------------------------------------------------
# Stan code for paper
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")
# devtools::use_data(sim_parms)
# source("~/github/scirt/R/cIRF_functions.R")
# source("~/github/scirt/R/IRF_functions.R")

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("~/github/scirt/R/stan_formatting.R")
# ------------------------------------------------------------
# AMT example data
# ------------------------------------------------------------

setwd("~/Dropbox/Academic/Projects/CA/Data/response_matrices")
temp_parms <- read.csv("calibration_parms.csv", row.names = 1)

# Drop DIF items ?
#dif_items <- "045|065"
#temp_parms <- temp_parms[!grepl(dif_items, row.names(temp_parms)),]

# Item names without versions suffix for easy extraction
items <- paste0(row.names(temp_parms), collapse = "|")

# Individual and collaborative versions
ind_parms <- col_parms <- temp_parms
row.names(ind_parms) <- paste0(row.names(ind_parms), "_IND")
row.names(col_parms) <- paste0(row.names(col_parms), "_COL")

# Final parameter set
parms <- rbind(ind_parms, col_parms)

collab <- read.csv("collaboration_2016_0.csv", check.names = F)
col_form <- format_resp(collab, row.names(col_parms), "COL")
ind_form <- format_resp(collab, row.names(ind_parms), "IND")

# Apply conjunctive scoring rule to collaborative form
odd <- seq(1, nrow(col_form), by = 2)
col_form[odd,] <- col_form[odd+1,] <- col_form[odd,]*col_form[odd+1,]
resp <- cbind(ind_form, col_form)


# ------------------------------------------------------------
# Stan step up
# ------------------------------------------------------------
RSC_logit <- "~/github/scirt/Halpin_Bergner_2017_code/RSC_logit.stan"
amt_data <- format_stan_data(resp, parms)
fit <- stan(file = RSC_logit, data = amt_data, iter = 2000, chains = 4)

# Extract parms
probs <- c(.005, .025, .05, .1, .9, .95, .975, .995, .99)
int_names <- c(paste0(rep(c("lower", "upper"), each = 4), c("_99", "_95", "_90", "_80",  "_80",  "_90",  "_95", "_99")), "one_sided_99")

theta1_hat <- summary(fit, pars = "theta1", probs = probs)$summary
theta2_hat <- summary(fit, pars = "theta2", probs = probs)$summary
u_hat <- summary(fit, pars = "u", probs = probs)$summary
log_lik <- summary(fit, pars = "log_lik", probs = probs)$summary

# ------------------------------------------------------------
# GOF  Plot
# ------------------------------------------------------------

# Set up df
gg <- data.frame(log_lik)
names(gg)[grep("X", names(gg))] <- int_names

gg$obs <- -2*l_RSC(col_form[odd,],
                u_hat[,"mean"],
                col_parms,
                theta1_hat[,"mean"],
                theta1_hat[,"mean"])



# Visual indicator for fit
gg$fit <- rep("<.95", times = K)
gg$fit[gg$obs > gg$one_sided_99] <- ">.95"
gg$fit <- ordered(gg$fit, c(">.95", "<.95"))
poor_fit <- which(gg$fit == ">.95")

# Sort
gg <- gg[order(gg$mean), ]
gg$group <- 1:nrow(gg)
head(gg)
# Plot
y_max <- max(gg$upper_99)

ggplot(gg[, ], aes(x = group, y = mean, group = fit)) +
    geom_errorbar(aes(ymin = lower_80, ymax = upper_80, color = fit), size = 2, width = 0) +
    geom_errorbar(aes(ymin = lower_95, ymax = one_sided_99, color = fit), size = .5, width = 0, alpha = 1) +
    geom_point(aes(x = group, y = obs, pch = fit, size = fit)) +
    scale_shape_manual(values = c(4, 20)) +
    scale_size_manual(values = c(4, 1)) +
    scale_color_manual(values = c("#132B43", "#56B1F7")) +
    #scale_color_manual(values = c("grey10", "grey70")) +
    xlab("Groups") +
    ylab("-2 * log-likelihood") +
    ylim(c(0, y_max)) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
          text = element_text(size=15),
          legend.title=element_blank())


# ------------------------------------------------------------
# Who doesnt fit? large ability DIF, but no apparant pattern on W
# People who used division of labour, or partners logged out
# ------------------------------------------------------------

collab$group_id[odd[poor_fit]]
apply(col_form[odd[poor_fit], ], 1, sum, na.rm = T)
apply(ind_form[odd[poor_fit], ], 1, sum, na.rm = T)
apply(ind_form[odd[poor_fit]+1, ], 1, sum, na.rm = T)

abs(theta1_hat[poor_fit, "mean"] - theta2_hat[poor_fit, "mean"]) %>% mean
abs(theta1_hat[-poor_fit, "mean"] - theta2_hat[-poor_fit, "mean"]) %>% mean
p(u_hat[poor_fit, "mean"])


# ------------------------------------------------------------
# Results plots
# ------------------------------------------------------------

# Set up df
gg <- data.frame(u_hat)
names(gg)[grep("X", names(gg))] <- int_names

# Indicator for match on ability
gg$Ability <- "Within 1/2 SD"
gg$Ability[abs(theta1_hat[,"mean"]  - theta2_hat[,"mean"] ) > .5] <- "Not Within 1/2 SD"
gg$Ability <- ordered(gg$Ability, c("Within 1/2 SD", "Not Within 1/2 SD"))

# Drop poor fitting groups ?
gg <- gg[-poor_fit,]

# Sort
gg <- gg[order(gg$mean), ]
gg$group <- 1:nrow(gg)

# Plot
ggplot(gg[, ], aes(x = group, y = mean, group = Ability)) +
    geom_errorbar(aes(ymin = lower_80, ymax = upper_80, color = Ability), size = 2, width = 0) +
    geom_errorbar(aes(ymin = lower_90, ymax = upper_90, color = Ability), size = .5, width = 0, alpha = 1) +
    geom_point(aes(size = Ability, color = Ability)) +
    scale_color_manual(values = c("#132B43", "#56B1F7")) +
    #scale_color_manual(values = c("grey10", "grey70")) +
    scale_shape_discrete(solid = F) +
    ylab("Estimated RSC weight (logit)") +
    xlab("Group")  +
    geom_abline(slope = 0, intercept = 0, col = "black", size = 1) +
    scale_size_manual(values = c(2, 2)) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), text = element_text(size=15))

# ------------------------------------------------------------
# Other Results stuff
# ------------------------------------------------------------

# Prop strong synergy
table(gg$Ability)
(gg$lower_80[gg$Ability == "Within 1/2 SD"] > 0) %>% mean
(gg$lower_90[gg$Ability == "Within 1/2 SD"] > 0) %>% mean


# Marginal reliability
(var(gg[,"mean"]) - mean(gg[,"sd"]^2))/var(gg[,"mean"])


# Plot zoom in
ggplot(gg[-(1:100),], aes(x = group, y = mean, group = Ability)) +
    geom_errorbar(aes(ymin = lower_80, ymax = upper_80, color = Ability), size = 6, width = 0) +
    geom_errorbar(aes(ymin = lower_90, ymax = upper_90, color = Ability), size = .5, width = 0, alpha = 1) +
    geom_point(aes(size = Ability, color = Ability)) +
    scale_color_manual(values = c("#132B43", "#56B1F7")) +
    #scale_color_manual(values = c("grey10", "grey70")) +
    scale_shape_discrete(solid = F) +
    ylab("Estimated RSC weight (logit)") +
    xlab("Group")  +
    geom_abline(slope = 0, intercept = 0, col = "black", size = 1) +
    scale_size_manual(values = c(2, 2)) +
    #theme(base_size = 15) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), text = element_text(size=15))


# Correlations with ability
ability <- c(
  apply(cbind(theta1_hat[-poor_fit,"mean"], theta2_hat[-poor_fit,"mean"]), 1, min),
  apply(cbind(theta1_hat[-poor_fit,"mean"], theta2_hat[-poor_fit,"mean"]), 1, max),
  theta1_hat[-poor_fit,"mean"] + theta2_hat[-poor_fit,"mean"]/2,
  abs(theta1_hat[-poor_fit,"mean"] - theta2_hat[-poor_fit,"mean"]))

RSC <- rep(u_hat[-poor_fit,"mean"], times = 4)
temp <- c("min", "max", "mean", "|dif|")
temp <- ordered(temp, temp)
Summary <- rep(temp, each = (length(RSC)/4))
gg <- data.frame(ability, RSC, Summary)
head(gg)

p1 <- ggplot(gg, aes(x = RSC, y = ability)) +
    geom_point(size = 2.5, col = "#56B1F7") +
    #geom_point(size = 2.5, col = "grey20") +
    scale_shape_discrete(solid=F) +
    xlab("Estimated RSC weight (logit)") +
    ylab("Individual Ability")  +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), text = element_text(size=15))

p1 + facet_wrap( ~Summary, ncol = 2)


# Some sampling distributions of weights
u_poor_fit <- paste0("u[", poor_fit, "]")
u_random <- paste0("u[", sample.int(K, 25), "]")

stan_hist(fit, pars = u_poor_fit, include = TRUE, col = "white", fill = "grey20") + xlim(-10, 10)

stan_hist(fit, pars = u_random, include = TRUE, col = "white", fill = "grey60") + xlim(-10, 10)
