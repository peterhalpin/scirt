# ------------------------------------------------------------
# Data simulation for paper
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")

library("ggplot2")
library("dplyr")
source("~/github/cirt/R/cIRF_functions.R")
source("~/github/cirt/R/IRF_functions.R")
source("~/github/cirt/R/wa_functions.R")

# ------------------------------------------------------------
# Data simulation for weighted addtive model
# ------------------------------------------------------------

# Data generating parameters
J <- 1000 # n respondents
I <- 100 # n items
i <- 20 # n items for short form
K <- J/2 # n groups
e <- .05 # beta prior parm
set.seed(101)

# Individual test
ind_alpha <- runif(I, .65, 2.5)
ind_beta <- sort(rnorm(I, mean = 0, sd = 1.3))
ind_parms <- data.frame(ind_alpha, ind_beta)
ind_names <- paste0("item", 1:I, "_IND")
ind_names_short <- ind_names[sample(I, i)] # short forms
row.names(ind_parms) <- ind_names
names(ind_parms) <- c("alpha", "beta")

# Group test
col_alpha <- runif(I, .65, 2.5)
col_beta <- sort(rnorm(I, mean = 0, sd = 1.3))
col_parms <- data.frame(col_alpha, col_beta)
col_names <- paste0("item", 1:I, "_COL")
col_names_short <- col_names[sample(I, i)] # short form
row.names(col_parms) <- col_names
names(col_parms) <- c("alpha", "beta")
parms <- rbind(col_parms, ind_parms)

# Respondents
odd <- seq(1, J, by = 2)
theta <- rnorm(J)
theta1 <- theta[odd] # odds
theta2 <- theta[odd + 1] # evens

# RSC parameter
w <- rbeta(K, 1 + e, 1 + e)

# Generate group data
col_data <- data_gen(1, w, col_parms, theta1, theta2)
ind_data <- data_gen(1, rep(.5, J), ind_parms, theta, theta)
data <- cbind(col_data[rep(1:K, each = 2),], ind_data[, -c(1:5)])
# ------------------------------------------------------------
# Parameter estimation for 4 conditions
# ------------------------------------------------------------

# long col, long ind
ll_items <- c(col_names, ind_names)
ll_parms <- parms[ll_items, ]
ml_ll <-  est_WA(data[ll_items], ll_parms, method = "ml",  SE = "exp", parallel = F)
map_ll <- est_WA(data[ll_items], ll_parms, method = "map", SE = "exp", parallel = F)


# long col, short ind
ls_items <- c(col_names, ind_names_short)
ls_parms <- parms[ls_items, ]
ml_ls <-  est_WA(data[ls_items], ls_parms, method = "ml",  SE = "exp", parallel = F)
map_ls <- est_WA(data[ls_items], ls_parms, method = "map", SE = "exp", parallel =F)

# short col, long ind
sl_items <- c(col_names_short, ind_names)
sl_parms <- parms[sl_items, ]
ml_sl <-  est_WA(data[sl_items], sl_parms, method = "ml",  SE = "exp", parallel = F)
map_sl <- est_WA(data[sl_items], sl_parms, method = "map", SE = "exp", parallel = F)

# short col, short ind
ss_items <- c(col_names_short, ind_names_short)
ss_parms <- parms[ss_items, ]
ml_ss <-  est_WA(data[ss_items], ss_parms, method = "ml",  SE = "exp", parallel = F)
map_ss <- est_WA(data[ss_items], ss_parms, method = "map", SE = "exp", parallel = F)


# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

gg <- rbind(ml_ll, map_ll, ml_sl, map_sl, ml_ls, map_ls, ml_ss, map_ss)
gg$ind_form <- rep(c("Individual long", "Indvidual short"), each = K*4)
gg$col_form <- rep(c("Group long", "Group short", "Group long", "Group short"), each = K*2)
gg$Method <- rep(rep(c("ML", "MAP"), each = K), times = 4)
gg$dgp_w <- rep(w, times = 8)
gg$bias <- gg$w - gg$dgp_w
gg$rmse <- sqrt((gg$w - gg$dgp_w)^2)

tapply(gg$rmse, INDEX = list(gg$col_form, gg$ind_form), FUN = mean)

gg$sample <- 0
gg$sample[sample(nrow(gg), nrow(gg)/10)] <- 1

p <- ggplot(gg[ ], aes(x = dgp_w, y = w, group = Method)) +
    geom_point(data = gg[gg$sample == 1,], size = 3, aes(shape = Method, color = Method)) +
    geom_smooth(se = F, size = .8, aes(linetype = Method, color = Method)) +
    scale_color_manual(values = c("grey10", "grey10")) +
    scale_shape_discrete(solid=F) +
    xlab("Data generating values") +
    ylab("Estimate")  +
    geom_abline(slope = 1, intercept = 0, col = "grey50", size = 1.2) +
    theme_bw(base_size = 15)


p + facet_grid(ind_form ~ col_form)


q <- ggplot(gg[ ], aes(x = w, y = w_se, group = Method)) +
    geom_point(data = gg[gg$sample == 1,], size = 3, aes(shape = Method, color = Method)) +
    geom_smooth(se = F, size = .8, aes(linetype = Method, color = Method)) +
    scale_color_manual(values = c("grey10", "grey10")) +
    scale_shape_discrete(solid=F) +
    xlab("Estimate") +
    ylab("Standard error")  +
    theme_bw(base_size = 15)

q + facet_grid(ind_form ~ col_form)
