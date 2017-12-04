# ------------------------------------------------------------
# Data simulation for paper
# ------------------------------------------------------------

# devtools::install_github("peterhalpin/cirt")
# library("cirt")
# devtools::use_data(sim_parms)
# source("~/github/cirt/R/cIRF_functions.R")
# source("~/github/cirt/R/IRF_functions.R")

# ------------------------------------------------------------
# Data simulation for weighted addtive model
# ------------------------------------------------------------

# Data generating parameters
n_obs <- 1000 # n respondents
n_items <- 200 # n items
i <- 20 # n items for short form
K <- n_obs/2 # n groups
sigma <- 2 # beta prior parm
set.seed(101)

# Individual test
ind_alpha <- runif(n_items, .65, 2.5)
ind_beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
ind_parms <- data.frame(ind_alpha, ind_beta)
ind_names <- paste0("item", 1:n_items, "_IND")
ind_names_short <- ind_names[sample(n_items, i)] # random short form
row.names(ind_parms) <- ind_names
names(ind_parms) <- c("alpha", "beta")

# Group test
col_alpha <- runif(n_items, .65, 2.5)
col_beta <- sort(rnorm(n_items, mean = 0, sd = 1.3))
col_parms <- data.frame(col_alpha, col_beta)
col_names <- paste0("item", 1:n_items, "_COL")
col_names_short <- col_names[sample(n_items, i)] # random short form
row.names(col_parms) <- col_names
names(col_parms) <- c("alpha", "beta")

# Respondents
odd <- seq(1, n_obs, by = 2)
theta <- rnorm(n_obs)
theta1 <- theta[odd] # odds
theta2 <- theta[odd + 1] # evens

# RSC parameter
w <- rnorm(K, 0, sigma)

# Generate selected group data
col_data <- data_gen(1, w, col_parms, theta1, theta2)
ind_data <- data_gen(1, rep(.5, n_obs), ind_parms, theta, theta)

col_test <- col_data[rep(1:K, each = 2), -c(1:5)]
ind_test <- ind_data[, -c(1:5)]

# ------------------------------------------------------------
# Parameter estimation for 4 conditions
# ------------------------------------------------------------

select_items <- function(n_items, info){
  temp <- apply(info, 1, order, decreasing = T) %>% t
  q <- cbind(rep(1:dim(info)[1], times = n_items), unlist(c(temp[,1:20])))
  out <- info*NA
  out[q] <- 1
  out
}

ind_info1 <- Info(ind_parms, theta1)
ind_info2 <- Info(ind_parms, theta2)
temp1 <- select_items(i, ind_info1)
temp2 <- select_items(i, ind_info2)
ind_selected_items <- ind_test
ind_selected_items[odd, ] <- temp1
ind_selected_items[odd+1, ] <- temp2

col_info <- Info_RSC(w, col_parms, theta1, theta2)
temp <- select_items(i, col_info)
col_selected_items <- col_test
col_selected_items[odd, ] <- temp
col_selected_items[odd+1, ] <- temp

# selected group test, selected individual test
ll_data <- cbind(ind_test*ind_selected_items, col_test*col_selected_items)
ll_parms <- rbind(ind_parms, col_parms)
ml_ll <-  est_RSC(ll_data, ll_parms, method = "ML")
map_ll <- est_RSC(ll_data, ll_parms, method = "MAP")

plot(w, ml_ll$w)
plot(w, map_ll$w)

# random individual test, selected group test
sl_data <- cbind(ind_test[ind_names_short], col_test*col_selected_items)
sl_parms <- rbind(ind_parms[ind_names_short, ], col_parms)
ml_sl <-  est_RSC(sl_data, sl_parms, method = "ML")
map_sl <- est_RSC(sl_data, sl_parms, method = "MAP")

# selected individual test, random group test
ls_data <- cbind(ind_test*ind_selected_items, col_data[col_names_short])
ls_parms <- rbind(ind_parms, col_parms[col_names_short,  ])
ml_ls <-  est_RSC(ls_data, ls_parms, method = "ML")
map_ls <- est_RSC(ls_data, ls_parms, method = "MAP")

# random group test, random individual test
ss_data <- cbind(ind_test[ind_names_short], col_data[col_names_short])
ss_parms <- rbind(ind_parms[ind_names_short, ], col_parms[col_names_short,  ])
ml_ss <-  est_RSC(ss_data, ss_parms, method = "ML")
map_ss <- est_RSC(ss_data, ss_parms, method = "MAP")


# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

gg <- rbind(ml_ll, map_ll, ml_ls, map_ls, ml_ls, map_ls, ml_ss, map_ss)
gg$ind_form <- rep(c("Individual selected", "Indvidual random"), each = K*4)
gg$col_form <- rep(c("Group selected", "Group random"), each = K*2) %>% rep(times = 2)
gg$Method <- rep(c("ML", "MAP"), each = K) %>% rep(times = 4)
gg$dgp_w <- rep(logit(w), times = 8)
gg[1:200,]
gg$sample <- 0
gg$sample[sample(nrow(gg), nrow(gg)/10)] <- 1

gg$w[abs(gg$w) > max(gg$dgp_w)] <- max(gg$dgp_w)

p <- ggplot(gg, aes(x = dgp_w, y = w, group = Method)) +
    geom_point(data = gg[gg$sample == 1,], size = 3, aes(shape = Method, color = Method)) +
    geom_smooth(se = F, size = .8, aes(linetype = Method, color = Method)) +
    scale_color_manual(values = c("grey10", "grey10")) +
    scale_shape_discrete(solid=F) +
    xlab("Data generating values") +
    ylab("Estimate")  +
    geom_abline(slope = 1, intercept = 0, col = "grey50", size = 1.2) +
    theme_bw(base_size = 15)

p + facet_grid(ind_form ~ col_form)


q <- ggplot(gg, aes(x = w, y = w_se, group = Method)) +
    geom_point(data = gg[gg$sample == 1,], size = 3, aes(shape = Method, color = Method)) +
    geom_smooth(se = F, size = .8, aes(linetype = Method, color = Method)) +
    scale_color_manual(values = c("grey10", "grey10")) +
    scale_shape_discrete(solid=F) +
    xlab("Estimate") +
    ylab("Standard error")  +
    theme_bw(base_size = 15)

q + facet_grid(ind_form ~ col_form)
