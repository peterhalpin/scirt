# --------------------------------------------------------------------
# Helper functions to convert data to named list without NAs for Stan
# -------------------------------------------------------------------

df_to_long <- function(df, names = NULL){
  temp <- which(!is.na(df), arr.ind = T)
  j <- temp[,2]
  k <- temp[,1]
  y <- df[temp]
  N <- length(y)
  out <- list(j, k, y, N)
  if (is.null(names)) names <- c("j", "k", "y", "N")
  names(out) <- names
  out
}

format_stan_data <- function(data, parms, sigma, item_names = NULL){

  # RSC_logit.stan parm names
  ind1_parm_names <- c("jj_ind1", "kk_ind1", "y_ind1", "N_ind1")
  ind2_parm_names <- c("jj_ind2", "kk_ind2", "y_ind2", "N_ind2")
  col_parm_names <- c("jj_col", "kk_col", "y_col", "N_col")
  item_parm_names <- c("alpha_ind", "beta_ind", "alpha_col", "beta_col")

  if (is.null(item_names)){item_names <- row.names(parms)}
  parms <- parms[item_names, ]
  ind_names <- row.names(parms)[grep("IND", row.names(parms))]
  col_names <- row.names(parms)[grep("COL", row.names(parms))]

  parms <- list(parms[ind_names, "alpha"], parms[ind_names, "beta"], parms[col_names, "alpha"], parms[col_names, "beta"])

  names(parms) <- item_parm_names

  odd <- seq(1, nrow(data), by = 2)
  ind1 <- df_to_long(data[odd, ind_names], names = ind1_parm_names)
  ind2 <- df_to_long(data[odd+1, ind_names], names = ind2_parm_names)
  col <- df_to_long(data[odd, col_names], names = col_parm_names)
  const <- list(K = length(odd), J_ind = length(ind_names), J_col = length(col_names), sigma = sigma)

  c(const, ind1, ind2, col, parms)
}
