em_parallel <- function(i, sim_data, parms) {
  ind <- sim_data$samples == i
  EM(models,
    sim_data[ind, grep("item", names(sim_data))],
    parms,
    sim_data$theta1[ind],
    sim_data$theta2[ind],
    sorted = T)
}

pv_priors <- function(temp_em){
  n_reps <- length(temp_em)
  temp <- lapply(temp_em, function(x) c(x$prior, x$se^2)) %>% unlist %>% matrix(nrow = n_reps, ncol = n_models*2, byrow = T) %>% data.frame
  mean <- apply(temp, 2, mean, na.rm = T)
  var <- apply(temp, 2, var, na.rm = T)
  pv_se <- sqrt(mean[5:8] + (1 + 1/n_reps) * var[1:4])
  pv_increase <- (1 + 1/n_reps) * var[1:4]/mean[5:8]
  rbind(mean[1:4], sqrt(mean[5:8]), sqrt(var[1:4]), pv_se, pv_increase)
}

pv_posteriors <- function(temp_em, rep_order){
  temp <- lapply(em_pv_b, function(x) x$posterior) %>% {do.call(rbind, .)}
  temp[rep_order, ]
}


boot_em <- function(sim_data, parms, parallel = T) {
  models <- c("Ind", "Min", "Max", "AI")
  n_reps <- max(sim_data$samples)
  n_models <- length(models)
  out <- data.frame(matrix(0, nrow = n_reps, ncol = n_models*2))
  names(out) <- paste0(rep(c("prior", "se"), each = n_models), 1:n_models)

  fun <- function(i) {
    ind <- sim_data$samples == i
    temp <- sim_data[ind, grep("item", names(sim_data))]
    EM(models, temp, parms, sim_data$theta1[ind], sim_data$theta2[ind], sorted = T)
  }

  if (parallel) {
    out <- t(out) # sigh
    out[] <- parallel::mclapply(1:n_reps, fun) %>% unlist
    out <- t(out)
  } else {
    for(i in 1:n_reps) {
      out[i, ] <- fun(i)
      message(cat("Processing replication ", i, " of ", n_reps, "..."),"\r", appendLF = FALSE)
      flush.console()
    }
  }
  out
}


boot_cp <- function(sim_data, mix_prop, parms) {
  n_reps <- max(sim_data$samples)
  models <- c("Ind", "Min", "Max", "AI")
  item_names <- paste(row.names(parms) , collapse = "|")
  resp <- sim_data[grep(item_names, names(sim_data))]
  components <- likelihood(models, resp, parms, sim_data$theta1, sim_data$theta2, sorted = T, Log = F)
  post  <- posterior(components, mix_prop)
  temp_cp <- lapply(1:n_reps, function(x) class_probs(post[sim_data$samples == x,], sim_data$model[sim_data$samples == x]))
  mean <- Reduce(`+`, temp_cp) / n_reps
  se <- Reduce(`+`, lapply(temp_cp, function(x) (x - mean)^2 / (n_reps-1))) %>% sqrt
  out <- list(mean, se)
  names(out) <- c("mean", "se")
  out
}


item_delta_gg <- function(test_length, sim_data, mix_prop, parms, parallel = T){
  n_reps <- max(sim_data$samples)
  models <- c("Ind", "Min", "Max", "AI")
  Models <- ordered(models, models)
  n_models <- length(models)
  deltas <- item_delta(parms, sim_data$theta1, sim_data$theta2) /.25
  delta_order <- apply(deltas, 1, order) %>% t
  resp <- sim_data[grep("item", names(sim_data))]
  n_items <- ncol(resp)
  n_obs <- nrow(resp)

  fun <- function(j){
    screen <- deltas*NA
    i <- j - test_length + 1
    ind <- cbind(rep(1:n_obs, times = test_length), c(delta_order[,i:j]))
    screen[ind] <- 1
    components <- likelihood(models, resp*screen, parms, sim_data$theta1, sim_data$theta2, Log = F)
    post <- posterior(components, mix_prop)

    temp_cp <- lapply(1:n_reps, function(x) class_probs(post[sim_data$sample == x,], sim_data$model[sim_data$sample == x]))

    cp <- lapply(temp_cp, function(x) (diag(x))) %>% unlist %>% matrix(nrow = 4, ncol = n_reps) %>% t

    data.frame("prob" = c(cp),
      "test" = rep(i, n_reps*n_models),
      "delta" = rep(mean(deltas*screen, na.rm = T), n_reps*n_models))
  }

  if (parallel) {
    gg <- parallel::mclapply(test_length:n_items, fun) %>% {do.call(rbind, .)}
  } else {
    gg <- data.frame(prob = NA, test = NA, delta = NA)
    for (j in test_length:n_items) {
      gg <- rbind(gg, fun(j))
      message(cat("Processing test ", j - test_length + 1, " of ", n_items - test_length + 1, "..."),"\r", appendLF = FALSE)
      flush.console()
    }
    gg <- gg[-1, ]
  }
  gg$model <- rep(rep(Models, each = n_reps), n_items - test_length + 1 )
  gg
}
