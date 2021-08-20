spim_restart_load <- function(restart, date) {
  d <- restart
  step <- sircovid::sircovid_date(date)
  i <- match(step, d$state$time)
  if (is.na(i)) {
    pos <- as.character(sircovid::sircovid_date_as_date(d$state$time))
    stop(sprintf("Can't restart at date '%s', must be one of %s",
                 date, paste(pos, collapse = ", ")))
  }
  d$state$state <- d$state$state[, , i, drop = TRUE]
  d$info$date_restart <- date

  d
}


spim_restart_filter <- function(data, pars, control, initial, date_restart,
                                seed = NULL) {
  ## We need to drop all data before the restart date:
  date_restart <- sircovid::sircovid_date(date_restart)
  data <- data[data$date > date_restart, ]

  steps_per_day <- pars$model(pars$initial())$steps_per_day
  initial_step <- date_restart * steps_per_day
  data <- mcstate::particle_filter_data(data, "date", steps_per_day,
                                        1)
  mcstate::particle_filter$new(
    sircovid:::carehomes_particle_filter_data(data),
    sircovid::carehomes,
    control$n_particles,
    if (control$compiled_compare) NULL else sircovid::carehomes_compare,
    sircovid::carehomes_index,
    initial,
    control$n_threads,
    seed)
}


## Create a suitable 'initial' function after inflating strains
spim_restart_initial_inflate_strain <- function(pars, restart, multistrain) {
  if (multistrain) {
    info <- sircovid::carehomes$new(pars$model(pars$initial()), 0, 1)$info()

    ## Add empty strains to the state:
    state <- sircovid::inflate_state_strains(
      restart$state$state, restart$info$info, info)

  } else {
    state <- restart$state$state
  }

  time <- sircovid::sircovid_date(restart$info$date_restart)

  spimalot::spim_restart_initial(state, time, multistrain)
}


spim_restart_join_parent <- function(fit, parent, data, restart_date) {
  ## First, fix first step; see
  ## https://github.com/mrc-ide/mcstate/issues/55
  fit$samples$trajectories$step <- fit$samples$trajectories$step[-1L]
  fit$samples$trajectories$date <- fit$samples$trajectories$date[-1L]
  fit$samples$trajectories$predicted <- fit$samples$trajectories$predicted[-1L]
  fit$samples$trajectories$state <- fit$samples$trajectories$state[, , -1L]

  ## Then stitch together. Work out what we keep from the parent fit:
  restart_date <- sircovid::sircovid_date(restart_date)
  i <- which(parent$trajectories$date <= restart_date)

  ## Filter trajectories:
  for (v in c("step", "date", "predicted")) {
    fit$samples$trajectories[[v]] <- c(
      parent$trajectories[[v]][i], fit$samples$trajectories[[v]])
  }
  fit$samples$trajectories$state <- mcstate::array_bind(
    parent$trajectories$state[, , i, drop = FALSE],
    fit$samples$trajectories$state)

  join_rt <- function(parent, new, i) {
    ret <- Map(function(a, b) rbind(a[i, , drop = FALSE], b[-1L, ]),
               parent, new)
    class(ret) <- class(new)
    ret
  }

  fit$rt <- join_rt(parent$rt, fit$rt, i)
  fit$ifr_t <- join_rt(parent$ifr_t, fit$ifr_t, i)


  i_data <- which(data$full$date <= restart_date)
  fit$data$full <- rbind(data$full[i_data, ], fit$data$full)
  fit$data$fitted <- rbind(data$fitted[i_data, ], fit$data$fitted)

  ## Extra metadata for the fit:
  fit$samples$info$restart <- TRUE
  fit$samples$info$restart_date <- restart_date
  i_predicted <- which(fit$samples$trajectories$predicted)
  if (length(i_predicted) == 0) {
    i_restart <- seq(max(i) + 1, length(fit$samples$trajectories$predicted))
  } else {
    i_restart <- seq(max(i) + 1, min(i_predicted) - 1)
  }
  fit$samples$info$time_index <- list(
    parent = i,
    restart = i_restart,
    predicted = i_predicted)

  fit
}


spim_restart_rerun_pars <- function(region, parameters, fixed, assumptions) {
  rerun_info <- read_csv(paste0("parameters/", 
                                tolower(assumptions), "/info.csv"))
  rerun_proposal <- read_csv(paste0("parameters/", 
                                    tolower(assumptions), "/proposal.csv"))
  
  rerun_info <- rerun_info[rerun_info$region == region, ]
  match_info <- match(rerun_info$name, parameters$info$name)
  parameters$info[match_info, ] <- rerun_info
  
  rerun_proposal <- rerun_proposal[rerun_proposal$region == region, ]
  match_proposal <- match(rerun_proposal$name, parameters$proposal$name)
  parameters$proposal[match_proposal, c(1, 2, match_proposal + 2)] <-
    rerun_proposal
  
  parameters$proposal[parameters$proposal$name %in% fixed, -c(1, 2)] <- 0
  parameters$proposal[, fixed] <- 0
  
  parameters
}
