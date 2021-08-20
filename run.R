main <- function() {
  orderly::orderly_run("vaccine_fits_data")

  kernel_scaling <- 0.2

  regions <- c("north_west",
               "north_east_and_yorkshire",
               "midlands",
               "east_of_england",
               "london",
               "south_west",
               "south_east")
  assumptions <- c("central", "optimistic", "pessimistic")
  for (a in assumptions) {
    for (r in regions) {
      orderly::orderly_run("vaccine_fits_regional",
                           parameters = list(region = r,
                                             assumptions = a,
                                             kernel_scaling = kernel_scaling,
                                             short_run = TRUE),
                           use_draft = TRUE)
    }
    orderly::orderly_run("vaccine_fits_combined",
                         parameters = list(assumptions = a,
                                           kernel_scaling = kernel_scaling,
                                           short_run = TRUE),
                         use_draft = TRUE)
    for (r in regions) {
      orderly::orderly_run("vaccine_restart_fits_regional",
                           parameters = list(region = r,
                                             assumptions = a,
                                             kernel_scaling = kernel_scaling,
                                             rerun = TRUE,
                                             short_run = TRUE),
                           use_draft = TRUE)
    }
    orderly::orderly_run("vaccine_restart_fits_combined",
                         parameters = list(assumptions = a,
                                           kernel_scaling = kernel_scaling,
                                           short_run = TRUE),
                         use_draft = TRUE)
  }

  n_threads <- spimalot::spim_control_cores()
  n_par <- 3
  restart_date <- c("march", "june", "july")
  for (a in assumptions) {
    for (r in restart_date) {
      for (s in c(FALSE, TRUE)) {
        if (s && r != "july") {
          ## Sensitivity only valid for restart date of july
          next
        }
        pars <- list(restart_date = r,
                     assumptions = a,
                     sensitivity = s,
                     kernel_scaling = kernel_scaling,
                     n_par = n_par,
                     n_threads = n_threads,
                     rrq = FALSE,
                     short_run = TRUE)
        orderly::orderly_run("vaccine_simulation",
                             parameters = pars,
                             use_draft = TRUE)
      }
    }
  }

  orderly::orderly_run("vaccine_simulation_plots",
                       parameters = list(kernel_scaling = kernel_scaling,
                                         n_par = n_par,
                                         short_run = TRUE),
                       use_draft = TRUE)
  orderly::orderly_run("vaccine_simulation_plots_sens",
                       parameters = list(kernel_scaling = kernel_scaling,
                                         n_par = n_par,
                                         short_run = TRUE),
                       use_draft = TRUE)
}

main()
