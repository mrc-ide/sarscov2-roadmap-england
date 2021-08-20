source("global_util.R")

version_check("spimalot", "0.2.53")
version_check("sircovid", "0.11.29")

date <- "2021-07-31"
date_restart <- "2021-03-08"
date_new_restart <- c("2021-03-09", "2021-06-21", "2021-07-19")

model_type <- "BB"
spimalot::spim_check_model_type(model_type)
multistrain <- TRUE

trim_deaths <- 0
trim_pillar2 <- 0
burnin <- 1000
forecast_days <- 0

control <- spimalot::spim_control(short_run, chains, date_new_restart,
                                  n_mcmc = n_mcmc, burnin = burnin,
                                  forecast_days = forecast_days)

restart <- spim_restart_load(readRDS("restart.rds"), date_restart)
parent_prior <- read_csv("prior_parent.csv")
parent_prior <- parent_prior[parent_prior$region == region, ]

data_vaccination <- read_csv("data/vaccination.csv")
data_rtm <- read_csv("data/rtm.csv")
data_serology <- read_csv("data/serology.csv")

## 1. Model data.
## Note that the parent fits will have dropped off the restart after running
## the data on date_restart. Hence in these fits, we run against data *after* 
## (and not including) date_restart
data <- spimalot::spim_data(date, region, model_type, data_rtm, data_serology,
                            trim_deaths, trim_pillar2, full_data = FALSE,
                            fit_to_variants = multistrain)
data <- data[data$date > sircovid::sircovid_date(restart$info$date_restart), ]
data_full <- spimalot::spim_data(date, region, model_type, data_rtm,
                                 data_serology, trim_deaths, trim_pillar2,
                                 full_data = TRUE,
                                 fit_to_variants = multistrain)
data_full <- data_full[data_full$date >
                         sircovid::sircovid_date(restart$info$date_restart), ]


## 2. Model parameters.
region <- restart$info$region # same as 'region' parameter
beta_date <- restart$info$beta_date
vaccination <- restart$data$vaccination

stopifnot(restart$info$model_type == model_type)

parameters_spim <- readRDS("parameters/spim.rds")

efficacy <- parameters_spim$vaccine_efficacy$central
uptake <- parameters_spim$vaccine_uptake[, "central"]
mean_days_between_doses <- parameters_spim$vaccine_mean_time_between_doses
end_date <- as.Date(date) + forecast_days
set.seed(1) # this has stochastic tiebreaks
vaccination <- spimalot::spim_vaccination_data(
  date, region, uptake, end_date,
  mean_days_between_doses, efficacy,
  data_vaccination)


VOC <- "delta"
## Create 2 strain efficacy structure:
vaccination$efficacy <- sircovid::modify_severity(
  vaccination$efficacy,
  parameters_spim$strain[[VOC]]$strain_vaccine_efficacy[[assumptions]],
  parameters_spim$strain[[VOC]]$strain_severity_modifier[[assumptions]])
cross_immunity <- 
  parameters_spim$strain[[VOC]]$strain_cross_immunity[[assumptions]]


data_inputs <- spimalot::spim_fit_process_data(NULL, data_rtm, data,
                                               data_full, vaccination)

parameters <- restart$pars
i <- max(which(sircovid::sircovid_date(beta_date) <=
                 sircovid::sircovid_date(restart$info$date_restart)))

## Retain original prior only for first beta on date_restart
beta_restart <- sprintf("beta%d", seq(i, length(beta_date)))
parameters$prior[match(beta_restart[-1], parameters$prior$name), ] <-
  parent_prior[match(beta_restart[-1], parent_prior$name), ]

  
parameters <- spimalot:::spim_add_par(parameters, "strain_transmission_2",
                                      1, 0, 3,
                                      proposal_variance = 0.02,
                                      prior = list(type = "null"))
  
parameters <- spimalot::spim_add_par(parameters, "strain_seed_date",
                                     sircovid::sircovid_date("2021-04-01"),
                                     sircovid::sircovid_date("2021-03-10"),
                                     sircovid::sircovid_date("2021-05-31"),
                                     proposal_variance = 5,
                                     prior = list(type = "null"),
                                     discrete = TRUE)

waning_rate <- parameters_spim$waning_rate[[assumptions]]

## We fix the past betas that have no impact on the restart fits and
## the start_date parameter
i <- max(which(sircovid::sircovid_date(beta_date) <=
                 sircovid::sircovid_date(restart$info$date_restart)))
fixed <- c(sprintf("beta%d", seq_len(i - 1)), "start_date")
if (rerun) {
  parameters <- spim_restart_rerun_pars(region, parameters, fixed, assumptions)
}

pars_full <- spimalot::spim_pars(date, region, model_type, multistrain,
                                 beta_date, vaccination, parameters,
                                 cross_immunity,
                                 kernel_scaling = kernel_scaling,
                                 waning_rate = waning_rate)


pars <- pars_full$fix(pars_full$initial()[fixed])
attr(pars, "inputs") <- attr(pars_full, "inputs")

initial <- spim_restart_initial_inflate_strain(pars, restart, multistrain)

filter <- spim_restart_filter(data, pars, control$particle_filter,
                              initial, date_restart)
samples <- spimalot::spim_pmcmc(pars, filter, control$pmcmc)


dat <- spimalot::spim_fit_process(samples, parameters, data_inputs,
                                  control$forecast, random_sample = FALSE)

new_restart <- dat$restart
dat$restart <- NULL

## add in previous data for the trajectories.
dat <- spim_restart_join_parent(dat, restart$parent, restart$data, date_restart)

dir.create("outputs", FALSE, TRUE)
saveRDS(dat, "outputs/fit.rds")
saveRDS(new_restart, "outputs/restart.rds")

message("Creating plots")
write_pdf(
  "outputs/pmcmc_traceplots.pdf",
  spimalot::spim_plot_fit_traces(dat$pmcmc),
  width = 16, height = 9)
