source("global_util.R")

version_check("spimalot", "0.2.53")
version_check("sircovid", "0.11.29")

model_type <- "BB"
spimalot::spim_check_model_type(model_type)

trim_deaths <- 0
trim_pillar2 <- 0
date <- "2021-07-31"
date_restart <- "2021-03-08"
burnin <- 1000
forecast_days <- 0

control <- spimalot::spim_control(
  short_run, chains, date_restart, n_mcmc = n_mcmc,
  burnin = burnin, forecast_days = forecast_days)

parameters <- 
  list(info = read_csv(paste0("parameters/", assumptions, "/info.csv")),
       prior = read_csv("parameters/prior.csv"),
       proposal = read_csv(paste0("parameters/", assumptions, "/proposal.csv")))
class(parameters) <- "spim_pars_pmcmc"

data_vaccination <- read_csv("data/vaccination.csv")
data_rtm <- read_csv("data/rtm.csv")
data_serology <- read_csv("data/serology.csv")

## Set beta dates aside from spimalot
beta_date <- c("2020-03-16", # PM advises WFH, against non-essential travel etc
               "2020-03-23", # PM announces full lockdown
               "2020-03-25", # Lockdown into full effect
               "2020-05-11", # Initial easing of lockdown
               "2020-06-15", # Non-essential shops can open
               "2020-07-04", # Restaurants, pubs etc can open
               "2020-08-01", # "Eat out to help out" scheme starts
               "2020-09-01", # Schools reopen
               "2020-09-14", # "Rule of six" introduced
               "2020-10-14", # Tiered system introduced
               "2020-10-31", # Lockdown announced
               "2020-11-05", # Lockdown 2 starts
               "2020-12-02", # Lockdown 2 ends
               "2020-12-18", # School Christmas holidays
               "2020-12-25", # Last day of holidays season relaxation
               "2021-01-05", # Lockdown 3 starts
               "2021-03-08", # Step 1: schools reopen
               "2021-04-01", # School holidays
               "2021-04-19", # Step 2: (04-12) + schools return (04-19)
               "2021-05-17", # Step 3: indoors hospitality
               "2021-06-21", # Original planned date of Step 4
               "2021-07-03", # Euro 2020 quarter finals
               "2021-07-11", # Euro 2020 Final
               "2021-07-19"  # Step 4: full lift
               )

parameters_spim <- readRDS("parameters/spim.rds")
efficacy <- parameters_spim$vaccine_efficacy$central
uptake <- parameters_spim$vaccine_uptake[, "central"]
mean_days_between_doses <-
  parameters_spim$vaccine_mean_time_between_doses
end_date <- as.Date(date) + forecast_days
vaccination <- spimalot::spim_vaccination_data(
  date, region, uptake, end_date,
  mean_days_between_doses, efficacy,
  data_vaccination)

waning_rate <- parameters_spim$waning_rate[[assumptions]]

pars <- spimalot::spim_pars(date, region, model_type, multistrain = FALSE,
                            beta_date = beta_date, vaccination = vaccination, 
                            parameters = parameters,
                            kernel_scaling = kernel_scaling,
                            waning_rate = waning_rate)

data <- spimalot::spim_data(date, region, model_type, data_rtm, data_serology,
                            trim_deaths, trim_pillar2, full_data = FALSE)
## This data set includes series we do not fit to
data_full <- spimalot::spim_data(date, region, model_type, data_rtm,
                                 data_serology, trim_deaths, trim_pillar2,
                                 full_data = TRUE)

data_inputs <- spimalot::spim_fit_process_data(NULL, data_rtm, data,
                                               data_full, vaccination)

filter <- spimalot::spim_particle_filter(data, pars, control$particle_filter)
samples <- spimalot::spim_pmcmc(pars, filter, control$pmcmc)

dat <- spimalot::spim_fit_process(samples, parameters, data_inputs,
                                  control$forecast, random_sample = FALSE)
restart <- dat$restart
dat$restart <- NULL

dir.create("outputs", FALSE, TRUE)
saveRDS(dat, "outputs/fit.rds")
saveRDS(restart, "outputs/restart.rds")

message("Creating plots")
write_pdf(
  "outputs/pmcmc_traceplots.pdf",
  spimalot::spim_plot_fit_traces(dat$pmcmc),
  width = 16, height = 9)
