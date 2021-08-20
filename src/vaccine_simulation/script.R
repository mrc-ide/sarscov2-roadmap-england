source("global_util.R")

if (rrq) {
  version_check("rrq", "0.4.4")
}
version_check("spimalot", "0.2.57")

stopifnot(restart_date %in% c("march", "june", "july"))

if (sensitivity && restart_date != "july") {
  stop("restart_date!='july' and sensitivity=TRUE cannot be run together")
}

parameters <- readRDS("spim_parameters.rds")
parameters <- c(parameters, parameters$strain[["delta"]])
s3_sd <- 0.25
s4_sd <- 0.20
schools_modifier <- 0.25
end_date <- "2022-06-01"
regions <- sircovid::regions("england")

combined <- readRDS(sprintf("combined_%s.rds", restart_date))
combined_simulate <- combined$simulate

june_date <- "2021-06-21"
july_date <- "2021-07-19"
date <- combined$date

## get fitted Rt values from final date
combined_final <- readRDS("combined_final.rds")
sircovid_restart_dates <- combined_final$simulate$date
restart_dates <- sircovid::sircovid_date_as_date(sircovid_restart_dates)
england_restart_rt <-
  apply(combined_final$simulate$Rt_general[, "england", , ], c(2, 3), mean)

eps <- england_restart_rt[, "strain_2"] / england_restart_rt[, "strain_1"]

rm(combined_final)

#####################################################
## From here only work with 'combined'
#####################################################
population <- spimalot::spim_population(combined, ignore_uk = TRUE)
prop_infected <- spimalot::spim_prop_infected(combined, population)

schools_schedule <- "schools_schedule.csv"

june_july <- seq(as.Date("2021-06-01"), as.Date("2021-07-31"), 1)
beta_mult <- spimalot:::calc_seasonality(
  sircovid::sircovid_date(june_july),
  sircovid::sircovid_date("2020-02-15"), 0.1)

avg_rt <- colMeans(england_restart_rt[restart_dates %in% june_july, ] /
                   beta_mult)

alpha_rt_s3 <- avg_rt[["strain_1"]]
alpha_rt_s4 <- alpha_rt_s3 * parameters$prop_full_lift
delta_rt_s3 <- avg_rt[["strain_2"]]
delta_rt_s4 <- delta_rt_s3 * parameters$prop_full_lift

beta_mult <- spimalot:::calc_seasonality(
  sircovid::sircovid_date("2021-06-21"),
  sircovid::sircovid_date("2020-02-15"), 0.1)
alpha_rt_s3_june <- unname(
  england_restart_rt[restart_dates == "2021-06-21", "strain_1"] / beta_mult
)
beta_mult <- spimalot:::calc_seasonality(
  sircovid::sircovid_date("2021-07-19"),
  sircovid::sircovid_date("2020-02-15"), 0.1)
alpha_rt_s3_july <- unname(
  england_restart_rt[restart_dates == "2021-07-19", "strain_1"] / beta_mult
)

data.frame(
  assumption = assumptions,
  alpha_S3 = round(alpha_rt_s3, 2),
  alpha_S4 = sprintf("{%s}",
                       paste0(round(rev(alpha_rt_s4), 2), collapse = ",")),
  delta_S3 = round(delta_rt_s3, 2),
  delta_S4 = sprintf("{%s}",
                       paste0(round(rev(delta_rt_s4), 2), collapse = ","))
) %>% write.csv("rt_values.csv", row.names = FALSE)

rm(delta_rt_s3, delta_rt_s4)


if (restart_date == "july") {
  stopifnot(date == july_date)
  ## Scenarios:
  ## 1. Unlock July
  ## 2. Unlock July gradual

  ## Create npis for Step 3
  npi_step3 <- data.frame(
    nation = "england", npi = "step_3", Rt = alpha_rt_s3_july,
    Rt_sd = s3_sd, adherence = c("low", "central", "high")
  )

  ## Create npis for July Unlock
  npi_july_unlock <- data.frame(nation = "england", npi = "full_lift",
                                Rt = alpha_rt_s4, Rt_sd = s4_sd,
                                adherence = c("low", "central", "high"))

  ## Create NPI key with gradual scenario (19 Jul is schools open)
  npis <- rbind(npi_step3, npi_july_unlock)
  npi_key <- spimalot::spim_prepare_npi_key(
    "open", schools_modifier, "england", npi_key = npis,
    gradual_start = "step_3_schools_open",
    gradual_end = "full_lift_schools_open",
    gradual_steps = 11
  )

  ## Create lifting schedule for gradual and non-gradual
  july19_schedule <- data.frame(rbind(
    c(date = "2021-05-17", npi = "step_3"),
    c(date = july_date, npi = "full_lift")
  ), nation = "england")

   gradual_schedule <- data.frame(
     date = round(c(as.Date("2021-05-17"),
              seq(as.Date("2021-07-19"), as.Date("2021-10-01"),
                  length.out = 11))),
     npi = c("step_3", sprintf("p%s_full_lift", 1:10), "full_lift"),
     nation = "england")

  ## Create Rt schedule for two scenarios
  rt_schedule <- rbind(
    spimalot::spim_create_rt_scenario(schools_schedule, "england", "July-19",
                                      july19_schedule),
    spimalot::spim_create_rt_scenario(schools_schedule, "england", "JulyGradual",
                                      gradual_schedule)
  )

  ## Get strain transmission for later
  strain_transmission <- sprintf(
    "%dperc", 5 * round(eps[restart_dates == date] * 20)
  )
  rm(npi_step3, npis)
} else {

  ## Create npis for June unlock with and w/o VOC
  npi_june_unlock <- data.frame(nation = "england", npi = "full_lift",
                                Rt = alpha_rt_s4,
                                Rt_sd = s4_sd,
                                adherence = c("low", "central", "high"))

  june21_schedule <- data.frame(date = june_date, npi = "full_lift",
                                nation = "england")

  ## Get strain transmission for later
  strain_transmission <- sprintf(
    "%dperc", 5 * round(eps[restart_dates == june_date] * 20)
  )

  if (restart_date == "june") {
    stopifnot(date == june_date)

    ## Create npis for Step 3
    npi_step3 <- data.frame(
      nation = "england", npi = "step_3", Rt = alpha_rt_s3_june,
      Rt_sd = s3_sd, adherence = c("low", "central", "high")
    )

    ## Scenarios:
    ##  1. Unlock on June 21 with VOC
    ##  2. Gradual unlock on June 21 with VOC
    npis <- rbind(npi_step3, npi_june_unlock)

    npi_key <- spimalot::spim_prepare_npi_key(
      "open", schools_modifier, "england", npi_key = npis,
      gradual_start = "step_3_schools_open",
      gradual_end = "full_lift_schools_open",
      gradual_steps = 11
   )

   gradual_schedule <- data.frame(
     date = round(c(as.Date("2021-05-17"),
            seq(as.Date("2021-06-21"), as.Date("2021-09-01"),
                                  length.out = 11))),
     npi = c("step_3", sprintf("p%s_full_lift", 1:10), "full_lift"),
     nation = "england")

    rt_schedule <- rbind(
      spimalot::spim_create_rt_scenario(schools_schedule, "england", "June-21",
                                        june21_schedule),
      spimalot::spim_create_rt_scenario(schools_schedule, "england", "JuneGradual",
                                        gradual_schedule)
    )
  } else {
    ## Scenarios:
    ##  1. Unlock on June 21 without VOC

    ## Pull the true Rt values from the combined object
    rt_fits_dates <- seq(as.Date("2021-03-09"), as.Date("2021-06-20"), 1)
    rt <- england_restart_rt[, "both"]
    rt <- rt[sircovid_restart_dates %in% sircovid::sircovid_date(rt_fits_dates)]
    beta_mult <- spimalot:::calc_seasonality(
      sircovid::sircovid_date(rt_fits_dates),
      sircovid::sircovid_date("2020-02-15"), 0.1)
    rt <- rt / beta_mult

    ## Create a schedule of daily Rt values
    fits_schedule <- data.frame(date = rt_fits_dates,
                                npi = paste0("fits_rt_", rt_fits_dates),
                                nation = "england")
    fits_schedule <- rbind(fits_schedule,
                          data.frame(date = june_date,
                                      npi = "full_lift",
                                      nation = "england"))

    ## Manually create NPI key with central adherence only
    npi_fits <- data.frame(
      nation = "england",
      npi = paste0("fits_rt_", rt_fits_dates),
      Rt = rt, Rt_sd = s3_sd, adherence = "central")
    npi_fits <-
      rbind(npi_fits,
            npi_fits %>% mutate(adherence = "low"),
            npi_fits %>% mutate(adherence = "high"))

    ## Create NPI key for June unlock (21 Jun is schools open)
    june_npi_key <- spimalot::spim_prepare_npi_key(
      "open", schools_modifier, "england",
      npi_key = npi_june_unlock
    )

    ## Bind fitted Rts to June unlock key
    npi_key <- rbind(june_npi_key, npi_fits)

    ## Create Rt schedule
    rt_schedule <- spimalot::spim_create_rt_scenario(
      schools_schedule, "england", "AlphaOnly", fits_schedule)
    ## Remove schools open/closed suffix for fitted Rts
    which <- grepl("fits_rt_", rt_schedule$npi)
    rt_schedule$npi[which] <- gsub("_schools_open", "", rt_schedule$npi[which])

    rm(npi_june_unlock, june_npi_key, npi_fits)
  }
}

## Clean up
rm(england_restart_rt, sircovid_restart_dates, restart_dates, eps,
   s3_sd, s4_sd, schools_modifier, schools_schedule, alpha_rt_s3_july,
   alpha_rt_s3_june, beta_mult)
#####################################################
## From here everything identical for VOC and no VOC
#####################################################

rt_future <- spimalot::spim_prepare_rt_future(npi_key, date, end_date,
                                              schedule = rt_schedule)
write.csv(rt_future, "Rt_future.csv")

combined <- spimalot::spim_simulate_prepare(combined, n_par,
                                            regions,
                                            inflate_strain = FALSE,
                                            inflate_booster = TRUE)

if (sensitivity) {
  expand_grid <- spimalot::spim_expand_grid(
    adherence_to_baseline_npis = c("low", "central", "high"),
    strain_transmission = strain_transmission,
    strain_vaccine_efficacy = c("pessimistic", "central", "optimistic"),
    strain_cross_immunity = c("pessimistic", "central", "optimistic"),
    waning_rate = c("pessimistic", "central", "optimistic")
  )
} else {
  expand_grid <- spimalot::spim_expand_grid(
    adherence_to_baseline_npis = c("low", "central", "high"),
    strain_transmission = strain_transmission,
    strain_vaccine_efficacy = assumptions
  )
}

run_grid <- spimalot::spim_run_grid(
  unique(rt_future$scenario),
  expand_grid = expand_grid,
  force_central = FALSE,
  set_strain_params = !sensitivity ## match WR and CI to VE if not sensitivity
)
rm(expand_grid, strain_transmission)

  if (sensitivity) {
  # still match severity
  # remove june 21
  run_grid <- run_grid %>%
    dplyr::mutate(strain_severity_modifier = strain_vaccine_efficacy)
} else {
  run_grid <- run_grid %>%
    mutate(
      analysis = sprintf(
        "VOC %s VE%s",
        case_when(
          strain_vaccine_efficacy == "pessimistic" ~ "Low",
          strain_vaccine_efficacy == "optimistic" ~ "High",
          strain_vaccine_efficacy == "central" ~ "Central",
        ),
        case_when(
          adherence_to_baseline_npis == "high" ~ " - low R after full lift",
          adherence_to_baseline_npis == "low" ~ " - high R after full lift",
          TRUE ~ ""
        )
      )
    )
}

#####################################################
## From here everything identical
#####################################################


if ("voc" %in% run_grid$strain_initial_proportion) {
  stop("`strain_initial_proportion` must be 'no_voc'")
}

write_csv(run_grid, "run_grid_final.csv")

ignore <- c("RUN", "analysis", "scenario", "adherence_to_baseline_npis",
            "full_run")

parameters$rt_future <- split(rt_future, rt_future$key)

base <- list(
  end_date = as.Date(end_date),
  seed = NULL,
  n_threads = n_threads,
  output_keep = c("deaths", "admitted", "diagnoses", "infections", "hosp",
                  "icu", "deaths_hosp", "sero_pos_1", "sero_pos_2"),
  output_rt = TRUE,
  output_time_series = TRUE,
  output_vaccination = TRUE,
  output_state_by_age = TRUE,
  output_weight_rt = TRUE,
  rt_type = "Rt_general",
  vaccine_delay_multiplier = 1,
  strain_seed_date = parameters$strain_seed_date,
  voc_seeded = restart_date != "march" ## no VOC seeded in March
)

vars_names <- c("seasonality",
                "rt_future",
                "vaccine_daily_doses",
                "vaccine_booster_daily_doses",
                "vaccine_efficacy",
                "vaccine_booster_efficacy",
                "vaccine_eligibility",
                "vaccine_uptake",
                "vaccine_lag_groups",
                "vaccine_lag_days",
                "strain_transmission",
                "strain_seed_rate",
                "strain_vaccine_efficacy",
                "strain_initial_proportion",
                "strain_vaccine_booster_efficacy",
                "strain_cross_immunity",
                "strain_severity_modifier",
                "waning_rate")
vars <- parameters[vars_names]

args <- spimalot::spim_simulate_args(run_grid, vars, base, ignore,
                                     regions, TRUE)

if (rrq) {
  ## will error if rrq parameter is TRUE but rrq worker not set-up
  message(sprintf("Running %s scenarios with rrq", nrow(run_grid)))
  res <- spimalot::spim_simulate_rrq(args, combined,
                                     spimalot::spim_rrq_controller())
} else {
  message(sprintf("Running %s scenarios locally", nrow(run_grid)))
  res <- spimalot::spim_simulate_local(args, combined)
}

dir.create("outputs", FALSE, TRUE)

eng_regions <- sircovid::regions("england")
incidence_states <- c("deaths", "infections", "diagnoses_admitted", "deaths_hosp")

if (restart_date == "march") {
  reset_date <- "2021-06-20"
} else {
  reset_date <- NULL
}
res_england <- lapply(res, spimalot::spim_simulate_process_output, "england",
                      eng_regions, incidence_states, reset_states = TRUE,
                      remove_dates_to = reset_date)

message("Creating summary")
summary_england <- lapply(res_england, spimalot::spim_simulate_create_summary)

## From here we work with "tidy" (long) versions ready for plotting
summary_tidy <- spimalot::spim_simulate_tidy_states(summary_england, run_grid,
                                                    combined)
summary_tidy$population_infected <- list(population = population,
                                         prop_infected = prop_infected)
saveRDS(summary_tidy, "outputs/summary.rds")

## do the same with VOC combined fits
summary_restart_tidy <-
  spimalot::spim_simulate_process_output(combined_simulate, "england",
                                        eng_regions, incidence_states,
                                        reset_states = FALSE, rm.rtUK = TRUE,
                                        remove_dates_to = reset_date) %>%
  spimalot::spim_simulate_create_summary() %>%
  spimalot::tidy_state_one(list(type = "fit"))

saveRDS(summary_restart_tidy, "outputs/summary_restart.rds")

## do the same with no voc fits
complete_second_doses <- dplyr::filter(summary_tidy$n_doses,
                                      state == "second_dose_inc",
                                      region == "england") %>%
  dplyr::group_by(across(-c(value, group))) %>%
  dplyr::summarise(mean = sum(value)) %>%
  dplyr::filter(mean < 5e3) %>%
  dplyr::filter(date == min(date))

write.csv(complete_second_doses, "outputs/complete_second_doses.csv")


dir.create("figs", FALSE, TRUE)

g <- summary_tidy$state %>%
  dplyr::filter(region == "england") %>%
  dplyr::mutate(scenario = as.factor(scenario),
                analysis  = as.factor(analysis))

g_restart <- summary_restart_tidy$state %>%
  dplyr::filter(region == "england")

if (restart_date != "march") {
  png("figs/check_Rt.png", width = 2400, height = 1000, res = 200)
  f <- function(g) {
    g %>%
      dplyr::filter(state %in% c("Rt_general_both", "eff_Rt_general_both")) %>%
      tidyr::pivot_wider(names_from = quantile)
  }

  p <- f(g) %>%
    ggplot(aes(x = date, y = `50%`, colour = scenario)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_ribbon(alpha = 0.3,
                aes(ymin = `2.5%`,
                    ymax = `97.5%`,
                    fill = scenario)) +
    geom_line() +

    geom_ribbon(data = f(g_restart), alpha = 0.3,
                aes(x = date,
                    ymin = `2.5%`,
                    ymax = `97.5%`), inherit.aes = FALSE) +
    geom_line(data = f(g_restart), aes(x = date, y = `50%`), inherit.aes = FALSE) +
    facet_grid(cols = vars(analysis), rows = vars(state),  labeller = label_wrap_gen(width=7)) +
    scale_x_date(date_breaks = "1 month") +
    geom_vline(xintercept = rt_schedule$date, lty = 2, color = "gray")
  plot(p)
  dev.off()
} else {
  g_restart <- summary_restart_tidy$state %>%
    dplyr::filter(region == "england")

  png("figs/check_Rt.png", width = 2400, height = 1000, res = 200)
  f <- function(g) {
    g$state <- gsub("_both", "", g$state)
    g %>%
      dplyr::filter(state %in% c("Rt_general", "eff_Rt_general")) %>%
      tidyr::pivot_wider(names_from = quantile)
  }

  p <- f(g) %>%
    ggplot(aes(x = date, y = `50%`, colour = scenario)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_ribbon(alpha = 0.3, color = "black",
                aes(ymin = `2.5%`,
                    ymax = `97.5%`,
                    fill = scenario)) +
    geom_line(color = "black") +

    geom_ribbon(data = f(g_restart), alpha = 0.3,
                aes(x = date,
                    ymin = `2.5%`,
                    ymax = `97.5%`), inherit.aes = FALSE, color = "red") +
    geom_line(data = f(g_restart), aes(x = date, y = `50%`),
              inherit.aes = FALSE, color = "red") +

    facet_grid(cols = vars(analysis), rows = vars(state),  labeller = label_wrap_gen(width=7)) +
    scale_x_date(date_breaks = "1 month") +
    geom_vline(xintercept = rt_schedule$date, lty = 2, color = "gray")
  plot(p)
  dev.off()
}

png("figs/check_state.png", width = 2400, height = 1000, res = 200)

  f <- function(g) {
    g %>%
      dplyr::filter(state %in% c("diagnoses_admitted_inc",
                                "deaths_inc", "deaths_hosp_inc",
                                "hosp"),
                    group == "all") %>%
      tidyr::pivot_wider(names_from = quantile)
  }

  p <- f(g) %>%
    ggplot(aes(x = date)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_ribbon(alpha = 0.3,
                aes(ymin = `2.5%`,
                    ymax = `97.5%`,
                    fill = scenario)) +

    geom_line(aes(y = `50%`, color = scenario)) +
      geom_ribbon(data = f(g_restart), alpha = 0.3,
                  aes(x = date,
                      ymin = `2.5%`,
                      ymax = `97.5%`), inherit.aes = FALSE) +
      geom_line(data = f(g_restart), aes(x = date, y = `50%`), inherit.aes = FALSE) +
    facet_grid(rows = vars(state), cols = vars(analysis), scales = "free",
              labeller = label_wrap_gen(width=5))
  plot(p)
  dev.off()

if (restart_date == "july" && !sensitivity) {
  central_analysis <- unique(run_grid$analysis)[2]
  png("figs/check_state_by_age.png", width = 2400, height = 1000, res = 200)
    p <- summary_tidy$state_by_age %>%
      dplyr::filter(region == "england",
                    analysis == central_analysis,
                    scenario == "July-19",
                    state %in% c("infections_inc", "diagnoses_admitted_inc",
                                 "deaths_inc", "deaths")) %>%
      ggplot(aes(x = date, y = value, fill = vaccine_status)) +
      geom_area() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_grid(vars(state), vars(group), scales = "free", labeller = label_wrap_gen(width=10))
    plot(p)
  dev.off()

  doses_g <- summary_tidy$n_doses %>%
    dplyr::filter(region == "england",
                  scenario == "July-19") %>%
    tidyr::pivot_wider(names_from = state, names_prefix = "state_") %>%
    dplyr::mutate(state_total_dose_inc = state_first_dose_inc +
                    state_second_dose_inc +
                    state_booster_dose_inc) %>%
    tidyr::pivot_longer(starts_with("state_"), names_to = "state") %>%
    tidyr::pivot_wider(names_from = group, names_prefix = "group_") %>%
    dplyr::mutate(group_total = rowSums(dplyr::across(starts_with("group")),
                                        na.rm = TRUE)) %>%
    tidyr::pivot_longer(starts_with("group_"), names_to = "group")

  pop_df <- data.frame(group = unique(doses_g$group),
                      pop = c(population$england,
                              total = sum(population$england)))
  doses_g <- doses_g %>%
    dplyr::left_join(pop_df) %>%
    dplyr::mutate(prop = value / pop)

  png("figs/check_doses.png", width = 2400, height = 2000, res = 200)
  p <- doses_g %>%
    ggplot(aes(x = date, y = value, colour = group)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_line() +
    facet_grid(rows = vars(state), cols = vars(analysis), scales = "free", labeller = label_wrap_gen(width=10))
  plot(p)
  dev.off()

  png("figs/check_total_doses.png", width = 2400, height = 1000, res = 200)
  p <- doses_g %>%
    dplyr::filter(group == "group_total",
                  state == "state_total_dose_inc",
                  region == "england") %>%
    ggplot(aes(x = date, y = value * 7, colour = group)) +
    theme_bw() +
    ylab("Weekly doses") +
    geom_line() +
    facet_wrap(vars(analysis))
  plot(p)
  dev.off()

  png("figs/check_uptake.png", width = 2400, height = 1000, res = 200)
  p <- doses_g %>%
    dplyr::filter(state %in% c("state_first_dose", "state_second_dose",
                              "state_booster_dose")) %>%
    ggplot(aes(x = date, y = prop, colour = group)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_line() +
      facet_grid(rows = vars(state), cols = vars(analysis), scales = "free", labeller = label_wrap_gen(width = 10))
  plot(p)
  dev.off()

report_date <- date
doses_g %>%
  dplyr::filter(
    state %in% c("state_first_dose", "state_second_dose"),
    date == report_date
  ) %>%
  dplyr::select(group, state, analysis, prop) %>%
  arrange(state, analysis, group, prop) %>%
  write.csv("outputs/uptake_simulation_date.csv", row.names = FALSE)

doses_g %>%
  dplyr::filter(
    state %in% c("state_first_dose", "state_second_dose"),
    date == "2022-06-01"
  ) %>%
  dplyr::select(group, state, analysis, prop) %>%
  arrange(state, analysis, group, prop) %>%
  write.csv("outputs/uptake_2022-06-01.csv", row.names = FALSE)
} else {
  write.csv("voc=FALSE | sensitivity=TRUE",
            "outputs/uptake_simulation_date.csv")
  write.csv("voc=FALSE | sensitivity=TRUE", "outputs/uptake_2022-06-01.csv")
}

## plot Rt distribution
rt_dist_labels <- c("July-19",
                    "July-19 [High R]",
                    "July-19 [Low R]",
                    "Schools open",
                    "Schools closed")

npi_pars <- lapply(unique(npi_key$nation), function(i) {
  key <- npi_key %>% dplyr::filter(nation == i) %>% dplyr::select(-nation)
  npi_pars <- mapply( function(mean, sd) {
    dist <- distr6::dstr("Lognormal", mean = mean, sd = sd)
    unlist(c(q2.5 = dist$quantile(0.025),
            q97.5 = dist$quantile(0.975),
            meanlog = dist$parameters("meanlog")$value,
            sdlog = dist$parameters("sdlog")$value))
  }, mean = key$Rt, sd = key$Rt_sd)
  npi_pars <- cbind(key, signif(t(npi_pars), 3))
  npi_pars$region <- i
  npi_pars
})
npi_pars <- suppressWarnings(dplyr::bind_rows(npi_pars))
write.csv(npi_pars, "outputs/npi_pars.csv")

