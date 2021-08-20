source("global_util.R")

version_check("spimalot", "0.2.53")

## Load real data
date_death_change <- as.Date("2021-07-31") - 45
agg_data <- read.csv("uk_rtm.csv") %>%
  dplyr::mutate(
    date = as.Date(date),
    deaths_hosp = death3,
    all_admission = phe_admissions,
    deaths_carehomes = dplyr::case_when(
      date < date_death_change ~ as.integer(ons_death_carehome),
      date >= date_death_change ~ as.integer(death_chr)
    ),
    deaths_comm = dplyr::case_when(
      date < date_death_change ~ as.integer(ons_death_noncarehome),
      date >= date_death_change ~ as.integer(death_comm)
    )
  )
deaths <- cbind(agg_data$deaths_hosp, agg_data$deaths_comm, agg_data$deaths_carehomes)
agg_data$deaths <- rowSums(deaths, na.rm = TRUE)
agg_data$deaths[apply(deaths, 1, function(x) all(is.na(x)))] <- NA
agg_data <- agg_data %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(
    deaths = sum(deaths, na.rm = TRUE),
    deaths_hosp = sum(deaths_hosp, na.rm = TRUE),
    all_admission = sum(all_admission, na.rm = TRUE),
    strain_non_variant = sum(n_symp_non_delta_variant, na.rm = TRUE),
    strain_tot = sum(n_symp_delta_variant, na.rm = TRUE) + sum(n_symp_non_delta_variant, na.rm = TRUE)
  )


## Load everything, combine, then clean
summaries <- list(
  readRDS("summary_opt_march.rds"),
  readRDS("summary_opt_june.rds"),
  readRDS("summary_opt_july.rds"),
  readRDS("summary_cen_march.rds"),
  readRDS("summary_cen_june.rds"),
  readRDS("summary_cen_july.rds"),
  readRDS("summary_pes_march.rds"),
  readRDS("summary_pes_june.rds"),
  readRDS("summary_pes_july.rds")
)

summary <- setNames(vector("list", length(summaries[[1]])),
                    names(summaries[[1]]))
for (i in seq_along(summary)) {
  if (names(summary)[[i]] == "population_infected") {
    ## all identical
    summary[[i]] <- summaries[[1]]$population_infected
  } else {
    summary[[i]] <- lapply(summaries, "[[", i) %>% dplyr::bind_rows()
  }
}
rm(summaries)

## Load everything, combine, then clean
fits <- list(
  optimistic = readRDS("summary_restart_opt.rds"),
  central = readRDS("summary_restart_cen.rds"),
  pessimistic = readRDS("summary_restart_pes.rds")
)
summary_fits <- setNames(vector("list", length(fits[[1]])),
                         names(fits[[1]]))
for (i in seq_along(summary_fits)) {
  summary_fits[[i]] <-
      lapply(names(fits), function(x) {
        fits[[x]][[i]] %>% dplyr::mutate(assumptions = x)
      }) %>% dplyr::bind_rows()
}
# save space
summary_fits$n_doses <- NULL
rm(fits)

## For fits combine single strain up to 8 March to restart after 8 March
population <- summary$population_infected$population$england
prop_infected <- summary$population_infected$prop_infected


## Set plot variables
policy_dates <- as.Date(c(
  "2021-03-29", "2021-04-12", "2021-05-17", "2021-06-21",
  "2021-07-19", "2021-07-24", "2021-09-01", "2021-10-23",
  "2021-11-01", "2021-12-18", "2022-01-04", "2022-02-12", "2022-02-21", "2022-04-02",
  "2022-04-19", "2022-05-28", "2022-06-06", "2022-07-23"
))

pop_eligible <- population *  c(rep(0, 3), 2 / 5, rep(1, 15))
prop_eligible <- sum(pop_eligible) / sum(population)


## Light cleaning
# reset cumulative deaths to zero in fits
fits_state <- summary_fits$state %>%
  dplyr::mutate(value = if_else(state == "deaths", 0, value))

summary$population_infected <- NULL
summary <- lapply(summary, clean_summary, TRUE, TRUE, TRUE, TRUE)
summary <- lapply(summary, function(x) {
  x$scenario <- factor(x$scenario,
    levels = c(
      "July-19 [Low R]", "July-19", "July-19 [High R]",
      "JulyGradual [Low R]", "JulyGradual", "JulyGradual [High R]",
      "June-21 [Low R]", "June-21", "June-21 [High R]",
      "JuneGradual [Low R]", "JuneGradual", "JuneGradual [High R]",
      "AlphaOnly [Low R]", "AlphaOnly", "AlphaOnly [High R]"
    )
  )
  x
})
## save some space
summary$summary_state <- NULL
summary_state <- summary$state

summary_fits <- lapply(summary_fits, clean_summary, FALSE, TRUE, FALSE, FALSE)
fits_state <- summary_fits$state


cols <- spimalot::spim_scenario_cols(
  unique(summary_state$scenario)[c(7, 13, setdiff(1:15, c(7, 13)))] %>% forcats::fct_drop(),
  weight = 0.5)

parent_fits <- readRDS("parent_fits.rds")$simulate %>%
  spimalot::spim_simulate_process_output(
    "england", sircovid::regions("england"),
    c("deaths", "infections", "diagnoses_admitted"),
    reset_states = FALSE, rm.rtUK = TRUE) %>%
  spimalot::spim_simulate_create_summary() %>%
  spimalot::tidy_state_one(list(type = "fit")) %>%
  lapply(clean_summary, FALSE, TRUE, FALSE, FALSE)
parent_fits_state <- parent_fits$state

## Render rmd
deaths_to_18_jul <- fits_state %>% dplyr::filter(
  assumptions == "central",
  date == "2021-07-18",
  state == "Total additional COVID-19 deaths") %>%
dplyr::select(`50%`, `2.5%`, `97.5%`)

deaths_to_20_jun <- fits_state %>% dplyr::filter(
  assumptions == "central",
  date == "2021-06-20",
  state == "Total additional COVID-19 deaths") %>%
dplyr::select(`50%`, `2.5%`, `97.5%`)

deaths_21jun_18jul <- deaths_to_18_jul - deaths_to_20_jun

fits_central_combined <- readRDS("fits_central_combined.rds")
rmarkdown::render("paper_numbers.Rmd")
rm(fits_central_combined)

dir.create("outputs", FALSE, TRUE)

peak_states <- summary_state %>%
  dplyr::filter(grepl("Daily|Total", state)) %>%
  dplyr::group_by(analysis, scenario, state) %>%
  dplyr::filter(`50%` == max(`50%`, na.rm = TRUE),
                (grepl("Total", state) & date == as.Date("2022-06-01")) | !grepl("Total", state)) %>%
  dplyr::group_by(analysis, scenario, state, `50%`) %>%
  dplyr::filter(`97.5%` == max(`97.5%`)) %>%
  dplyr::select(analysis, scenario, state, `2.5%`, `50%`, `97.5%`) %>%
    dplyr::mutate(
      `2.5%` = dplyr::if_else(grepl("July", scenario) & state == "Total additional COVID-19 deaths",
                              `2.5%` + deaths_21jun_18jul[["2.5%"]], `2.5%`),
      `50%` = dplyr::if_else(grepl("July", scenario) & state == "Total additional COVID-19 deaths",
                              `50%` + deaths_21jun_18jul[["50%"]], `50%`),
      `97.5%` = dplyr::if_else(grepl("July", scenario) & state == "Total additional COVID-19 deaths",
                              `97.5%` + deaths_21jun_18jul[["97.5%"]], `97.5%`)
    ) %>%
  dplyr::mutate(value = sprintf("%s (%s, %s)", round(`50%` / 100) * 100,
                                round(`2.5%` / 100) * 100, round(`97.5%` / 100) * 100),
                analysis = as.character(analysis),
                scenario = as.character(scenario)) %>%
  dplyr::ungroup() %>%
  dplyr::select(scenario, analysis, state, value) %>%
  dplyr::arrange(scenario, analysis, state)

peak_states %>% write.csv("outputs/peak_states.csv", row.names = FALSE)
peak_states %>% stargazer::stargazer(out = "outputs/peak_states.tex", summary = FALSE, rownames = FALSE)

## Plots and csv
dir.create("figs", FALSE, TRUE)

cat("\nPlotting - Figure 1")

dat <- spimalot::spim_combined_load("dat", "england")

fig1 <- vaccine_figure_1(parent_fits_state, fits_state, summary_state, agg_data,
                         cols, policy_dates, dat)

jpeg("figs/figure1.jpg", width = 31, height = 20, unit = "cm", res = 300)
plot(fig1)
dev.off()

for (i in seq(6)) {
  jpeg(sprintf("figs/figure1_%s.jpg", LETTERS[[i]]), width = 15, height = 15,
       unit = "cm", res = 300)
  plot(fig1[[i]])
  dev.off()
}


cat("\nPlotting - Figure 2")

fig2 <- vaccine_figure_2(summary, summary_fits, parent_fits, dat, agg_data,
                         prop_eligible, population)

jpeg("figs/figure2.jpg",
    width = 31, height = 18, unit = "cm", res = 300)
fig2
dev.off()

jpeg("figs/figure2_A.jpg", width = 31, height = 9, unit = "cm", res = 300)
fig2[[1]]
dev.off()

jpeg("figs/figure2_B.jpg", width = 31, height = 9, unit = "cm", res = 300)
fig2[[2]]
dev.off()


cat("\nPlotting - Figure 3")

fig3 <- vaccine_figure_3(fits_state, summary_state, policy_dates)
jpeg("figs/figure3.jpg", width = 38, height = 21, unit = "cm", res = 300)
fig3
dev.off()

for (j in seq(4)) {
  for (i in seq(3)) {
    jpeg(sprintf("figs/figure3_%s%s.jpg", LETTERS[[j]], i), width = 15, height = 15,
        unit = "cm", res = 300)
    plot(fig3[[j]][[i]])
    dev.off()
  }
}

cat("\nPlotting - Figure SI Trajectories")

jpeg("figs/figureSI_traj.jpg", width = 31, height = 21, unit = "cm", res = 300)
vaccine_figure_SI_traj(fits_state, summary_state, policy_dates) %>% plot()
dev.off()

cat("\nPlotting - Figure SI Vax")

jpeg("figs/figureSI_vax.jpg",  width = 31, height = 16, unit = "cm", res = 300)
vaccine_figure_SI_state(summary$state_by_age, summary_fits$state_by_age) %>%
  plot()
dev.off()

cat("\nPlotting - Figure SI Proportions")

jpeg("figs/figureSI_prop.jpg", width = 24, height = 18, unit = "cm", res = 300)
plot(spimalot::spim_multivariant_rt_plot(dat, as.Date("2021-07-31"), manuscript = TRUE) /
  (vaccine_figure_SI_prop(dat) + plot_layout(nrow = 2, ncol = 4))) + plot_annotation(tag_levels = "A")
dev.off()
