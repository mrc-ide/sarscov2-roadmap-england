source("global_util.R")

version_check("spimalot", "0.2.53")
version_check("sircovid", "0.11.29")

date <- "2021-07-31"

dat <- spimalot::spim_combined_load("regional_results", "england")

dir.create("outputs", FALSE, TRUE)
dir.create("figs", FALSE, TRUE)
dir.create("main_figs", FALSE, TRUE)

saveRDS(dat$data, "outputs/aggregated_data.rds")

saveRDS(dat$rt$england, "regional_results/Rt_england.rds")
saveRDS(dat$rt, "regional_results/Rt_all.rds")

saveRDS(dat$onward, "outputs/combined.rds")
saveRDS(dat$ifr_t, "outputs/ifr_t.rds")

spimalot::spim_pars_pmcmc_save(dat$parameters, "outputs/parameters")

write_png("figs/forest_plot.png",
          width = 2400, height = 1600, res = 200,
          spim_plot_forest(dat, plot_type = "non_betas"))

write_png("figs/forest_plot_betas.png", width = 2400, height = 1600, res = 200,
          spim_plot_forest(dat, plot_type = "betas"))

write_png("figs/data_fits.png", width = 2400 / 5 * 7, height = 1800,
          res = 200,
          spimalot::spim_plot_trajectories(
            dat, sircovid::regions("england"),
            c("deaths_hosp", "deaths_carehomes", "deaths_comm", "icu",
              "general", "hosp", "all_admission"),
            with_forecast = FALSE, add_betas = FALSE))

write_png("figs/serology_euroimmun.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_serology(
            dat, sircovid::regions("england"), 1, 40))

write_png("figs/serology_roche_n.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_serology(
            dat, sircovid::regions("england"), 2, 40))

write_png("figs/pillar2.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_pillar2_positivity(
            dat, sircovid::regions("england"), date_min = as.Date("2020-05-15"),
            ymax = 40))

write_png("figs/react.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_react(
            dat, sircovid::regions("england"), date_min = as.Date("2020-05-15"),
            ymax = 3))

write_png("figs/incidence.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_incidence(
            dat, c(sircovid::regions("england"), "england")))

write_png("figs/incidence_per_1000.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_incidence(
            dat, c(sircovid::regions("england"), "england"), per_1000 = TRUE))

write_png("figs/Rt_all.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"), "Rt_all"))

write_png("figs/Rt_eff_all.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"), "eff_Rt_all"))

write_png("figs/Rt_eff_general.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"),
            "eff_Rt_general"))

write_png("figs/Rt_general.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"), "Rt_general"))

write_png("figs/beta.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"), "beta"))

write_png("figs/IFR_t_all.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_ifr_t(
            dat, c(sircovid::regions("england"), "england"),
            "IFR_t_all"))

write_png("figs/IFR_t_all_no_vacc.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_ifr_t(
            dat, c(sircovid::regions("england"), "england"),
            "IFR_t_all_no_vacc"))

write_png("figs/IFR_t_general.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_ifr_t(
            dat, c(sircovid::regions("england"), "england"),
            "IFR_t_general"))

write_png("figs/IFR_t_general_no_vacc.png", width = 2400, height = 1200,
          res = 200,
          spimalot::spim_plot_ifr_t(
            dat, c(sircovid::regions("england"), "england"),
            "IFR_t_general_no_vacc"))

write_png("figs/ALOS.png", width = 2400, height = 1200,
          res = 200,
          spimalot::spim_plot_alos(
            dat, c(sircovid::regions("england"), "england"), 5, 25))

## plotting admissions demography
mean_admissions <- spimalot::spim_extract_admissions_by_age(dat)
write_csv(mean_admissions, "outputs/mean_admissions_by_age.csv")

write_png("figs/admissions_demo.png", width = 2400, height = 1000, res = 200,
          spimalot::spim_plot_admissions_by_age(dat, "england"))

## add (zoomed in) plots of SPI-M-relevant trajectories
write_png("main_figs/regions.png", width = 2400 / 5 * 7, height = 1800,
          res = 200,
          spimalot::spim_plot_trajectories(
            dat, sircovid::regions("england"),
            c("deaths_hosp", "icu", "general", "all_admission"),
            date_min = as.Date("2020-12-01"),
            with_forecast = FALSE, add_betas = TRUE))

write_png("main_figs/pillar2_positivity.png", width = 2400, height = 1200,
          res = 200,
          spimalot::spim_plot_pillar2_positivity(
            dat, sircovid::regions("england"),
            date_min = as.Date("2020-12-01"),
            ymax = 40, add_betas = TRUE))

write_png("main_figs/prevalence.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_react(
            dat, sircovid::regions("england"),
            date_min = as.Date("2020-12-01"),
            ymax = 3, add_betas = TRUE))

