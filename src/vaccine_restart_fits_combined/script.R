source("global_util.R")

version_check("spimalot", "0.2.53")
version_check("sircovid", "0.11.29")

date <- "2021-07-31"
date_restart <- "2021-03-08"

dat <- spimalot::spim_combined_load("regional_results", "england")

dir.create("outputs", FALSE, TRUE)
dir.create("figs", FALSE, TRUE)
dir.create("main_figs", FALSE, TRUE)

saveRDS(dat$data, "outputs/aggregated_data.rds")
saveRDS(dat$rt$england, "regional_results/Rt_england.rds")
saveRDS(dat$rt, "regional_results/Rt_all.rds")
saveRDS(dat$variant_rt$england, "regional_results/multivariant_Rt_england.rds")
saveRDS(dat$variant_rt, "regional_results/multivariant_Rt_all.rds")
saveRDS(dat$ifr_t, "outputs/ifr_t.rds")
saveRDS(dat$onward, "outputs/combined.rds")
spimalot::spim_pars_pmcmc_save(dat$parameters, "outputs/parameters")

onward_2021_03_09 <- combined_onward_restart(dat$onward, "2021-03-09",
                                             "regional_results",
                                             "england")
saveRDS(onward_2021_03_09, "outputs/combined_2021-03-09.rds")
onward_2021_06_21 <- combined_onward_restart(dat$onward, "2021-06-21",
                                             "regional_results",
                                             "england")
saveRDS(onward_2021_06_21, "outputs/combined_2021-06-21.rds")
onward_2021_07_19 <- combined_onward_restart(dat$onward, "2021-07-19",
                                             "regional_results",
                                             "england")
saveRDS(onward_2021_07_19, "outputs/combined_2021-07-19.rds")

## outputs a CSV with current Rt values for variants
write_csv(spimalot::spim_extract_variants_rt(dat, "Rt_general"),
          "outputs/current_rt_multivariant.csv")

write_csv(spimalot::spim_extract_variants_rt(dat, "eff_Rt_general"),
          "outputs/current_Reff_multivariant.csv")

write_png("figs/forest_plot.png",
          width = 2400, height = 1600, res = 200,
          spimalot::spim_plot_forest(dat, sircovid::regions("england"), "all"))

write_png("figs/data_fits.png", width = 2400 / 5 * 7, height = 1800,
          res = 200,
          spimalot::spim_plot_trajectories(
            dat, sircovid::regions("england"),
            c("deaths_hosp", "deaths_carehomes", "deaths_comm", "icu",
              "general", "hosp", "all_admission"),
            with_forecast = FALSE, add_betas = FALSE))

write_png("figs/variant.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_variant(
            dat, sircovid::regions("england"),
            date_min = as.Date(date_restart)))

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
            dat, c(sircovid::regions("england"), "england"),  ymax = 3,
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

write_png("main_figs/positivity.png", width = 2400, height = 1200,
          res = 200,
          spimalot::spim_plot_pillar2_positivity(
            dat, sircovid::regions("england"),
            date_min = as.Date("2020-12-01"),
            ymax = 40, add_betas = TRUE))

write_png("main_figs/eff_Rt_general.png",
          width = 2400, height = 1200, res = 200,
          spimalot::spim_multivariant_rt_plot(dat, date,
                                              last_beta_days_ago = 12,
                                              rt_type = "eff_Rt_general"))

write_png("main_figs/Rt_general.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_multivariant_rt_plot(dat, date,
                                              last_beta_days_ago = 12,
                                              rt_type = "Rt_general"))
