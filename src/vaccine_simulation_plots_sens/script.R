states <- c("Total additional COVID-19 deaths", "Total hospital admissions",
            "Total infections")
dir.create("figs", FALSE, TRUE)

for (i in c("cen", "opt", "pes")) {
  message(sprintf("Loading summary (%s)", i))
  summary <- readRDS(sprintf("summary_%s.rds", i))$summary_state

  message(sprintf("Cleaning summary (%s)", i))
  summary <- clean_summary(summary, TRUE, TRUE, TRUE, TRUE) %>%
    dplyr::filter(state %in% states) %>%
    dplyr::mutate(
      Scenario = factor(
          scenario,
          c("July-19 [Low R]", "July-19", "July-19 [High R]",
            "JulyGradual [Low R]", "JulyGradual", "JulyGradual [High R]")
        ),
      "Cross immunity" = strain_cross_immunity,
      "Waning rate" = waning_rate,
      "Vaccine effectiveness" = strain_vaccine_efficacy
    ) %>%
    dplyr::select(-dplyr::contains("_"), -scenario)

  cols <- spimalot::spim_scenario_cols(unique(summary$Scenario))
  cols <- cols[levels(summary$Scenario)]

  message(sprintf("Calculating shaps (%s)", i))
  shaps <- spimalot::spim_simulation_shaps(summary)

  message(sprintf("Ploting shaps (%s)", i))

  ## combined
  jpeg(sprintf("figs/varimp_%s.jpg", i), width = 8, height = 8,
      units = "in", res = 300)
  plot(plot_shap_varimps(shaps, states, ncol = 1, nrow = 3))
  dev.off()

  ## separate
  jpeg(sprintf("figs/varimp_deaths_%s.jpg", i), width = 8, height = 8,
      units = "in", res = 300)
  plot(plot_shap_varimp(shaps, states[[1]]))
  dev.off()

  jpeg(sprintf("figs/varimp_adm_%s.jpg", i), width = 8, height = 8,
      units = "in", res = 300)
  plot(plot_shap_varimp(shaps, states[[2]]))
  dev.off()

  jpeg(sprintf("figs/varimp_inf_%s.jpg", i), width = 8, height = 8,
      units = "in", res = 300)
  plot(plot_shap_varimp(shaps, states[[3]]))
  dev.off()

  if (i == "cen") {

    shaps$state <- factor(shaps$state, levels = c("Total infections", "Total hospital admissions", "Total additional COVID-19 deaths", ""))
    jpeg("figs/shaps.jpg", width = 15, height = 7,
        units = "in", res = 300)
    suppressMessages(plot_shap_grid(
      shaps,
      c("pessimistic", rep("JulyGradual [Low R]", 3)),
      TRUE, "none"
    )) %>% plot()
    dev.off()

    jpeg("figs/boxplots.jpg", width = 15, height = 7,
        units = "in", res = 300)
    p <- summary %>%
      dplyr::mutate(
        "Waning rate" = factor(`Waning rate`,
                              levels = c("optimistic", "central", "pessimistic")),
        "Cross immunity" = factor(`Cross immunity`,
                                  levels = c("optimistic", "central",
                                            "pessimistic")),
        "Vaccine effectiveness" = factor(
          dplyr::recode(`Vaccine effectiveness`,
                        "optimistic" = "VE: optimistic",
                        "central" = "VE: central",
                        "pessimistic" = "VE: pessimistic"),
          levels = c("VE: optimistic", "VE: central","VE: pessimistic"))
      ) %>%
      ggplot_boxplot(cols) +
        facet_grid(vars(`Waning rate`), vars(`Cross immunity`),
                   labeller = label_both)
    plot(p)
    dev.off()
  }
}
