fix_data_aggregation <- function(agg_data) {
  regions <- sircovid::regions("england")
  f <- function(what) {
    rowSums(sapply(agg_data[regions], function(x) x$full[[what]]), na.rm = TRUE)
  }
  what <- c("icu", "general", "hosp",  "admitted","diagnoses", "all_admission")
  agg_data$england$full[, what] <- sapply(what, f)
  agg_data
}


clean_summary <- function(summary, fcts, eng_only, highR_scen, lowR_scen) {

  states <- c("diagnoses_admitted_inc" = "Daily hospital admissions",
              "deaths_inc" = "Daily deaths",
              "deaths_hosp_inc" = "Daily hospital deaths",
              "hosp" = "Hospital bed occupancy",
              "infections_inc" = "Daily infections",
              "infections" = "Total additional COVID-19 infections",
              "deaths" = "Total additional COVID-19 deaths")

  if (fcts) {
    summary <- summary %>%
      dplyr::mutate(analysis = factor(analysis, unique(analysis)),
                    scenario = factor(scenario, unique(scenario)))
  }

  if (eng_only) {
    summary <- summary %>%
      dplyr::filter(region == "england")
  }

  if (highR_scen && any(grepl("low", summary$adherence_to_baseline_npis))) {
    summary <- summary %>% mutate(
      scenario = if_else(grepl("low", .[["adherence_to_baseline_npis"]]),
                         paste(as.character(.[["scenario"]]), "[High R]"),
                         as.character(.[["scenario"]])),
      scenario = factor(scenario, levels = unique(scenario)),
      analysis = if_else(analysis == "High R after full lift",
                         "Central", as.character(analysis)),
      analysis = gsub(" - high R after full lift", "", analysis),
      analysis = factor(analysis, levels = unique(analysis)),
    )
  }


  if (lowR_scen && any(grepl("high", summary$adherence_to_baseline_npis))) {
    summary <- summary %>% mutate(
      scenario = if_else(grepl("high", .[["adherence_to_baseline_npis"]]),
                         paste(as.character(.[["scenario"]]), "[Low R]"),
                         as.character(.[["scenario"]])),
      scenario = factor(scenario, levels = unique(scenario)),
      analysis = if_else(analysis == "Low R after full lift",
                         "Central", as.character(analysis)),
      analysis = gsub(" - low R after full lift", "", analysis),
      analysis = factor(analysis, levels = unique(analysis)),
    )
  }

  summary %>%
    dplyr::mutate(state = as.character(state),
                  state = if_else(state %in% names(states),
                                  states[state], state)) %>%
    tidyr::pivot_wider(names_from = quantile)
}


mean_ci <- function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975)))
