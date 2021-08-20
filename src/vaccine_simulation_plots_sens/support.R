clean_summary <- function(summary, fcts, eng_only, highR_scen, lowR_scen) {

  states <- c("deaths" = "Total additional COVID-19 deaths",
              "admitted" = "Total hospital admissions",
              "diagnoses" = "Total COVID-19 diagnoses",
              "infections" = "Total infections",
              "diagnoses_admitted" = "Total COVID-19 admissions",
              "peak_hosp" = "Peak hospital bed occupancy",
              "peak_hosp_date" = "Peak hospital bed occupancy date")

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