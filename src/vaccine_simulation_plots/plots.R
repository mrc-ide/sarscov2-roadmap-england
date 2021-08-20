ggplot_Rt <- function(df, dat) {

  states <- c(
    "Rt_general_both" = "R(t) excluding immunity",
    "eff_Rt_general_both" = "Effective R(t)"
  )

  df <- df %>% dplyr::filter(
                 analysis == "VOC Central VE",
                 grepl("general_both", state),
                 scenario == "JulyGradual"
               ) %>%
    dplyr::mutate(
      state = if_else(state %in% names(states), states[state], state),
                R = state,
      lb = `2.5%`, c = `50%`, ub = `97.5%`, alpha = 0.4) %>%
    dplyr::select(date, R, lb, c, ub, alpha)

  # Get relevant betas to current date and filter out school holidays
  betas <- data.frame(
    dates = as.Date(tail(dat$samples[[1]]$info$beta_date, 12)),
    label = c(
      "End of 2nd\nLockdown",
      "School Holidays",
      "Holiday \nRestrictions",
      "Start of 3rd\nLockdown",
      "Roadmap\nStep 1",
      "School Holidays",
      "Roadmap\nStep 2",
      "Roadmap\nStep 3",
      "Step 4\nDelayed",
      "Euro 2020\nQtr Final",
      "Euro 2020\nFinal",
      "Roadmap\nStep 4")
  ) %>%
    dplyr::filter(!stringr::str_detect(label, "School")) %>%
    dplyr::arrange(dates) %>%
    dplyr::mutate(label_y = c(4.5, 5.5, 4.5, 5, 4.5, 5, 3.7, 5.6, 4.8, 3.8))

  rts <- rbind(
    t(apply(dat$rt$england$eff_Rt_general[-1, ], 1, quantile,
            probs = c(0.025, 0.5, 0.975))) %>% data.frame() %>%
    dplyr::mutate(R = "Effective R(t)"),
    t(apply(dat$rt$england$Rt_general[-1, ], 1, quantile,
            probs = c(0.025, 0.5, 0.975))) %>% data.frame() %>%
    dplyr::mutate(R = "R(t) excluding immunity")
  ) %>%
    `colnames<-`(c("lb", "c", "ub", "R")) %>%
    data.frame(date = sircovid::sircovid_date_as_date(dat$rt$england$date[-1, 1])) %>%
    tibble() %>%
    dplyr::filter(date <= "2021-07-19") %>%
    dplyr::mutate(
      alpha = 0.1
    )

  plot_df <- rbind(rts, df) %>% dplyr::filter(date >= as.Date("2020-12-01"))


  ggplot(plot_df) +
    geom_line(aes(x = date, y = c, group = R, color = R)) +
    geom_ribbon(aes(x = date, ymin = lb, ymax = ub, group = R, fill = R), alpha = 0.4) +
    theme_minimal() +
    scale_x_date(date_breaks = "months", date_labels = "%d-%b") +
    annotate(geom = "rect",
             xmin = as.Date(c("2020-12-18", "2021-04-01", "2021-07-24", "2021-10-23", "2021-12-18")),
             xmax = as.Date(c("2021-01-05", "2021-04-19", "2021-08-31", "2021-10-31", "2021-12-31")),
             ymin = -Inf, ymax = 6, fill = "grey", alpha = 0.4) +
    geom_hline(yintercept = 1, lty = 2, col = "black") +
    geom_segment(aes(x = x, xend = x, y = 0, yend = 6),
                 data.frame(x = as.Date(betas[, 1])), lty = 3, col = "red4") +
    geom_label(aes(x = dates, label = label, y = label_y), data = betas,
               hjust = 0.5, size = 3,
               vjust = 0.5, family = "serif",
               label.padding = unit(0.15, "lines")) +
    ylab("R(t)") +
    scale_x_date("Date (2020-21)", limits = as.Date(c("2020-11-14", "2021-12-31")),
                 date_breaks = "months",
                 date_labels = "%d-%b", expand = c(0, 0.1)) +
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(),
          text = element_text(family = "serif", size=10),
          legend.title = element_blank(),
          legend.position = "bottom",
          axis.ticks = element_line(),
          plot.margin = unit(rep(0, 4), units = "cm")) +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(breaks = 0:6, labels = 0:6, position = "right", expand = c(0, 0)) +
    ggthemes::scale_colour_colorblind() +
    ggthemes::scale_fill_colorblind() +

  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = 6, ymax = 6.5),
    data.frame(xmin = as.Date(c("2020-12-01", "2021-03-08", "2021-07-19")),
               xmax = as.Date(c("2021-03-08", "2021-07-19", "2021-12-31"))),
    color = "black", fill = rev(khroma::color("vibrant")(3)[1:3])
  ) + geom_text(
        aes(x = xmin + (xmax - xmin) / 2, y = 6.25),
        label = c("Single strain model", "Multi-strain model", "Forward projection"),
        data.frame(xmin = as.Date(c("2020-12-01", "2021-03-08", "2021-07-19")),
                   xmax = as.Date(c("2021-03-08", "2021-07-19", "2021-12-31"))),
        color = "white", family = "serif", size = 5
      )
}

ggplot_trajectories <- function(df, State, fits_df = NULL,
                                filter = "Central",
                                ylab = State,
                                upper_ylim = NULL,
                                real_data = NULL,
                                policy_dates = NULL,
                                title = NULL,
                                legend.position = "n",
                                y_axis_labels = TRUE,
                                cols = NULL, xlim = NULL, log = FALSE,
                                label_x = "2021-11-19", label_nudge = 100) {

  min_date <- as.Date("2021-06-01")
  max_date <- as.Date("2021-12-31")
  df <- df %>% dplyr::filter(analysis == filter,
                             state == State,
                             date <= max_date, date >= min_date)

  ylab <- State
  if (log) {
    ylab <- sprintf("%s (log10)", ylab)
  }

  p <-
    ggplot(df, aes(x = date)) +
    geom_vline(xintercept = policy_dates,
               color = "lightgray", linetype = "dashed") +
    theme_minimal() +
    geom_ribbon(aes(y = `50%`, ymin = `2.5%`, ymax = `97.5%`, group = scenario,
                    fill = scenario), lwd = 0.1, alpha = 0.2) +
    geom_line(aes(y = `50%`, group = scenario, color = scenario), lwd = 1,
              key_glyph = "rect") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.minor.y = element_line(), panel.grid.major.y = element_line(),
          axis.line.y = element_line(), axis.line.x.bottom = element_line(),
          axis.ticks = element_line(),
          legend.position = legend.position,
          legend.key = element_rect(size = 1, colour = "white"),
          legend.key.size = unit(4, "mm"),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.title = element_blank(),
          text = element_text(family = "serif", size = 10),
          axis.title = element_text(family = "serif", size = 9.5)) +
    labs(y = ylab, color = "scenario", fill = "scenario", x = "Date (2021-22)") +
    scale_x_date(limits = c(min_date, max_date), expand = c(0, 0),
                 date_labels = "%d-%b", breaks = "1 month") +
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
    scale_y_continuous(limits = c(0, upper_ylim),
                       breaks = seq.int(0, upper_ylim, length.out = 3),
                       expand = c(0, 0), labels = scales::comma_format())

  if (!is.null(fits_df)) {
    fits_df <- fits_df %>%
      dplyr::filter(state == State, date >= min_date)
    p <- p +
      geom_ribbon(data = fits_df, aes(x = date, y = `50%`, ymin = `2.5%`,
                                      ymax = `97.5%`), lwd = 0.1, alpha = 0.4, fill = "grey70",
                  inherit.aes = FALSE) +
      geom_line(data = fits_df, aes(x = date, y = `50%`),
                lwd = 1, color = "grey70", key_glyph = "rect",
                inherit.aes = FALSE)
  }

  if (!is.null(real_data)) {
    fitted_data <- real_data %>%
      dplyr::filter(date < as.Date("2021-07-19")) %>%
      data.frame()
    new_data <- real_data %>%
      dplyr::filter(date >= as.Date("2021-07-19")) %>%
      data.frame()

    p <- p +
      geom_point(aes(x = x, y = y),
                 data = data.frame(x = fitted_data[, 2], y = fitted_data[, 1]),
                 pch = 23, fill = grey(0.9), color = grey(0.1), cex = 1) +
      geom_point(aes(x = x, y = y),
                 data = data.frame(x = new_data[, 2], y = new_data[, 1]),
                 pch = 23, fill = "orange", color = grey(0.1), cex = 1) +
      geom_hline(
        yintercept = max(
          real_data[real_data$date >= as.Date("2020-12-01"), 1],
          na.rm = TRUE), lty = 2
      ) +
      geom_label(aes(x = x, y = y, label = sprintf("Winter peak: %s", y)), data = data.frame(
                                                                             x = as.Date(label_x),
                                                                             y = max(real_data[real_data$date >= as.Date("2020-12-01"), 1],
                                                                                     na.rm = TRUE) - label_nudge), size = 3, family = "serif")
  }


  if (!y_axis_labels) {
    p <- p + theme(axis.text.y = element_blank())
  }

  if (!is.null(title)) {
    p <- p + labs(title = title)
  }

  if (log) {
    p <- p + scale_y_log10()
  }

  p
}


ggplot_vacc <- function(df, restart_df, parent_df, agg_data, p_elig, pop) {

  cols <- hcl.colors(6, "Bluyl")
  min_date <- as.Date("2020-12-01")
  max_date <- as.Date("2021-12-31")

  ## Before 8 March is Alpha only
  parent_df <- parent_df %>%
    dplyr::filter(date < "2021-03-08") %>%
    dplyr::select(date, state, mean)

  alpha_prop <- agg_data %>%
    dplyr::filter(date >= "2021-07-19") %>%
    dplyr::mutate(alpha_prop = strain_non_variant / strain_tot) %>%
    dplyr::select(date, alpha_prop) %>%
    rbind(data.frame(date = seq.int(max(agg_data$date) + 1, max(df$date), 1),
                     alpha_prop = 0))

  df <- df %>%
    dplyr::filter(strain == "strain_2", analysis == "VOC Central VE",
                  scenario == "JulyGradual") %>%
    dplyr::left_join(alpha_prop) %>%
    dplyr::mutate(
      mean = dplyr::if_else(strain == "strain_1", mean * alpha_prop,
                            mean * (1 - alpha_prop))
    ) %>%
    dplyr::group_by(date, state) %>%
    dplyr::summarise(mean = sum(mean, na.rm = TRUE)) %>%
    dplyr::ungroup()

  ## Weight everything in between
  alpha_prop <- agg_data %>%
    dplyr::filter(date > "2021-03-08", date < "2021-07-19") %>%
    dplyr::mutate(alpha_prop = strain_non_variant / strain_tot) %>%
    dplyr::select(date, alpha_prop)

  restart_df <- restart_df %>%
    dplyr::filter(assumptions == "central", date < "2021-07-19") %>%
    dplyr::left_join(alpha_prop) %>%
    dplyr::mutate(
      mean = dplyr::if_else(strain == "strain_1", mean * alpha_prop,
                            mean * (1 - alpha_prop))
    ) %>%
    dplyr::group_by(date, state) %>%
    dplyr::summarise(mean = sum(mean, na.rm = TRUE)) %>%
    dplyr::ungroup()

  plot_df <- rbind(parent_df, restart_df, df)

  states <- c(ever_vaccinated = "Unprotected despite vaccination",
              protected_against_death = "Protected vs death",
              protected_against_severe_disease = "Protected vs severe disease",
              protected_against_infection = "Protected vs infection after vaccination",
              ever_infected = "Protected vs infection after infection and vaccination",
              ever_infected_unvaccinated = "Protected vs infection after infection")

  plot_df <- plot_df %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(mean = if_else(state %in% c("ever_infected",
                                              "ever_infected_unvaccinated"),
                                 mean,
                                 mean + mean[state == "ever_infected_unvaccinated"]),
                           prop = mean / pop,
                  state = factor(state, levels = names(states),
                                 labels = states)) %>%
    dplyr::arrange(date, state) %>%
    dplyr::select(state, date, prop)

  ggplot(data = plot_df, aes(x = date)) +
    geom_ribbon(aes(x = date, ymin = 0, ymax = prop, fill = state),
                colour = "transparent")  +
    geom_vline(aes(xintercept = max(restart_df$date)), col = "grey30", lty = 2) +
    theme_minimal() +
    scale_x_date("Date (2020-21)", limits = as.Date(c("2020-11-14", "2021-12-31")),
                 date_breaks = "months",
                 date_labels = "%d-%b", expand = c(0, 0.1)) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(breaks = seq.int(0, 1, length.out = 6),
                       labels = seq.int(0, 100, length.out = 6),
                       expand = c(0.01, 0), limits = c(0, 1),
                       position = "right") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.y = element_line(),
          text = element_text(family = "serif", size = 10),
          axis.line.x.bottom = element_line(),
          axis.ticks = element_line(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(family = "serif", size = 8.5)) +
    labs(y = "Proportion of population (%)") +
    guides(fill = guide_legend(nrow = 1, ncol = 6))
}

plot_state_by_age <- function(df, fits_df) {

  groups <- c(age_0 = "Age <30", age_30 = "Age 30-49", age_50 = "Age 50-74",
              age_75 = "Age 75+")
  n_age <- length(unique(groups))
  vaccine_statuses <- c(unvaccinated = "No vaccine protection",
                        partial_protection = "Single dose protection",
                        full_protection = "Two dose protection")

  df <- df %>%
    dplyr::filter(vaccine_status != "booster",
                  region == "england",
                  scenario == "JulyGradual") %>%
    dplyr::select(analysis, date, state, group, vaccine_status, mean)
  fits_df <- fits_df %>%
    dplyr::filter(region == "england", date <= as.Date("2021-07-19")) %>%
    tidyr::expand_grid(analysis = unique(df$analysis)) %>%
    dplyr::select(analysis, date, state, group, vaccine_status, mean)

  plot_df <- dplyr::bind_rows(df, fits_df) %>%
    dplyr::filter(state %in% c("Daily infections",
                               "Daily hospital admissions",
                               "Daily deaths")) %>%
    dplyr::mutate(group = factor(group, levels = names(groups), labels = groups,
                                 ordered = TRUE),
                  vaccine_status = factor(vaccine_status,
                                          levels = names(vaccine_statuses),
                                          labels = vaccine_statuses,
                                          ordered = TRUE),
                  age_vacc = paste(vaccine_status, group),
                  analysis = factor(analysis, levels = c("High VE", "Central VE", "Low VE")))

  plot_df %>%
    ggplot(aes(x = date, y = mean, fill = age_vacc)) +
    geom_area() +
    theme_bw() +
    theme(axis.line = element_line(),
          legend.position = "right",
          legend.text = element_text(size = 10, family = "serif"),
          text = element_text(size = 10, family = "serif")) +
    labs(x = "Date (2021)", y = "") +
    facet_grid(vars(state), vars(analysis), scales = "free",
               labeller = label_wrap_gen(20)) +
    scale_fill_manual(values = c(hcl.colors(n_age, "Purp", rev = TRUE),
                                 hcl.colors(n_age, "Burg", rev = TRUE),
                                 hcl.colors(n_age, "Teal", rev = TRUE))) +
    scale_x_date(limits = as.Date(c("2021-03-08", "2021-12-31")),
                 breaks = "2 months", date_labels = "%d-%b", expand = c(0, 0))
}


ggplot_full_trajectory <- function(parent, restart, simulation, data, State, cols,
                                   data_state, ylim, policy_dates) {

  fitted_data <- data %>%
    dplyr::filter(date < as.Date("2021-07-19")) %>%
    data.frame()
  new_data <- data %>%
    dplyr::filter(date >= as.Date("2021-07-19")) %>%
    data.frame()

  ggplot() +
    geom_line(aes(x = date, y = `50%`),
              data = restart %>%
                dplyr::filter(state == State, assumptions == "central"),
              color = "blue", lwd = 1
              ) +
    geom_ribbon(data = restart %>%
                  dplyr::filter(state == State, assumptions == "central"),
                aes(x = date, y = `50%`, ymin = `2.5%`,
                    ymax = `97.5%`), lwd = 0.1, alpha = 0.4, fill = "blue",
                inherit.aes = FALSE) +
    geom_line(aes(x = date, y = `50%`),
              data = parent %>%
                dplyr::filter(state == State),
              color = "blue", lwd = 1
              ) +
    geom_ribbon(data = parent %>%
                  dplyr::filter(state == State),
                aes(x = date, y = `50%`, ymin = `2.5%`,
                    ymax = `97.5%`), lwd = 0.1, alpha = 0.4, fill = "blue",
                inherit.aes = FALSE) +
    geom_line(aes(x = date, y = `50%`, group = scenario, color = scenario),
              data = simulation %>%
                dplyr::filter(state == State,
                              analysis == "VOC Central VE",
                              grepl("JulyGradual", scenario)), lwd = 1, key_glyph = "rect") +
    geom_ribbon(data = simulation %>%
                  dplyr::filter(state == State,
                                analysis == "VOC Central VE",
                                grepl("JulyGradual", scenario)),
                aes(x = date, y = `50%`, ymin = `2.5%`, , group = scenario, fill = scenario,
                    ymax = `97.5%`), lwd = 0.1, alpha = 0.2, inherit.aes = FALSE) +
    geom_point(aes_string(x = "date", y = data_state),
               data = fitted_data, pch = 23, fill = grey(0.9), color = grey(0.1), cex = 1) +
    geom_point(aes_string(x = "date", y = data_state),
               data = new_data, pch = 23, fill = "orange", color = grey(0.1), cex = 1) +
    scale_x_date("Date (2021-22)", limits = as.Date(c("2021-01-01", "2021-12-31")),
                 date_labels = "%d-%b", breaks = "month") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.minor.y = element_line(), panel.grid.major.y = element_line(),
          axis.line.y = element_line(), axis.line.x.bottom = element_line(),
          axis.ticks = element_line(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = c(0.7, 0.9),
          legend.key = element_rect(size = 2, colour = "white"),
          legend.key.size = unit(4, "mm"),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.title = element_blank(),
          text = element_text(family = "serif", size = 10)) +
    labs(y = State, color = "Scenario", fill = "Scenario") +
    scale_color_manual(values = cols[grepl("JulyGradual", names(cols))]) +
    scale_fill_manual(values = cols[grepl("JulyGradual", names(cols))]) +
    scale_y_continuous(limits = c(0, ylim),
                       breaks = seq.int(0, ylim, length.out = 3),
                       expand = c(0, 0), labels = scales::comma_format()) +
    geom_vline(xintercept = policy_dates, color = "lightgray",
               linetype = "dashed")

}

vaccine_figure_1 <- function(parent, restart, simulation, data, cols,
                             policy_dates, combined_data) {

  fig_1a <- ggplot_full_trajectory(parent, restart, simulation, data,
                                   "Daily hospital admissions", cols,
                                   "all_admission", 4400, policy_dates)
  fig_1b <- ggplot_trajectories(
    df = simulation %>% dplyr::filter(grepl("JulyGradual", scenario)),
    State = "Daily hospital admissions",
    filter = "VOC Central VE",
    upper_ylim = 25000,
    policy_dates = policy_dates, cols = cols[grepl("JulyGradual", names(cols))],
    real_data = data[, c("all_admission", "date")], log = TRUE
  ) +
    scale_x_date("Date (2021-22)",
                 limits = as.Date(c("2021-06-01", "2021-12-31")),
                 date_labels = "%d-%b", date_breaks = "months") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  fig_1c <- spimalot:::spim_plot_voc_proportion(combined_data, "2021-03-08",
                                                "london") +
    scale_x_date(date_labels = "%d-%b", breaks = "1 month") +
    theme(text = element_text(size = 10, family = "serif"),
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = "Date (2021)", title = "London") +
    geom_vline(xintercept = policy_dates, color = "lightgray",
               linetype = "dashed")

  parent <- parent %>% dplyr::select(date, state, `2.5%`, `50%`, `97.5%`)
  parent <- t(apply(
    dat$samples$england$trajectories$state["deaths_hosp_inc", , ], 2,
    quantile, c(0.025, 0.5, 0.975)
  )) %>%
    data.frame(date = sircovid::sircovid_date_as_date(dat$samples$england$trajectories$date),
               state = "Daily hospital deaths") %>%
    `colnames<-`(c("2.5%", "50%", "97.5%", "date", "state")) %>%
    rbind(parent)

  fig_1d <- ggplot_full_trajectory(parent, restart, simulation, data,
                                   "Daily hospital deaths", cols,
                                   "deaths_hosp", 1000, policy_dates) +
    theme(legend.position = "n")

  fig_1e <- ggplot_trajectories(
    df = simulation %>% dplyr::filter(grepl("JulyGradual", scenario)),
    State = "Daily hospital deaths",
    filter = "VOC Central VE",
    upper_ylim = 2000,
    policy_dates = policy_dates, cols = cols[grepl("JulyGradual", names(cols))],
    real_data = data[, c("deaths_hosp", "date")], log = TRUE
  ) +
    scale_x_date("Date (2021-22)",
                 limits = as.Date(c("2021-06-01", "2021-12-31")),
                 date_labels = "%d-%b", date_breaks = "months") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  fig_1f <- spimalot:::spim_plot_seeding_date(combined_data) +
    scale_x_date(date_labels = "%d-%b", breaks = "3 days") +
    theme(text = element_text(size = 10, family = "serif")) +
    labs(x = "Date (2021)") +
    geom_vline(xintercept = policy_dates, color = "lightgray",
               linetype = "dashed") +
    scale_x_date(date_labels = "%d-%b", date_breaks = "4 days")

                               (fig_1a + fig_1b + fig_1c + fig_1d + fig_1e + fig_1f) +
                                 plot_annotation(tag_levels = "A")
}

vaccine_figure_2 <- function(simulation, restart, parent, dat, agg, prop, pop) {
  fig2_a <- ggplot_Rt(simulation$state, dat)
  fig2_b <- ggplot_vacc(simulation$n_protected,
                        dplyr::filter(restart$n_protected,
                                      assumptions == "central"),
                        parent$n_protected, agg,
                        100 - prop, pop = sum(pop))

  fig2_a + fig2_b + plot_layout(heights = c(1.6, 1), nrow = 2, ncol = 1) +
    plot_annotation(tag_levels = "A")
}

vaccine_figure_3 <- function(fits_state, summary_state, policy_dates) {

  fits_cen <- fits_state %>% dplyr::filter(assumptions == "central")
  fits_pes <- fits_state %>% dplyr::filter(assumptions == "pessimistic")
  fits_opt <- fits_state %>% dplyr::filter(assumptions == "optimistic")

  summary_voc <- summary_state %>% dplyr::filter(grepl("Gradual", scenario))
  fig3cols <- cols[c(
    "JuneGradual [Low R]", "JuneGradual", "JuneGradual [High R]",
    "JulyGradual [Low R]", "JulyGradual", "JulyGradual [High R]"
  )]
  y_inf <- 800000
  y_hosp <- 14000
  y_ddeaths <- 1300
  y_cdeaths <- 60000

  fig3_a2 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_opt,
    State = "Daily infections",
    filter = "VOC High VE",
    upper_ylim = y_inf,
    policy_dates = policy_dates, cols = fig3cols,
    legend.position = c(0.6, 0.8),
    title = "Optimistic immunity scenario",
    y_axis_labels = TRUE
  )
  fig3_a3 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_cen,
    State = "Daily infections",
    filter = "VOC Central VE",
    upper_ylim = y_inf,
    policy_dates = policy_dates, cols = fig3cols,
    title = "Central immunity scenario",
    ylab = NULL, y_axis_labels = FALSE
  )
  fig3_a4 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_pes,
    State = "Daily infections",
    filter = "VOC Low VE",
    upper_ylim = y_inf,
    policy_dates = policy_dates, cols = fig3cols,
    title = "Pessimistic immunity scenario",
    ylab = NULL, y_axis_labels = FALSE
  )


  fig3_b2 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_opt,
    State = "Daily hospital admissions",
    filter = "VOC High VE",
    upper_ylim = y_hosp,
    policy_dates = policy_dates, cols = fig3cols,
    y_axis_labels = TRUE,
    real_data = agg_data[, c("all_admission", "date")],
    label_nudge = -4000, label_x = "2021-11-25"
  )
  fig3_b3 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_cen,
    State = "Daily hospital admissions",
    filter = "VOC Central VE",
    upper_ylim = y_hosp,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE,
    real_data = agg_data[, c("all_admission", "date")],
    label_nudge = -4000, label_x = "2021-11-25"
  )
  fig3_b4 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_pes,
    State = "Daily hospital admissions",
    filter = "VOC Low VE",
    upper_ylim = y_hosp,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE,
    real_data = agg_data[, c("all_admission", "date")],
    label_nudge = -4000, label_x = "2021-11-25"
  )


  fig3_c2 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_opt,
    State = "Daily deaths",
    filter = "VOC High VE",
    upper_ylim = y_ddeaths,
    policy_dates = policy_dates, cols = fig3cols,
    y_axis_labels = TRUE,
    real_data = agg_data[, c("deaths", "date")]
  )
  fig3_c3 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_cen,
    State = "Daily deaths",
    filter = "VOC Central VE",
    upper_ylim = y_ddeaths,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE,
    real_data = agg_data[, c("deaths", "date")]
  )
  fig3_c4 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_pes,
    State = "Daily deaths",
    filter = "VOC Low VE",
    upper_ylim = y_ddeaths,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE,
    real_data = agg_data[, c("deaths", "date")]
  )



  fig3_d2 <- ggplot_trajectories(
    df = summary_voc,
    State = "Total additional COVID-19 deaths",
    filter = "VOC High VE",
    upper_ylim = y_cdeaths,
    policy_dates = policy_dates, cols = fig3cols,
    y_axis_labels = TRUE
  )
  fig3_d3 <- ggplot_trajectories(
    df = summary_voc,
    State = "Total additional COVID-19 deaths",
    filter = "VOC Central VE",
    upper_ylim = y_cdeaths,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE
  )
  fig3_d4 <- ggplot_trajectories(
    df = summary_voc,
    State = "Total additional COVID-19 deaths",
    filter = "VOC Low VE",
    upper_ylim = y_cdeaths,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE
  )

  ((fig3_a2 | fig3_a3 | fig3_a4) + plot_layout(tag_level = "new")) /
    ((fig3_b2 | fig3_b3 | fig3_b4) + plot_layout(tag_level = "new")) /
    ((fig3_c2 | fig3_c3 | fig3_c4) + plot_layout(tag_level = "new")) /
    ((fig3_d2 | fig3_d3 | fig3_d4) + plot_layout(tag_level = "new")) +
    plot_layout(guides = "collect") + plot_annotation(tag_levels = c("A", "1")) &
    theme(legend.position = "right", text = element_text(family = "serif", size = 10))
}


vaccine_figure_SI_traj <- function(fits_state, summary_state, policy_dates) {

  fits_cen <- fits_state %>% dplyr::filter(assumptions == "central")
  fits_pes <- fits_state %>% dplyr::filter(assumptions == "pessimistic")
  fits_opt <- fits_state %>% dplyr::filter(assumptions == "optimistic")

  reg_scen <- "(AlphaOnly)|(July-19 \\[Low R\\])|(July-19 \\[High R\\])|(June-19 \\[Low R\\])|(June-19 \\[High R\\])"
  summary_AlphaOnly <- summary_state %>% dplyr::filter(grepl("AlphaOnly", scenario))
  summary_voc <- summary_state %>% dplyr::filter(!grepl(reg_scen, scenario))
  fig3cols <- cols[c("AlphaOnly", "AlphaOnly [Low R]", "AlphaOnly [High R]",
                     "June-21", "JuneGradual", "JuneGradual [Low R]",
                     "JuneGradual [High R]", "July-19", "JulyGradual",
                     "JulyGradual [Low R]", "JulyGradual [High R]")]

  y_inf <- 1600000
  y_hosp <- 28000
  y_ddeaths <- 2500
  y_cdeaths <- 70000

  fig3_a1 <- ggplot_trajectories(
    df = summary_AlphaOnly,
    fits_df = fits_cen,
    State = "Daily infections",
    filter = "VOC Central VE",
    upper_ylim = y_inf,
    policy_dates = policy_dates,
    cols = fig3cols,
    title = "Alpha Only (Central)"
  )
  fig3_a2 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_opt,
    State = "Daily infections",
    filter = "VOC High VE",
    upper_ylim = y_inf,
    policy_dates = policy_dates, cols = fig3cols,
    legend.position = c(0.6, 0.8),
    title = "VOC (Optimistic)",
    ylab = NULL, y_axis_labels = FALSE
  )
  fig3_a3 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_cen,
    State = "Daily infections",
    filter = "VOC Central VE",
    upper_ylim = y_inf,
    policy_dates = policy_dates, cols = fig3cols,
    title = "VOC (Central)",
    ylab = NULL, y_axis_labels = FALSE
  )
  fig3_a4 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_pes,
    State = "Daily infections",
    filter = "VOC Low VE",
    upper_ylim = y_inf,
    policy_dates = policy_dates, cols = fig3cols,
    title = "VOC (Pessimistic)",
    ylab = NULL, y_axis_labels = FALSE
  )

  fig3_b1 <- ggplot_trajectories(
    df = summary_AlphaOnly,
    fits_df = fits_cen,
    State = "Daily hospital admissions",
    filter = "VOC Central VE",
    upper_ylim = y_hosp, cols = fig3cols,
    policy_dates = policy_dates,
    real_data = agg_data[, c("all_admission", "date")],
    label_x = "2021-11-01"
  )
  fig3_b2 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_opt,
    State = "Daily hospital admissions",
    filter = "VOC High VE",
    upper_ylim = y_hosp,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE,
    real_data = agg_data[, c("all_admission", "date")],
    label_x = "2021-11-01"
  )
  fig3_b3 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_cen,
    State = "Daily hospital admissions",
    filter = "VOC Central VE",
    upper_ylim = y_hosp,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE,
    real_data = agg_data[, c("all_admission", "date")],
    label_x = "2021-11-01"
  )
  fig3_b4 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_pes,
    State = "Daily hospital admissions",
    filter = "VOC Low VE",
    upper_ylim = y_hosp,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE,
    real_data = agg_data[, c("all_admission", "date")],
    label_x = "2021-11-01"
  )

  fig3_c1 <- ggplot_trajectories(
    df = summary_AlphaOnly,
    fits_df = fits_cen,
    State = "Daily deaths",
    filter = "VOC Central VE",
    upper_ylim = y_ddeaths, cols = fig3cols,
    policy_dates = policy_dates,
    real_data = agg_data[, c("deaths", "date")],
    label_x = "2021-11-01"
  )
  fig3_c2 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_opt,
    State = "Daily deaths",
    filter = "VOC High VE",
    upper_ylim = y_ddeaths,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE,
    real_data = agg_data[, c("deaths", "date")],
    label_x = "2021-11-01"
  )
  fig3_c3 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_cen,
    State = "Daily deaths",
    filter = "VOC Central VE",
    upper_ylim = y_ddeaths,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE,
    real_data = agg_data[, c("deaths", "date")],
    label_x = "2021-11-01"
  )
  fig3_c4 <- ggplot_trajectories(
    df = summary_voc,
    fits_df = fits_pes,
    State = "Daily deaths",
    filter = "VOC Low VE",
    upper_ylim = y_ddeaths,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE,
    real_data = agg_data[, c("deaths", "date")],
    label_x = "2021-11-01"
  )


  fig3_d1 <- ggplot_trajectories(
    df = summary_AlphaOnly,
    State = "Total additional COVID-19 deaths",
    filter = "VOC Central VE",
    upper_ylim = y_cdeaths, cols = fig3cols,
    policy_dates = policy_dates
  )
  fig3_d2 <- ggplot_trajectories(
    df = summary_voc,
    State = "Total additional COVID-19 deaths",
    filter = "VOC High VE",
    upper_ylim = y_cdeaths,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE
  )
  fig3_d3 <- ggplot_trajectories(
    df = summary_voc,
    State = "Total additional COVID-19 deaths",
    filter = "VOC Central VE",
    upper_ylim = y_cdeaths,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE
  )
  fig3_d4 <- ggplot_trajectories(
    df = summary_voc,
    State = "Total additional COVID-19 deaths",
    filter = "VOC Low VE",
    upper_ylim = y_cdeaths,
    policy_dates = policy_dates, cols = fig3cols,
    ylab = NULL, y_axis_labels = FALSE
  )

  ((fig3_a1 | fig3_a2 | fig3_a3 | fig3_a4) + plot_layout(tag_level = "new")) /
    ((fig3_b1 | fig3_b2 | fig3_b3 | fig3_b4) + plot_layout(tag_level = "new")) /
    ((fig3_c1 | fig3_c2 | fig3_c3 | fig3_c4) + plot_layout(tag_level = "new")) /
    ((fig3_d1 | fig3_d2 | fig3_d3 | fig3_d4) + plot_layout(tag_level = "new")) +
    plot_layout(guides = "collect") + plot_annotation(tag_levels = c("A", "1")) &
    theme(legend.position = "right", text = element_text(family = "serif", size = 10),
          axis.text.x = element_text(angle = 30, hjust = 1))
}


vaccine_figure_SI_state <- function(summary, fits) {
  summary %>%
    dplyr::mutate(analysis = gsub("VOC ", "", analysis)) %>%
    plot_state_by_age(
      fits %>% dplyr::filter(assumptions == "central")
    ) + labs(fill = "Vaccination status / Age group")
}


vaccine_figure_SI_prop <- function(dat) {
  patchwork::wrap_plots(
    lapply(sircovid::regions("england"),
           function(x) spimalot:::spim_plot_voc_proportion(
                         dat, "2021-03-08", x))) &
    scale_x_date(date_labels = "%d-%b", breaks = "1 month") &
    theme(text = element_text(size = 10, family = "serif"),
          axis.text.x = element_text(angle = 30, hjust = 1)) &
    labs(x = "Date (2021)")
}
