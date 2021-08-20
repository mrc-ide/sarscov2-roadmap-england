plot_shap_grid <- function(shaps, x_lab_pos, add_text = TRUE,
                           mean = c("lvl", "global", "none")) {

  mean <- match.arg(mean)

  shaps <-
    shaps %>%
    arrange(Var, desc(state))

  ord <- aggregate(value ~ Var, data = shaps, function(x) mean(abs(x)))
  shaps <- shaps[order(factor(shaps$Var, levels = as.character(ord[order(ord$value, decreasing = TRUE), 1]))), ]

  tmp_shaps <- rbind(
    shaps,
    data.frame(value = 0,
               Var = rep(c("Vaccine effectiveness", "Cross immunity",
                           "Waning rate", "Scenario"), each = 3),
               lvl = rep(c("pessimistic", "central", "optimistic"), 4),
               state = "")
  )

  tmp_shaps$Var <- factor(tmp_shaps$Var, levels = unique(tmp_shaps$Var))
  tmp_shaps$lvl <- factor(tmp_shaps$lvl, levels = unique(tmp_shaps$lvl))

  ints <- interaction(tmp_shaps$state, tmp_shaps$Var)
  ints <- factor(ints, levels = unique(ints))
  levels(ints) <- levels(ints)[c(1:3, 16, 4:6, 13, 7:9, 14, 10:12, 15)]
  tmp_shaps$ints <- ints
  tmp_shaps$color <- if_else(grepl("^\\.", tmp_shaps$ints), "white", "black")
                       tmp_shaps$height <- if_else(grepl("^\\.", tmp_shaps$ints), 1, 2)
                                             tmp_shaps <- tmp_shaps[tmp_shaps$ints != ".Waning rate", ]

  p <- ggplot() +
    geom_tile(aes(y = ints, x = lvl, color = color),
              data = tmp_shaps,
              fill = "white", lwd = 1.1) +
    scale_color_manual(values = c(white = "#FFFFFF", black = "#000000"),
                       guide = "none")

  shaps <- rbind(
    shaps,
    data.frame(value = 0,
               Var = rep(c("Vaccine effectiveness", "Cross immunity",
                           "Waning rate"), each = 3),
               lvl = rep(c("pessimistic", "central", "optimistic"), 3),
               state = "")
  )

  shaps$Var <- factor(shaps$Var, levels = unique(shaps$Var))
  shaps$lvl <- factor(shaps$lvl, levels = unique(shaps$lvl))

  ints <- interaction(shaps$state, shaps$Var)
  ints <- factor(ints, levels = unique(ints))
  shaps$ints <- ints
  sec <- vapply(strsplit(as.character(unique(ints)), ".", TRUE), "[[", character(1), 1)

  for (i in seq_along(unique(sec))) {
    if (i == 1) {
      guide_i <- guide_legend(order = i, override.aes = list(lwd = 0))
    } else {
      guide_i <- "none"
    }
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_tile(aes(y = ints, x = lvl, fill = value), lwd = 1.1,
                color = "black",
                data = shaps %>% dplyr::filter(state == sec[[i]])) +
      scale_fill_gradient2(
        low = "blue", high = "red", name = "", guide = guide_i,
        breaks = function(x) c(min(x), 0, max(x)),
        labels = c("min", 0, "max")
      )

    if (add_text) {
      p <- p + geom_text(aes(y = ints, x = lvl,
                             label = scales::comma(round(value))),
                         data = shaps %>% dplyr::filter(state == sec[[i]]))
    }
  }

  pri <- vapply(strsplit(as.character(unique(ints)), ".", TRUE),
                "[[", character(1), 2)

  ord$fontface <- "plain"
  if (mean == "global") {
    ord$value <- prettyNum(format(round(ord$value, 2), nsmall = 2), ",")
    pri <- paste0(pri, "\n", ord[match(pri, ord$Var), 2])
  } else if (mean == "lvl") {
    ord <- aggregate(value ~ Var + state, data = shaps,
                     function(x) mean(abs(x))) %>%
      arrange(Var, desc(state))
    ord <- ord %>%
      group_by(state) %>%
      mutate(fontface = if_else(value == max(value), "bold", "plain")) %>%
      ungroup()
    ord$value <- prettyNum(format(round(ord$value, 2), nsmall = 2), ",")
    pri <- paste0(pri, "\n", arrange(ord, Var, desc(state))$value)
  }

  nums <- as.numeric(factor(rev(unique(shaps$state))))
  ns <- length(unique(shaps$state))
  nl <- length(unique(shaps$Var))
  sec <- rep(c(sec[1:3], ""), 4)

  p <- p +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10)) +
    labs(x = "", y = "") +
    scale_y_discrete(labels = sec)

  if (mean == "global") {
    pos <- (ns / 2 + ns * (seq_len(nl) - 1))
    if (ns %% 2 == 0) pos <- pos + 0.5
  } else if (mean == "lvl") {
    pos <- seq_len(ns * nl)
    x <- rep(x, each = ns)
  } else if (mean == "none") {
    pos <- (ns / 2 + ns * (seq_len(nl) - 1))
    if (ns %% 2 == 0) pos <- pos + 0.5
  }

  p +
    geom_label(aes(x = x, y = y, label = label, fontface = font),
               data.frame(x = x_lab_pos[x_lab_pos == "pessimistic"],
                          y = 2, font = ord$fontface,
                          label = unique(pri)[1]), hjust = "left",
               family = "serif") +
    geom_label(aes(x = x, y = y, label = label, fontface = font),
               data.frame(x = x_lab_pos[x_lab_pos != "pessimistic"],
                          y = seq.int(6, 14, by = 4), font = ord$fontface[1],
                          label = unique(pri)[-1]), hjust = "right",
               family = "serif") +
    theme(text = element_text(size = 10,  family = "serif"))

}


plot_shap_varimp <- function(shaps, State = NULL) {

  if (!is.null(State)) {
    shaps <- shaps %>% dplyr::filter(state == State)
    xlab <- sprintf("Abs. expected difference for %s", tolower(State))
    xlab <- gsub("covid", "COVID", xlab)
  } else {
    xlab <- "Scaled absolute expected difference"
    shaps <- shaps %>%
      group_by(state) %>%
      mutate(value = abs(value) / max(abs(value)))
  }

  aggregate(value ~ Var, data = shaps, function(x) mean(abs(x))) %>%
    arrange(value) %>%
    mutate(Var = factor(Var, levels = Var)) %>%
    ggplot(aes(x = value, y = Var, fill = Var)) +
    geom_bar(stat = "identity") +
    geom_label(aes(label = scales::comma(value)), fill = "white", hjust = "inward") +
    labs(x = xlab) +
    scale_fill_manual(values = scales::seq_gradient_pal("#e7b2da", "#54029b")(seq(0, 1, length.out = 6))) +
    scale_x_continuous(label = scales::comma_format()) +
    theme_classic() +
    theme(legend.position = "n", axis.text = element_text(angle = 10, hjust = 1),
          text = element_text(size = 12,  family = "serif"),
          axis.title.y = element_blank())
}


plot_shap_varimps <- function(shaps, states, ...) {
  patchwork::wrap_plots(lapply(states, function(x) plot_shap_varimp(shaps, x)),
                        ...)
}


ggplot_boxplot <- function(summary_data, cols) {

  summary_data %>% dplyr::filter(region == "england",
                                 state == "Total additional COVID-19 deaths") %>%
    ggplot(aes(x = `Vaccine effectiveness`, y = `50%`, ymin = `2.5%`, ymax = `97.5%`,
               group = Scenario, fill = Scenario)) +
    geom_crossbar(
      position = position_dodge(width = 0.7, preserve = "single"),
      fatten = 2, width = 0.7,
      lwd = 0.7)  +
    scale_fill_manual(values = cols) +
    ylab("Total additional COVID-19 deaths") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          text = element_text(size = 12,  family = "serif")) +
    scale_y_log10(minor_breaks = 12, n.breaks = 6, labels = scales::comma)
}
