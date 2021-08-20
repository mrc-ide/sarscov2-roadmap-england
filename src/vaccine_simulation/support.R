## FIXME - ROUND TO NEAREST 100
output_table <- function(summary, levels = NULL) {

  summary %>%
    dplyr::filter(region == "england") %>%
    dplyr::mutate(
      value = dplyr::if_else(state == "peak_hosp_date",
                             as.character(sircovid::sircovid_date_as_date(value)),
                             as.character(value))) %>%
    dplyr::group_by(analysis, scenario, state) %>%
    dplyr::select(analysis, scenario, state, quantile, value) %>%
    tidyr::pivot_wider(names_from = quantile) %>%
    dplyr::mutate(value = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
    dplyr::select(analysis, scenario, state, value) %>%
    tidyr::pivot_wider(names_from = state) %>%
    dplyr::select(deaths, diagnoses_admitted, infections, peak_hosp,
                  peak_hosp_date) %>%
    dplyr::arrange(analysis)
}


export_rts <- function(sims, type, names, regions, ratio = FALSE,
                       range = NULL, times = c(1, 8)) {

  out <- lapply(sims, function(x) {
    out <- apply(x[[type]][, , times, ], c(2, 3, 4), mean)
    if (ratio) {
      colMeans(out[, , 1] / out[, , 2])
    } else {
      cbind(out[, , 1], out[, , 2])
    }
  })

  out <- data.frame(do.call(rbind, out))
  rownames(out) <- NULL

  if (ratio) {
    colnames(out) <- c("T1", "T8")
    out$analysis <- names
    out <- out[, c(3, 1:2)]
    out <- out[out$T1 <= range[[2]] & out$T1 >= range[[1]] &
      out$T8 <= range[[2]] & out$T8 >= range[[1]], ]
  } else {
    colnames(out) <- c("T1V1", "T8V1", "T1V2", "T8V2")
    out$analysis <- rep(names, each = length(regions))
    out$region <- regions
    out <- out[, c(5:6, 1:4)]
  }

  out
}
