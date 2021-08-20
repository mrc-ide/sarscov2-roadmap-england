combined_onward_restart <- function(onward, restart_date, 
                                    path, regions = "all") {
  regions <- sircovid::regions(regions)
  
  files <- file.path(path, regions, "restart.rds")
  
  msg <- !file.exists(files)
  if (any(msg)) {
    msg <- sprintf("  - %s", file.path(regions[msg], "restart.rds"))
    stop(sprintf("Missing files at '%s':\n%s",
                 path, paste(msg, collapse = "\n")),
         call. = FALSE)
  }
  
  restart <- Map(readRDS, files, regions)
  names(restart) <- regions
  
  sircovid_restart_date <- sircovid::sircovid_date(restart_date)
  
  onward$date <- restart_date
  onward$step <- sircovid_restart_date * onward$steps_per_day
  onward$state <- lapply(restart, function (x) 
    x$state$state[, , which(x$state$time == sircovid_restart_date)])
  
  onward
}