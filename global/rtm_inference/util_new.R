read_csv <- function(path, ...) {
  read.csv(path, ..., stringsAsFactors = FALSE, check.names = FALSE)
}


write_device <- function(dev, filename, code, ...) {
  dev(filename, ...)
  on.exit(dev.off())
  res <- force(code)
  if (inherits(res, "ggplot")) {
    print(res)
    NULL
  }
}


write_pdf <- function(filename, code, ...) {
  write_device(pdf, filename, code, ...)
}


write_png <- function(filename, code, ...) {
  write_device(png, filename, code, ...)
}


version_check <- function(package, version) {
  if (packageVersion(package) < version) {
    stop(sprintf(
      paste("Please update %s with: drat:::add('ncov-ic');",
            "install.packages('%s')"),
      package, package))
  }
}


write_csv <- function(...) {
  write.csv(..., row.names = FALSE)
}


`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


abind_quiet <- function(...) {
  suppressWarnings(abind::abind(...))
}


switch_levels <- function(x) {
  nms <- names(x[[1]]) %||% seq_along(x[[1]])
  y <- lapply(nms, function(z) lapply(x, "[[", z))
  names(y) <- nms
  y
}
