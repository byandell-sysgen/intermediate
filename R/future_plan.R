future_plan <- function(cores = 1) {
  # Set up future plan for package `furrr`
  if(cores == 1 | !future::supportsMulticore()) {
    future::plan(future::sequential)
  } else {
    if(cores < 0) cores <- 0
    if(cores == 0) cores <- future::availableCores() - 1
    future::plan(future::multicore, workers = cores)
  }
  invisible()
}
