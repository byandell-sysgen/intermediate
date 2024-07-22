#'
#' Mediation Scan
#' Scan set of mediators for a target.
#' 
#' @param target A numeric vector with gene/protein expression
#' @param mediator A matrix, each column is one gene/protein's expression
#' @param driver A matrix, haplotype probabilities at QTL we try to mediate
#' @param annotation A data frame with mediators' annotation with columns for facet and index
#' @param covar A matrix with additive covariates
#' @param intcovar A matrix of covariate interacting with driver
#' @param method A method to handle missing cases
#' @param fitFunction function to fit models
#' @param annotation_names names in annotation of columns for facet, index and, optionally, driver (default `c(facet = "chr", index = "pos", driver = NULL)`)
#' @param verbose If TRUE display information about the progress
#' @param cores use multiple cores if value is 0 or >1
#' @param ... additional parameters
#' 
#' @details 
#' For a given QTL haplotype probabilities `driver`` and target `target`,
#' the function sequentially tries to add each column of `mediator` matrix as a covariate
#' and calculates LR statistic. The low LR value indicates `driver` and
#' `target` are conditionally independent given `mediator`,
#' i.e. `mediator` is a mediator of causal relationship from `driver` to `target`.
#'
#' @examples
#' data(Tmem68)
#' 
#' target <- Tmem68$target
#' m <- match("Tmem68", Tmem68$annotation$symbol)
#' 
#' med_scan <- mediation_scan(target = target,
#'                       mediator = Tmem68$mediator,
#'                       driver = Tmem68$driver,
#'                       annotation = Tmem68$annotation,
#'                       covar = Tmem68$covar)
#'                       
#' summary(med_scan)
#' 
#' ggplot_mediation_scan(med_scan)
#' 
#' @export
#' @importFrom purrr transpose
#' @importFrom future availableCores multicore plan sequential
#'             supportsMulticore
#' @importFrom furrr future_map

mediation_scan <- function(target, 
                           mediator, 
                           driver, 
                           annotation,
                           covar=NULL,
                           intcovar=NULL,
                           method=c("double-LR-diff", "ignore", "LR-diff"), 
                           fitFunction = fitDefault,
                           annotation_names = c(facet = "chr", index = "pos", driver = NULL),
                           verbose=TRUE,
                           cores = 0, ...) {
  future_plan(cores)
  # Get common data.
  commons <- common_data(target, mediator, driver, covar, intcovar = intcovar,
                         ...)
  if(is.null(commons))
    return(NULL)
  
  target <- commons$target
  mediator <- commons$mediator
  driver <- commons$driver
  covar <- commons$covar_tar
  intcovar <- commons$intcovar
  common <- commons$common
  rm(commons)
  
  # Fit model without mediation.
  if(length(dim(driver)) == 2)
    driver_tar <- driver
  else
    driver_tar <- driver[,,1]
  loglik0 <- fitFunction(driver_tar, target, covar, intcovar, ...)$LR
  
  # check input
  stopifnot(all(tolower(annotation_names) %in% tolower(names(annotation))))
  method = match.arg(method)
  
  # Match up annotation with mediators
  stopifnot(all(!is.na(m <- match(colnames(mediator), annotation$id))))
  annotation <- annotation[m,]
  rownames(annotation) <- colnames(mediator)
  
  med_pur <- 
    purrr::transpose(
      list(mediator = as.data.frame(mediator),
           annotation = split(annotation, rownames(annotation))))
  LRfn <- 
    switch(
      method,
      ignore            = function(loglik, loglik0) loglik[1],
      "LR-diff"        = function(loglik, loglik0) loglik[2] - loglik[1],
      "double-LR-diff" = function(loglik, loglik0) loglik0 - (loglik[2] - loglik[1]))
  
  mapfn <- function(x, target, covar, driver, loglik0) {
    if(length(dim(driver)) > 2) {
      if(is.null(dcol <- x$annotation$driver_names))
        driver <- driver[,,1]
      else
        driver <- driver[,, dcol]
    }

    loglik <- fitFunction(driver, target, cbind(covar, x$mediator), intcovar, ...)$LR
    if(!is.null(x$mediator)) {
      if(is.matrix(x$mediator)) {
        na <- apply(x$mediator, 1, function(x) any(is.na(x)))
      } else {
        na <- is.na(x$mediator)
      }
      target[na] <- NA
    }
    loglik <- c(loglik,
                fitFunction(driver, target, covar, intcovar, ...)$LR)
    LRfn(loglik, loglik0)
  }
  output <- annotation
  # Compute Likelihood Ratio
  output$LR <- unlist(furrr::future_map(med_pur, mapfn,
    target, covar, driver, loglik0, .progress = TRUE))
  attr(output, "targetFit") <- loglik0
  attr(output, "annotation_names") <- annotation_names
  class(output) <- c("mediation_scan", "data.frame")
  return(output)
}
#' @param facets names of facets to subset
#' @param chrs chromosome names to subset
#' 
#' @export
#' @rdname mediation_scan
#' 
subset.mediation_scan <- function(x, facets=NULL, chrs = NULL, ...) {
  facet_name <- attr(x, "annotation_names")["facet"]
  if(!is.null(facets)) {
    new_x <- x[x[[facet_name]] %in% facets,]
    x <- modify_object(x, new_x)
  }
  if(!is.null(chrs)) {
    new_x <- x[x[["chr"]] %in% chrs,]
    x <- modify_object(x, new_x)
  }
  x
}
#' @export
#' @rdname mediation_scan
#' 
#' @param n maximum number of mediators to show (default 10)
#' @param minimal show only "symbol", annontation_names and "LR" if \code{TRUE}
#' 
#' @importFrom dplyr arrange select
#' @importFrom rlang .data
#' @importFrom utils head
#' 
summary.mediation_scan <- function(object, n = 10, minimal = FALSE, ...) {
  if(minimal) {
    annotation_names <- attr(object, "annotation_names")
    cols <- c("symbol", annotation_names, "LR")
    m <- match(cols, names(object), nomatch = 0)
    object <- object[, m, drop = FALSE]
  }
  utils::head(dplyr::arrange(object, .data$LR), n = n)
}