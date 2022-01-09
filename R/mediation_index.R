# Mediation tests across index
#
#' Test mediation across set of indexed drivers
#'
#' @param target vector or 1-column matrix with target values
#' @param mediator vector or 1-column matrix with mediator values
#' @param driver vector or matrix with driver values
#' @param annotation optional annotation data frame for mediators
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param driver_med driver array for mediators
#' @param annotation_names names in annotation of columns for facet, index and, optionally, driver (default `c(facet = "chr", index = "pos", driver = NULL)`)
#' @param ... additional parameters
#'
#' @importFrom purrr map transpose
#' @importFrom stringr str_replace
#' @importFrom dplyr arrange as_tibble bind_rows desc filter group_by left_join mutate one_of rename ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 aes autoplot element_blank facet_grid facet_wrap 
#' geom_hline geom_point geom_rug geom_vline ggplot guides guide_legend
#' ggtitle scale_color_manual scale_shape_manual theme xlab ylab
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @importFrom RColorBrewer brewer.pal
#' 
#' @examples
#' data(Tmem68)
#'  
#' target <- Tmem68$target
#' 
#' # Mediator is Tmem68.
#' m <- grep("Tmem68", Tmem68$annotation$symbol)
#' mediator <- Tmem68$mediator[,m]
#' annotation <- Tmem68$annotation[m,, drop = FALSE]
#' 
#' # Reduce to SNP for B6
#' driver_SNP <- cbind(B6 = Tmem68$driver[,2], rest = 1 - Tmem68$driver[,2])
#' 
#' med_index <- mediation_index(target = target,
#'                       mediator = mediator,
#'                       driver = Tmem68$driver,
#'                       annotation = annotation,
#'                       covar_tar = Tmem68$covar)
#' summary(med_index)
#' ggplot2::autoplot(med_index)
#' 
#' @export
#'
mediation_index <- function(target, mediator, driver = NULL,
                            annotation = NULL, covar_tar = NULL, covar_med = NULL,
                            driver_med = NULL,
                            annotation_names = c(facet = "chr", index = "pos"),
                            ...) {
  # Mediation test over interval
  
  # Mediation driver is typically a 3-D array with third dimension corresponding to mediator columns.
  # If only one driver, then do reduced version.
  driver <- as.matrix(driver)
  if(is.null(driver_med)) {
    dn <- dimnames(driver)
    dn[[3]] <- "mediator"
    driver_med <- array(driver, c(dim(driver), 1), dn)
  }
  if(3 != length(dim(driver_med)))
    return(NULL)
  nmed <- dim(driver_med)[3]
  
  # Propagate mediator over third dimension of driver_med.
  mediator <- as.matrix(mediator)[, rep(1, nmed), drop = FALSE]
  colnames(mediator) <- dimnames(driver_med)[[3]]

  # Propagate annotation over 
  annotation <- annotation[rep(1, nmed),, drop = FALSE]
  annotation$id <- colnames(mediator)
  annotation$driver_names <- colnames(mediator)
  stopifnot(length(annotation$id) == nmed)
  
  annotation[[annotation_names["index"]]] <- annotation$id
  if(is.null(annotation[[annotation_names["facet"]]]))
    annotation[[annotation_names["facet"]]] <- ""
  if(!is.data.frame(annotation)) {
    # Could happen if was NULL originally
    annotation <- dplyr::as_tibble(annotation)
  }
  
  #   run mediation test and find the best models (using BIC among the four models)
  out <- intermediate::mediation_test(
    target = target,
    mediator = mediator,
    annotation = annotation,
    covar_tar = covar_tar,
    covar_med = covar_med,
    driver = driver,
    driver_med = driver_med,
    annotation_names = annotation_names, ...)
  
  class(out) <- c("mediation_index", class(out))
  attr(out, "annotation_names") <- annotation_names
  out
}
#' @export
#' @rdname mediation_index
plot.mediation_index <- function(x, ...)
  ggplot_mediation_index(x, ...)
#' @export
#' @rdname mediation_index
autoplot.mediation_index <- function(x, ...)
  ggplot_mediation_index(x, ...)
#' @param x object of class \code{mediation_index}
#' @param response type of response from \code{c("pvalue","IC")}
#' @param alpha ggplot2 transparency parameter with default \code{0.5}
#' @param pattern_name name of pattern with default \code{"pattern"}
#' @param ... additional parameters
#' 
#' @export
#' @rdname mediation_index
ggplot_mediation_index <- function(x, response = c("pvalue","IC"), alpha = 0.5,
                                   pattern_name = "pattern", ...) {
  response <- match.arg(response)
  index_name <- as.vector(attr(x, "annotation_names")["index"])
  
  ## Somehow index_name is not propagating through mediation_test.
  ## It seems cmst_default uses chr and pos. Fix earlier?
  if(index_name != "index" & "index" %in% names(x$best)) {
    # Make sure we don't clash with column named index.
    x$best$index <- NULL
  }
  best <-
    dplyr::rename(
      x$best,
      index = index_name)
  
  p <- ggplot2::ggplot(best)
  switch(
    response,
    pvalue = {
      p <- p +
        ggplot2::aes(.data$index, -log10(.data$pvalue))
      },
    IC     = {
      p <- p +
        ggplot2::aes(.data$index, .data$IC) +
        ggplot2::ylab("BIC on log10 scale")
      })
  if(pattern_name %in% names(best)) {
    p <- p + 
      ggplot2::aes(col = get(pattern_name)) +
      ggplot2::guides(color = ggplot2::guide_legend(pattern_name))
  }
  if(pattern_name != "pattern" & "pattern" %in% names(best))
    p <- p + ggplot2::aes(pattern = .data$pattern)
  p <- p +
    ggplot2::aes(group = .data$triad, id = .data$id)
  if(!is.null(target_index <- x$params$target_index)) {
    if(target_index >= min(best$index) & target_index <= max(best$index))
      p <- p +
        ggplot2::geom_vline(xintercept = target_index, col = "gray")
  }
  p <- p +
    ggplot2::geom_point(alpha = alpha) +
    ggplot2::facet_wrap(~ .data$triad) +
    ggplot2::xlab(index_name)
  if(!is.null(x$map)) {
    tmp <- 
      dplyr::filter(
        data.frame(map = x$map),
        map >= min(best$index),
        map <= max(best$index))
    p <- p +
      ggplot2::geom_rug(aes(.data$map), data = tmp, inherit.aes = FALSE, col = "gray")
    
  }
  p
}
