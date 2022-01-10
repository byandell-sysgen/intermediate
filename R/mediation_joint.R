# Joint likelihood ratio for target and mediator
#
#' Fit joint likelihood for `target` and `mediator` given `driver`.
#'
#' @param target vector or 1-column matrix with target values
#' @param mediator matrix of mediators
#' @param driver vector or matrix with driver values
#' @param annotation A data frame with mediators' annotation with columns for `facet_name` and `index_name`
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param driver_med optional driver matrix for mediators
#' @param intcovar optional interactive covariates (assumed same for `mediator` and `target`)
#' @param fitFunction function to fit models with driver, target and mediator
#' @param annotation_names names in annotation of columns for facet, index and, optionally, driver (default `c(facet = "chr", index = "pos", driver = NULL)`)
#' @param ... additional parameters
#'
#' @importFrom purrr map transpose
#' @importFrom stringr str_replace
#' @importFrom dplyr arrange bind_rows desc filter group_by left_join mutate one_of rename ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 aes autoplot element_blank facet_grid facet_wrap 
#' geom_hline geom_point geom_vline ggplot 
#' ggtitle scale_color_manual scale_shape_manual theme xlab ylab
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @importFrom RColorBrewer brewer.pal
#' 
#' @return Data frame with `id` and `LR` as well as `annotation` columns.
#' 
#' @examples
#' data(Tmem68)
#'  
#' target <- Tmem68$target
#' 
#' # Find mediators with significant effect
#' med_LR <- mediator_LR(mediator = Tmem68$mediator,
#'                         driver = Tmem68$driver,
#'                         annotation = Tmem68$annotation,
#'                         covar_med = Tmem68$covar)
#' med_signif <- med_LR$id[med_LR$LR >= 5 * log(10)]
#' # Add info column.
#' med_LR$info <- paste("chr =", med_LR$chr)
#' # Rename mediator LR column to not conflict with mediation test LR.
#' med_LR <- dplyr::rename(med_LR, mediatorLR = "LR")
#' 
#' med_joint <- mediation_joint(target = target,
#'                       mediator = Tmem68$mediator[, med_signif, drop = FALSE],
#'                       driver = Tmem68$driver,
#'                       annotation = med_LR,
#'                       covar_tar = Tmem68$covar,
#'                       covar_med = Tmem68$covar)
#'                       
#' ggplot_mediation_joint(med_joint)
#' 
#' @export
#'
mediation_joint <- function(target, mediator, driver, annotation,
                          covar_tar=NULL, covar_med=NULL,
                          driver_med = NULL, intcovar = NULL,
                          fitFunction = fitDefault,
                          annotation_names = c(index = "pos"),
                          ...) {
  
  if(is.null(mediator))
    return(NULL)
  
  # If only one mediator, replicate it to match annotation.
  mediator <- as.matrix(mediator)
  if(ncol(mediator) == 1) {
    mediator <- mediator[, rep(1, nrow(annotation)), drop = FALSE]
    colnames(mediator) <- annotation$id
  }
  
  result <- mediation_test_internal(target, mediator, driver, annotation,
                                    covar_tar, covar_med,
                                    driver_med, intcovar,
                                    fitFunction, NULL,
                                    fit_joint,
                                    ...)
    
  out <- dplyr::left_join(
    dplyr::bind_rows(
      result,
      .id = "id"),
    annotation,
    by = "id")
  attr(out, "annotation_names") <- annotation_names
  class(out) <- c("mediation_joint", class(out))
  
  out
}

fit_joint <- function(object, driver, target, 
                      covar_tar, covar_med,
                      driver_med, intcovar,
                      fitFunction, testFunction,
                      common = TRUE, 
                      ...) {
  
  # Make sure we have driver or driver_med.
  driver_med <- get_driver_med(driver_med, object)
  if(is.null(driver)) {
    if(!is.null(driver_med))
      driver <- driver_med
    else {
      stop("must supply driver or driver_med")
    }
  }
  
  # Force x (= mediator column) to be matrix.
  mediator <- as.matrix(object[[1]])
  colnames(mediator) <- "mediator"
  rownames(mediator) <- rownames(driver)
  
  # Fit models
  fits <- med_fits(driver, target, mediator, fitFunction,
                   covar_tar, covar_med, driver_med,
                   intcovar, common = common,
                   fit_list = c("m.d_m","t.md_t"), ...)
  data.frame(LR = sum(fits$LR))
}

#' @export
#' @rdname mediation_joint
plot.mediation_joint <- function(x, ...)
  ggplot_mediation_joint(x, ...)
#' @export
#' @rdname mediation_joint
autoplot.mediation_joint <- function(x, ...)
  ggplot_mediation_joint(x, ...)
#' @param x object of class \code{mediation_joint}
#' @param lod plot lod if \code{TRUE}
#' @param xlab horizontal label
#' @param ylab vertical label
#' @param ... additional parameters
#' 
#' @export
#' @rdname mediation_joint
ggplot_mediation_joint <- function(x, lod = FALSE,
                                   xlab = index_name, ylab = ylab_name, ...) {
  index_name <- as.vector(attr(x, "annotation_names")["index"])
  if(index_name != "index" & "index" %in% names(x)) {
    # Make sure we don't clash with column named index.
    x$index <- NULL
  }
  x <- dplyr::rename(x, index = index_name)
  if(lod) {
    x <- dplyr::mutate(x, LR = .data$LR / log(10))
    ylab_name <- "LOD"
  } else {
    ylab_name <- "LR"
  }
  p <- ggplot2::ggplot(x) +
    ggplot2::aes(.data$index, .data$LR)
  if("pattern" %in% names(x))
    p <- p + ggplot2::aes(col = .data$pattern)
  p +
    ggplot2::geom_point() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
}
