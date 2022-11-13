# Mediation tests
#
#' Develop mediation models from driver, target and mediator
#'
#' @param target vector or 1-column matrix with target values
#' @param mediator matrix of mediators
#' @param driver vector or matrix with driver values
#' @param annotation A data frame with mediators' annotation with columns for facet and index
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param driver_med optional driver matrix for mediators
#' @param intcovar optional interactive covariates (assumed same for `mediator` and `target`)
#' @param test Type of CMST test.
#' @param fitFunction function to fit models with driver, target and mediator
#' @param annotation_names names in annotation of columns for facet, index and, optionally, driver (default `c(facet = "chr", index = "pos", driver = NULL)`)
#' @param ... additional parameters
#'
#' @importFrom purrr map transpose
#' @importFrom stringr str_replace
#' @importFrom dplyr any_of arrange bind_rows desc filter group_by left_join
#' mutate one_of rename ungroup
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggplot2 aes autoplot element_blank facet_grid facet_wrap 
#' geom_hline geom_point geom_vline ggplot 
#' ggtitle scale_color_manual scale_shape_manual theme xlab ylab
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang .data
#' 
#' @return List with elements (names could be improved):
#' \describe{
#'   \item{best}{best triad via causal model selection test per id
#'   \describe{
#'     \item{id}{identifier for mediator}
#'     \item{symbol}{symbol for mediator, such as gene symbol (\code{symbol})}
#'     \item{facet}{facet for mediator, such as chromosome (\code{chr})}
#'     \item{index}{index of mediator within facet, such as position (\code{pos})}
#'     \item{triad}{best model: one of causal, reactive, independent, undecided}
#'     \item{pvalue}{p-value for causal model selection test}
#'     \item{alt_test}{best alternative for model test}
#'     \item{undecided}{p-value for best model vs undecided}
#'     \item{LR_mediation}{LR score for target on driver adjusted by mediator}
#'     \item{LR_mediator}{LR score for mediator on driver}
#'     \item{LR_CMST}{LR score for causal model selection test (see test below)}
#'     \item{df}{degrees of freedom for best test (see test below)}
#'     \item{IC}{information criterion (see test below)}
#'     \item{...}{optional columns from annotation}
#'   }}
#'   \item{test}{Four causal test results per id
#'   \describe{
#'     \item{id}{mediator identifier}
#'     \item{ref_test}{reference model being tested}
#'     \item{LR}{LR score for causal model selection test of ref_test against alt_test}
#'     \item{df}{difference in model degrees of freedom}
#'     \item{IC}{information criterion (determined by \code{flavor} option (default \code{"B"} for AIC, or \code{"A"} for AIC))}
#'     \item{pvalue}{p-value to compare ref_test with alt_test}
#'     \item{alt_test}{closest alternative for ref_test}
#'   }}
#'   \item{driver}{list of driver names for target and mediator(s)}
#'   \item{normF}{Frobenius norm if using both target and mediator drivers}
#'   \item{params}{list of parameter settings for use by summary and plot methods}
#' }
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
#' 
#' med_test <- mediation_test(target = target,
#'                       mediator = Tmem68$mediator[, med_signif, drop = FALSE],
#'                       driver = Tmem68$driver,
#'                       annotation = med_LR,
#'                       covar_tar = Tmem68$covar,
#'                       covar_med = Tmem68$covar)
#' summary(med_test)
#' 
#' ggplot_mediation_test(med_test)
#'
#' @export
#' @importFrom dplyr arrange bind_rows everything filter left_join matches
#'             mutate rename select
#' @importFrom purrr map transpose
#'
mediation_test <- function(target, mediator, driver, annotation = NULL,
                          covar_tar=NULL, covar_med=NULL,
                          driver_med = NULL, intcovar = NULL,
                          test = c("wilcoxon","binomial","joint","normal"),
                          fitFunction = fitDefault,
                          annotation_names = c(facet = "chr", index = "pos", driver = "qtl_pos"),
                          ...) {
  
  ## Need to enable different covariates for different mediators.
  
  # Get common data. **FIX THIS**
#  commons <- common_data(target, mediator, driver, covar_tar, covar_med,
#                         intcovar = intcovar)
  
  if(is.null(mediator))
    return(NULL)
  
  if(!is.null(driver)) {
    scan_max <- fitFunction(driver, target, covar_matrix(covar_tar), ...)
  } else {
    scan_max <- NULL
  }
  LR_target <- scan_max$LR
  
  test <- match.arg(test)
  testFunction <- switch(test,
                   wilcoxon = wilcIUCMST,
                   binomial = binomIUCMST,
                   joint    = normJointIUCMST,
                   normal   = normIUCMST)
  
  if(!is.matrix(mediator)) {
    mediator <- as.matrix(mediator)
    if(is.null(colnames(mediator)))
      colnames(mediator) <- "mediator"
  }
  
  # Convert any blank driver names to V1, V2, ...
  driver <- driver_blank_names(driver)
  driver_med <- driver_blank_names(driver_med)

  # Return list of test, fit, fitsLR from causal model selection tests by mediator.
  result <- mediation_test_internal(target, mediator, driver, annotation,
                                    covar_tar, covar_med,
                                    driver_med, intcovar,
                                    fitFunction, testFunction,
                                    cmst_default,
                                    ...)
  if(is.null(result))
    return(NULL)
  
  facet_name <- as.vector(annotation_names["facet"])
  index_name <- as.vector(annotation_names["index"])
    
  # Transpose result. Make sure elemnts of result are data frames.
  result <-
    purrr::map(
      purrr::transpose(result),
      function(x) {
        if(is.data.frame(x[[1]])) {
          dplyr::bind_rows(x, .id = "id")
        } else {
          isnt <- !sapply(x, is.null)
          if(any(isnt))
            as.data.frame(t(as.data.frame(x[isnt])))
          else
            NULL
        }
      })
  
  if(!is.null(result$normF) && all(is.na(result$normF)))
    result$normF <- NULL
  
  # Driver names are names of the drivers if separate by mediator.
  # This only arises if there are separate drivers for target and mediator(s).
  result$driver_names <- {
    if(is.null(driver_med)) NULL
    else {
      if(is.array(driver_med)) dimnames(driver_med)[[3]]
      else names(driver_med)
    }
  }
  
  # Driver levels are levels of the driver.
  result$driver_levels <- {
    if(!is.null(driver))
      colnames(driver)
    else {
      if(is.null(driver_med)) NULL
      else {
        if(is.array(driver_med)) colnames(driver_med)
        else colnames(driver_med[[1]])
      }
    }
  }
  
  # Fix up annotation so it has id column and renames LR to LR_annot if present.
  # Note that LR_annotation should agree with LR_mediator
  if(is.null(annotation))
    annotation <- data.frame(id = colnames(mediator),
                             stringsAsFactors = FALSE)
  if("LR" %in% names(annotation)) {
    annotation <- dplyr::rename(annotation, LR_annotation = "LR")
  }
  
  decided <- dplyr::filter(
    result$test,
    .data$ref_test != "undecided")
  
  # Best mediator causal model among causal, reactive, independent.
  # Separately report undecided p-value (as "undecided" column).
  result$best <-
    dplyr::select(
      dplyr::arrange(
        dplyr::rename(
          dplyr::left_join(
            # Join causal model tests ($test) with mediation and mediator LR ($fit).
            # result$test information
            dplyr::left_join(
              # join undecided pvalue with other tests as its own column ("undecided").
              dplyr::rename(
                dplyr::select(
                  dplyr::filter(
                    result$test,
                    .data$ref_test == "undecided"),
                  .data$id, .data$pvalue),
                undecided = "pvalue"),
              dplyr::left_join(
                # Join best test among decided with annotation.
                dplyr::bind_rows(
                  purrr::map(
                    split(
                      decided,
                      decided$id),
                      # Pick best of decided models.
                      # In case of ties, pick in order causal, reactive, independent.
                      function(x) x[which.min(x$pvalue)[1],, drop = FALSE])),
                annotation,
               by = "id"),
              by = "id"),
            # result$fit information
            # Mediation and mediator LR information from result$fit.
            dplyr::rename(
              tidyr::pivot_wider(
                dplyr::filter(
                  dplyr::select(
                    result$fit,
                    .data$id, .data$response, .data$LR),
                  .data$response %in% c("mediation", "mediator")),
                names_from = "response", values_from = "LR"),
              LR_mediation = "mediation", # LR for target on driver given mediator
              LR_mediator = "mediator"), # LR for mediator on driver
            by = "id"),
          LR_CMST = "LR",
          triad = "ref_test"),
        .data$pvalue),
      .data$id, 
      dplyr::matches("symbol"),
      dplyr::matches(facet_name),
      dplyr::matches(index_name),
      .data$triad, .data$pvalue, .data$alt_test,
      .data$undecided, .data$LR_mediation, .data$LR_mediator, 
      dplyr::everything())
  
  result$params <-
    list(LR_target = LR_target,
         target_name = colnames(target),
         test = test)
  result$targetFit <- scan_max
  
  attr(result, "annotation_names") <- annotation_names
  class(result) <- c("mediation_test", class(result))
  result
}
mediation_test_internal <- function(target, mediator, driver, annotation,
                                    covar_tar, covar_med,
                                    driver_med, intcovar,
                                    fitFunction,
                                    testFunction,
                                    cmstfn = cmst_default,
                                    ...) {
                                    
  use_1_driver <- is.null(annotation$driver_names) | is.null(driver_med)
  if(use_1_driver & !is.null(driver_med))
    driver_med <- NULL
  
  # Get common data.
  commons <- common_data(target, mediator, driver,
                         covar_tar, NULL, intcovar = intcovar,
                         common = use_1_driver, ...)
  if(is.null(commons))
    return(NULL)
  
  target <- commons$target
  mediator <- commons$mediator
  driver <- commons$driver
  covar_tar <- commons$covar_tar
  intcovar <- commons$intcovar
  common <- commons$common
  rm(commons)
  
  if(!is.null(covar_med)) {
    m <- match(rownames(mediator), rownames(covar_med), nomatch = 0)
    covar_med <- covar_med[m,, drop = FALSE]
    # Convert to matrix.
    covar_med <- covar_matrix(covar_med)
    covar_med <- convert_matrix(covar_med,
                                paste0("covM", seq_len(ncol(covar_med))), 
                                rownames(mediator))
  }
  
  # If we have driver_med, reduce to same subset as others
  if(!is.null(driver_med)) {
    if(is.array(driver_med)) {
      m <- match(rownames(mediator), rownames(driver_med), nomatch = 0)
      driver_med <- driver_med[m,,, drop = FALSE]
    } else {
      if(is.list(driver_med)) {
        m <- match(rownames(mediator), rownames(driver_med[[1]]), nomatch = 0)
        driver_med <- lapply(driver_med, function(x) x[m,, drop = FALSE])
      } else {
        stop("driver_med is neither an array nor a list")
      }
    }
  }
  
  # Reorganize annotation and mediator data.
  # Need to make sure elements of mediator have same ids.
  if(!is.null(annotation)) {
    annotation <- dplyr::filter(
      annotation,
      .data$id %in% colnames(mediator))
    # Make sure annotation is in same order as mediator.
    m <- match(colnames(mediator), annotation$id)
    if(any(is.na(m))) {
      cat("mediator and annotation do not match\n", file = stderr())
      return(NULL)
    }
    driver_names <- annotation$driver_names
    if(!is.null(driver_names))
      driver_names <- driver_names[m]
  } else {
    driver_names <- NULL
  }
  if(is.null(driver_names))
    driver_names <- rep("", ncol(mediator))
  names(driver_names) <- colnames(mediator)

  # Workhorse: CMST on each mediator.
  mediator <- as.data.frame(mediator, make.names = FALSE)

  purrr::map(
    purrr::transpose(
      list(
        mediator = mediator,
        driver_names = driver_names)),
    cmstfn, driver, target, 
    covar_tar, covar_med,
    driver_med, intcovar,
    fitFunction, testFunction, common, ...)
}

#' @export
#' @param x object of class \code{mediation_test}
#' @param not_type biotypes to not include
#' @rdname mediation_test
#' 
subset.mediation_test <- function(x, not_type, ...) {
  if("biotype" %in% names(x$best)) {
    attrc <- class(x)
    x$best <- dplyr::filter(x$best, 
                            .data$biotype != not_type)
    class(x) <- attrc
  }
  x
}
#' @param object object of class \code{mediation_test}
#' @export
#' @rdname mediation_test
#' 
summary.mediation_test <- function(object, ..., lod = FALSE) {
  out <-
    dplyr::mutate(
      dplyr::arrange(
        object$best,
        pmin(.data$pvalue, .data$undecided), .data$id),
      pvalue = signif(.data$pvalue, 3),
      LR_mediation = signif(.data$LR_mediation, 3),
      LR_mediator = signif(.data$LR_mediator, 3))

  if(lod) {
    out <- 
      dplyr::rename(
        dplyr::mutate(
          out,
          mediation = .data$LR_mediation / log(10),
          LR_mediator = .data$LR_mediator / log(10)),
        lod_mediator = "LR_mediator")
  }
  
  dplyr::arrange(out, .data$pvalue)
}
driver_blank_names <- function(driver) {
  if(is.null(driver)) 
    return(NULL)
  
  if(is.array(driver)) {
    if(is.null(colnames(driver))) {
      colnames(driver) <- LETTERS[seq_len(ncol(driver))]
    }
    if(any((cold <- colnames(driver)) == "")) {
      colnames(driver)[cold] <- LETTERS[seq_len(ncol(driver))][cold]
    }
  }
  driver
}
