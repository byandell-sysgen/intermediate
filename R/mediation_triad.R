#' Triad scatter plot for target on mediator given driver
#' 
#' Explores how target is explained by mediator and driver.
#' If mediator explains most of the driver effect on target,
#' then driver lines will be indistinguishable.
#' However, if driver has more information than the mediator
#' about the target, lines will be separated.
#' Parallel lines indicate driver and mediator have additive
#' effects on target; otherwise relationship is more complicated. 
#' 
#' @param target vector or 1-column matrix with target values
#' @param mediator vector or 1-column matrix with mediator values
#' @param driver vector or matrix with driver values
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param fitFunction function to fit models with driver, target and mediator
#' @param points,lines type of points and lines, from \code{c("dyad","monad")}
#' @param sdp optional sum across driver pattern to collapse columns of driver matrix
#' @param ... additional arguments
#' 
#' @details 
#' Plot method shows target (vertical axis) against mediator (horizontal axis)
#' with plot symbol identifying individuals (rows) by some summary of the driver.
#' Regression lines highlight effects of driver columns.
#' This works best when columns of driver are on the same scale, and even better
#' if the row values across column entries are positive summing to 1.
#' 
#' Plot label may be for column with largest row value (\code{monad}),
#' or a concatenation of labels for the two largest row values (\code{dyad}).
#' Lines may correspond to \code{monad} or to \code{sdp}. 
#' 
#' Driver columns may be compressed to 3 columns using \code{sdp} to 
#' Simplify Dyad Pattern
#' (idea taken from genetics \code{sdp} = SNP distribution pattern). 
#' Basically, the \code{sdp} recoded in base 2 separates columns of
#' driver into two groups: with 8 columns, \code{sdp} = 2 = 00000010 separate column 2 from
#' the other columns, while \code{sdp} = 5 = 00000101 separates columns 1 and 3 from the others.
#' 
#' Summary method shows coefficients from model fit of target on mediator and driver.
#' 
#' @examples
#' data(Tmem68)
#' 
#' target <- cbind(Tmem68 = Tmem68$target)
#' 
#' # Pick strongest mediator that is not target.
#' m <- match("Nnt", Tmem68$annotation$symbol)
#' mediator <- Tmem68$mediator[, m, drop = FALSE]
#' colnames(mediator) <- "Nnt"
#' 
#' med_triad <- mediation_triad(target = target,
#'                       mediator = mediator,
#'                       driver = Tmem68$driver,
#'                       covar_tar = Tmem68$covar,
#'                       sdp = 2)
#' 
#' summary(med_triad)
#' 
#' ggplot_mediation_triad(med_triad, tname = "Tmem68", mname = "Nnt",
#'                        fitlines = "sdp")
#' 
#' colors <- c("gray", "blue", rep("gray", 6))
#' names(colors) <- LETTERS[1:8]
#' ggplot_mediation_triad(med_triad, tname = "Tmem68", mname = "Nnt", 
#'                        fitlines = "driver", col = colors)
#' 
#' @export
#' 
#' @importFrom stringr str_remove str_split
#' @importFrom ggplot2 aes autoplot facet_wrap geom_hline geom_line geom_smooth 
#' geom_text ggplot ggtitle scale_color_discrete xlab ylab
#' @importFrom broom tidy
#' @importFrom stats lm
#' @importFrom knitr kable
#' @importFrom rlang .data
#'  

# Want to refactor so that we specify point text (monad or dyad or ?)
# and lines (sdp or monad or dyad). What we have is too complicated.
# sdp = collapse columns of driver based on base 2
# monad = label of column with largest value
# dyad

mediation_triad <- function(target, mediator, driver,
                        covar_tar = NULL, covar_med = NULL,
                        fitFunction = fitDefault,
                        points = c("dyad", "monad"),
                        lines = ifelse(points != "dyad", points, "monad"),
                        ...) {
  
  # Convert any blank driver names to A, B, ...
  driver <- driver_blank_names(driver)
  
  # points and lines must by dyad, monad or a positive integer < 2^ncol(driver) 
  choices <- c("dyad", "monad")
  maxsdp <- 2 ^ ncol(as.matrix(driver))
  points = points[1]
  if(n <- pmatch(points, choices, nomatch = 0)) {
    points <- choices[n]
  } else {
    points <- as.integer(points)
    stopifnot(!is.na(points),
              points > 0,
              points < maxsdp)
  }
  if(n <- pmatch(lines, choices, nomatch = 0)) {
    lines <- choices[n]
  } else {
    lines <- as.integer(lines)
    stopifnot(!is.na(lines),
              lines > 0,
              lines < maxsdp)
  }
  
  # ****NEED TO REFACTOR BELOW ***
  # Fix triad_abline when no sex but interaction
  # Separate out summary by model, or create print option
  # Incorporate points and lines options
  
  # Only allow one mediator for triad.
  if(ncol(mediator) > 1) {
    mediator <- mediator[,1, drop = FALSE]
  }

  # Fit target and target|mediator models
  fit <- med_fits(driver, target, mediator,
                  fitFunction, covar_tar, covar_med, ...)
  
  dat <- triad_data(target, mediator, driver, 
                    covar_tar, covar_med, ...)
  
  for(i in c("t.d_t","t.md_t.m")) {
    tmp <- fit$coef[[i]][seq_len(ncol(driver))]
    dat[[i]] <- c(as.matrix(dat[names(tmp)]) %*% tmp)
  }
  
  # Need to account for covariates and sex.
  out <- list(data = dat,
              coef = fit$coef[["t.d_t"]], coef_med = fit$coef[["t.md_t.m"]],
              drivers = colnames(driver), med_name = colnames(mediator))
  
  class(out) <- c("mediation_triad", class(out))
  
  out
}
triad_data <- function(target, mediator, driver, 
                       covar_tar, covar_med,
                       sdp = NULL,
                       label_fn = pattern_label,
                       group_fn = pattern_sdp,
                       ...) {
  
  # Find common data.
  commons <- common_data(target, mediator, driver, 
                         covar_tar, covar_med)
  
  # Get covariations from covar_med that are not in covar_tar
  if(!is.null(covar_med)) {
    cov_names <- colnames(covar_med)[!(colnames(covar_med) %in% colnames(covar_tar))]
    commons$covar_med <- commons$covar_med[,cov_names, drop = FALSE]
  } else {
    commons$covar_med <- matrix(NA, length(commons$target), 0)
  }
  
  # Set names for target and mediator (can change in plot)
  for(i in c("target","mediator"))
    colnames(commons[[i]]) <- i
  
  # Set up point labels and groups.
  if(is.null(label_fn))
    label_fn <- function(driver)
      toupper(substr(colnames(driver), 1, 1))[apply(driver, 1, function(x) which.max(x)[1])]
  label <- label_fn(commons$driver)
  # Group labels.
  if(is.null(group_fn))
    group_fn = function(label, a, b) label
  group <- group_fn(label, sdp, colnames(commons$driver))
  
  dat <- data.frame(commons$driver, commons$target, commons$mediator)
  if(length(commons$covar_tar))
    dat <- data.frame(dat, commons$covar_tar)
  if(length(commons$covar_med))
    dat <- data.frame(dat, commons$covar_med)
  dat <- data.frame(dat, label = label, group = group)
  
  if(!is.null(dat$sex))
    dat$Sex <- c("Female", "Male")[1 + dat$sex]

  dat
}
#' @param x object of class \code{summary.mediation_triad}
#' 
#' @rdname mediation_triad
#' @export
#' 
print.summary.mediation_triad <- function(x, ...) {
  cat("Model fit with driver columns")
  print(knitr::kable(
    dplyr::select(
      dplyr::filter(x,
                    .data$model == "allele"),
      -.data$model)))
  cat("\nModel fit with driver groups")
  print(knitr::kable(
    dplyr::select(
      dplyr::filter(x,
                    .data$model == "driver"),
      -.data$model)))
}
#' @param object object of class \code{mediation_triad}
#' 
#' @rdname mediation_triad
#' @export
#' 
summary.mediation_triad <- function(object, ...) {
  out <- dplyr::bind_rows(
    driver = lm_tidy(object$data, "group"),
    allele = lm_tidy(object$data, paste(object$drivers, collapse = "+")),
    .id = "model")
  class(out) <- c("summary.mediation_triad", class(out))
  out
}
lm_tidy <- function(object, driver, inter = "+") {
  form <- formula(paste("target ~ 0 + mediator",
                        ifelse(match("Sex", colnames(object), nomatch = 0),
                               paste("* Sex", inter),
                               inter),
                        driver))
  dplyr::mutate_if(
    broom::tidy(stats::lm(form, object)),
    is.numeric,
    function(x) ifelse(is.na(x), 0, x))
}
triad_abline <- function(object, drivers,
  fitlines = c("driver","parallel","sdp", "sdp-parallel")) {
  fitlines = match.arg(fitlines)
  inter <- ifelse(fitlines %in% c("driver", "sdp"), "*", "+")
  
  fit <- lm_tidy(object, 
                 paste0("(", paste(drivers, collapse = "+"), ")"),
                 inter)
  if(length(Sexes <- grep("^Sex", fit$term))) { # Sex covariate
    fit_drivers <- fit$estimate[match(drivers, fit$term)]
    if(inter != "+") { # interaction
      Sexes2 <- Sexes[grep(":", fit$term[Sexes])]
      Sexes <- Sexes[-grep(":", fit$term[Sexes])]
      sexes <- stringr::str_remove(fit$term[Sexes], "^Sex")
    } else {
      sexes <- stringr::str_remove(fit$term[Sexes], "^Sex")
    }
    Slopes <- grep("mediator", fit$term)
    if(inter != "+") {
      # We have mediator, mediator:Sex, mediator:driver and mediator:driver:Sex
      Slopes2 <- stringr::str_split(fit$term[Slopes], ":")
      Slopes4 <- sapply(Slopes2, length) == 3
      Slopes3 <- sapply(Slopes2, function(x) length(grep("Sex", x))) > 0
      slopes4 <- fit$estimate[Slopes[Slopes4]] # mediator:driver:Sex
      slopes3 <- fit$estimate[Slopes[!Slopes3]][-1] # mediator:driver
      slopes2 <- fit$estimate[Slopes[Slopes3 & !Slopes4]] # mediator:Sex
      # mediator and sex
      slopes <- fit$estimate[Slopes[1]]
      slopes <- c(slopes, slopes + slopes2)
      # mediator and drive (and sex effect on driver)
      slopes3 <- c(slopes3, rep(slopes3, length(Sexes) - 1) + slopes4)
    } else {
      slopes <- fit$estimate[Slopes]
      # The `slopes[-1]` are offsets from `slopes[1]`.
      slopes <- c(slopes[1], slopes[1] + slopes[-1])
    }
    out <- data.frame(
      driver = rep(drivers, length(sexes)),
      Sex = rep(sexes, rep(length(fit_drivers), length(sexes))),
      intercept = rep(fit$estimate[Sexes], rep(length(drivers), length(sexes))) +
        rep(fit_drivers, length(sexes)),
      slope = rep(slopes, rep(length(fit_drivers), length(sexes))))
    # Adjust slope and intercept for interaction with `Sex`.
    if(inter != "+") {
      out$intercept <- out$intercept +
        c(rep(0, length(drivers)), fit$estimate[Sexes2])
      out$slope <- out$slope + slopes3
    }
  } else { # No Sex covariate
    fit_drivers <- fit$estimate[match(drivers, fit$term)]
    slopes <- fit$estimate[grep("mediator", fit$term)]
    slopes <- slopes[1] + slopes[-1]
    out <- data.frame(
      driver = drivers,
      intercept = fit_drivers,
      slope = rep(slopes, length(fit_drivers)))
  }
  out
}
ggplot_triad_abline <- function(data, drivers, fitlines, colors, dname,
                                scale_color = TRUE) {
  # Add fitted model line.
  dat <- triad_abline(data, drivers, fitlines)
  dat$col <- dat$driver
  drivers <- unique(dat$driver)
  if(length(drivers) == length(colors)) {
    if(!is.null(names(colors)))
      dat$col <- names(colors)[match(dat$col, drivers)]
    else
      names(colors) <- drivers
    dat$col <- factor(dat$col, names(colors))
  }
  ggplot2::geom_abline(
      ggplot2::aes(slope = .data$slope,
                   intercept = .data$intercept,
                   col = .data$col),
      data = dat)
}
#' @param x object of class \code{mediation_triad}
#' @param tname target name (default \code{"target"})
#' @param mname mediator name (default \code{"mediator"})
#' @param dname driver name (default \code{"driver"})
#' @param centerline horizontal line at value (default = \code{NULL}); set to \code{NA} for no line or \code{NULL} for mean
#' @param fitlines use driver-specific (\code{"driver"}--the default), parallel (\code{"parallel"}), or 3 SDP lines (\code{"sdp"}) if \code{sdp} in \code{\link{mediation_triad}} is not \code{NULL}
#' @param main main title (defautl \code{tname})
#' @param colors named colors to use if \code{fitlines} is \code{TRUE}
#' @param size size of text (default \code{2})
#' 
#' @rdname mediation_triad
#' @export
ggplot_mediation_triad <- function(x, 
                             tname = "target", mname = "mediator", dname = "driver",
                             centerline = NULL,
                             main = paste(tname, "by", mname, "and", dname),
                             colors = NULL,
                             size = 2,
                             fitlines = c("driver","parallel","sdp",
                                          "sdp-parallel"),
                             ...) {
  
  fitlines = match.arg(fitlines)
  
  p <- ggplot2::ggplot(x$data) +
    ggplot2::ggtitle(main)
  
  if("label" %in% names(x$data)) {
    p <- p + 
      ggplot2::aes(label = .data$label) +
      ggplot2::geom_text(size = size)
  } else {
    p <- p +
      ggplot2::geom_point(alpha = 0.2)
  }
  
  if("Sex" %in% names(x$data)) {
    p <- p +
      ggplot2::facet_wrap(~ .data$Sex)
  }
  
  # set up mediator and target.
  p <- p + 
    ggplot2::aes(.data$mediator, .data$target) +
    ggplot2::xlab(mname) +
    ggplot2::ylab(tname)

  if(is.null(centerline)) {
    centerline <- mean(x$data$target, na.rm = TRUE)
  }
  if(!is.na(centerline)) {
    p <- p +
      ggplot2::geom_hline(yintercept = centerline, linetype = "dashed",
                          col = "grey60")
  }

  if(fitlines %in% c("driver", "parallel")) {
    if(is.null(colors)) colors <- seq_along(unique(x$drivers))
    p <- p + ggplot_triad_abline(x$data, x$drivers, fitlines, colors, dname)
    if((length(x$drivers) == length(colors)))
      p <- p + ggplot2::scale_color_manual(name = dname, values = colors)
  } else { # `fitlines %in% c("sdp", "sdp-parallel")`
    # Set up data for SDP fit lines.
    data <- x$data
    data$driver <- data$group
    drivers <- sort(unique(data$driver))
    col <- seq_along(drivers)
    # Add dummy columns for SDP values.
    sdps <- model.matrix(~ driver, data)
    colnames(sdps) <- drivers
    sdps[,1] <- 1 - apply(sdps[,-1, drop=FALSE], 1, sum)
    data <- data.frame(data, sdps)
    
    p <- p + ggplot_triad_abline(data, drivers, fitlines, col, dname) +
      ggplot2::aes(col = .data$group) +
      ggplot2::scale_color_discrete(name = dname)
  }
  p
}
#' @param object object of class \code{mediation_triad}
#' 
#' @rdname mediation_triad
#' @export
autoplot.mediation_triad <- function(object, ...) {
  ggplot_mediation_triad(object, ...)
}
#' @rdname mediation_triad
#' @export
plot.mediation_triad <- function(x, ...) {
  ggplot_mediation_triad(x, ...)
}

