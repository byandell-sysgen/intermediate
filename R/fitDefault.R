#' Fit target relative to driver
#'
#' Fit a model for a target and get detailed results
#' about estimated coefficients and individuals contributions to the LOD score.
#' 
#' @details Perform regression under the null hypothesis of no driver,
#' and then taken as fixed and known with driver.
#'
#' @param driver A matrix of drivers, individuals x drivers
#' @param target A numeric vector of target values
#' @param addcovar An optional matrix of additive covariates.
#' @param intcovar An optional matrix of interactive covariates.
#' @param weights An optional vector of positive weights for the
#' individuals. As with the other inputs, it must have `names`
#' for individual identifiers.
#' @param ... additional parameters
#' 
#' @return A list containing
#' * `LR` - The overall likelihood ratio.
#' * `indLR` - Vector of individual contributions to the likelihood ratio.
#' * `df` - Model degrees of freedom.
#' Currently `weights` is ignored.
#' 
#' @export
#' @importFrom stats as.formula dnorm model.matrix
#' 
fitDefault <- function(driver,
                  target,
                  addcovar = NULL,
                  intcovar=NULL, weights=NULL,
                  ...) {
  
  no.na <- apply(as.matrix(target), 1, function(x) !any(is.na(x)))
  if(!is.null(driver))
    no.na <- no.na & apply(driver, 1, function(x) !any(is.na(x)))
  if(!is.null(addcovar))
    no.na <- no.na & apply(addcovar, 1, function(x) !any(is.na(x)))
  driver <- driver[no.na,]
  target <- cbind(target)[no.na,]
  addcovar <- addcovar[no.na,]
  intcovar <- intcovar[no.na,]
  weights <- weights[no.na]

  # Original code fit T|D,C but want T|D,C - T|C
  full <- fitDefault_internal(driver, target, addcovar, intcovar, weights, ...) 
  red  <- fitDefault_internal(NULL,   target, addcovar, intcovar, weights, ...) 
  
  # If LR is 0, then make sure individual contributions are 0.
  if(full$LR == 0)
    full$indLR <- rep(0, length(full$indLR))
  if(red$LR == 0)
    red$indLR <- rep(0, length(red$indLR))
  
  full$LR <- full$LR - red$LR
  full$indLR <- full$indLR - red$indLR
  full$df <- full$df - red$df
  full$RSS <- full$RSS - red$RSS
  full
}
fitDefault_internal <- function(driver,
                       target,
                       addcovar = NULL,
                       intcovar=NULL, weights=NULL,
                       ...,
                       tol = 1e-16) {
  
  # Construct design matrix X
  if(is.null(driver)) {
    driver <- matrix(1, nrow(as.matrix(target)), 1)
  }
  
  # Form model matrix from additive covariates.
  if(!is.null(addcovar)) {
    X <- addcovar
    if(is.null(dim(X))) {
      X <- as.matrix(X)
    }
  } else 
    X <- NULL
  
  # Add interactive covaraties.
  if(!is.null(intcovar)) {
    if(ncol(driver) > 1) {
      driverbyintcovar <- driver[,-1, drop = FALSE] * intcovar
      X <- cbind(driver, X, driverbyintcovar)
    } else {
      X <- cbind(driver, X)
    }
  } else {
    if(!is.null(addcovar)){
      X <- cbind(driver, X)
    } else {
      X <- driver
    }
  }
  
  # Calculate log likelihood components
  n <- length(target)
  qrX <- qr(X)
  dX <- qrX$rank
  resid <- qr.resid(qrX, target)
  RSS <- sum(resid ^ 2)
  if(RSS <= tol) {
    RSS <- 0
    LR <- 0
  } else {
    LR <- as.vector(- (n/2) * (log(RSS)))
  }

  list(LR = LR, #as.vector(- (n/2) * (1 + log(2 * pi) + log(RSS / n))),
       indLR = stats::dnorm(target, qr.fitted(qrX, target), sqrt(RSS / n), log = TRUE),
       coef = qr.coef(qrX, target),
       df = dX,
       RSS = RSS,
       resid = resid)
}

