#' @export
select <- function(object, ...) {
  UseMethod("select")
}

#' Select a tuning parameter value from a 'ppmnet' object
#'
#' Selects a point from along the fitted regularization path of a 'ppmnet'
#' object on the basis of a (composite) information criterion.
#'
#' @param object A fitted \code{ppmnet} object.
#' @param criterion The criterion by which to select the tuning parameter. One
#'        of \code{"AIC"}, \code{"BIC"}, or \code{"ERIC"}.
#' @param ... Ignored
#'
#' @return A list containing \code{lambda}, the selected tuning parameter;
#'         \code{coefs}, the coefficients of the selected model; \code{edf}, the
#'         effective degrees of freedom for the regularization path; and
#'         \code{criterion}; the computed criterion values for the
#'         regularization path.
#'
#' @examples
#' Qp <- spatstat::quadscheme(Xp)
#' fitp <- ppmnet(Qp, exdata)
#' select(fitp, "BIC")
#'
#' Qg <- spatstat::quadscheme.logi(Xg)
#' fitg <- ppmnet(Qg, exdata, interaction = spatstat::Geyer(5, 1),
#'                method = "logi", nlambda = 20)
#' select(fitg, "BIC")
#'
#' @aliases select
#' @importFrom stats as.formula
#' @export
select.ppmnet <- function(object, criterion = c("AIC", "BIC", "ERIC"), ...) {

  # Compute / extract required quantities
  ll <- logLik.ppmnet(object)
  lambda <- object$lambda
  m <- sum(is.data(object$Q))
  edf <- object$df
  if (!is.poisson(object$interaction)) {
    edf <- compute_edf(object)
  }

  # Compute (composite) information criterion
  criterion <- match.arg(criterion)
  penalty <- switch(criterion,
                    AIC  = 2 * edf,
                    BIC  = log(m) * edf,
                    ERIC = log(m / lambda) * edf)
  ic <- (-2 * ll) + penalty

  # Extact minimum
  ind <- which.min(ic)
  if (object$lambda[ind] == min(object$lambda)) {
    warning("minimum value of tuning parameter selected", call. = FALSE)
  }
  if (object$lambda[ind] == max(object$lambda)) {
    warning("maximum value of tuning parameter selected", call. = FALSE)
  }

  return(list(lambda = object$lambda[ind],
              coefs = coef.ppmnet(object)[, ind],
              edf = edf, criterion = ic))
}

compute_edf <- function(object, s = NULL, ...) {

  # Extract quadrature scheme and covariate matrix
  Q <- object$Q
  x <- object$x

  # Extract active set for each value of the tuning parameter
  nz <- predict.glmnet(object, s = s, type = "nonzero")

  # Identify terms associated with the interaction
  covnames <- colnames(x)
  vnames <- object$vnames
  inters <- which(covnames %in% vnames)

  # Convert covariate matrix into a data.frame for use in 'ppm' calls below
  covdf <- as.data.frame(x[, -inters])

  # For Poisson models, edf == df
  edf <- object$df

  # Follow regularization path, stopping if the model is non-Poisson
  state <- list()
  cat("Computing effective degrees of freedom for Gibbs models along",
      "regularization path...\n")
  for (i in 1:length(nz)) {
    if (inters %in% nz[[i]]) {

      # Extract coefficients and create trend formula
      coefs <- coef.ppmnet(object)[, i][c(1, nz[[i]] + 1)]
      coefnames <- names(coefs)
      trendvars <- coefnames[-which(coefnames %in% c("(Intercept)", vnames))]
      if (length(trendvars) == 0) {
        # Stationary: no spatial trend
        trend <- as.formula("~ 1")
      } else {
        trend <- paste(trendvars, collapse = " + ")
        trend <- as.formula(paste("~", trend))
      }

      # Construct 'ppm' object
      obj <- ppm(Q, trend, covariates = covdf,
                 interaction = object$interaction,
                 correction = object$correction,
                 method = object$method)

      #  Perform variance-covariance calculations
      vv <- vcov.ppm(obj, new.coef = coefs, what = "internals")

      # Get effective degrees of freedom
      if (object$method == "logi") {
        J <- vv$Sigma1log + vv$Sigma2log
        H <- vv$Slog
      } else {
        J <- vv$Sigma
        H <- vv$A1
      }
      JiH <- try(solve(H, J))
      edf[i] <- sum(diag(JiH))
    }
    state <- progressreport(i, length(nz), state = state)
  }
  edf
}
