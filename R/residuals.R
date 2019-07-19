#' Residuals for a 'ppmnet' object
#'
#' Computes residual measures for a regularization path of point process models
#' fit via penalized composite likelihood.
#'
#' @param object A fitted \code{ppmnet} object.
#' @param s Value(s) of the penalty tuning parameter at which residuals are
#'        to be computed. Default is the entire sequence used to fit the
#'        regularization path.
#' @param type Type of residuals to be computed. Options are \code{"raw"},
#'        \code{"inverse"}, and \code{"pearson"}.
#' @param ... Ignored
#'
#' @return If \code{s} is of length 1, a measure (of class \code{msr});
#'         otherwise, a vector-valued measure (of class \code{msr}), with each
#'         component corresponding to a model in the  regularization path.
#'
#' @examples
#' Qp <- quadscheme(Xp)
#' fit <- ppmnet(Qp, exdata)
#' res <- residuals(fit)
#' smo <- Smooth(unstack(res)[[40]])
#' plot(smo)
#'
#' @importFrom stats fitted
#' @export
residuals.ppmnet <- function(object, s = NULL,
                             type = c("raw", "inverse", "pearson"), ...) {

  # Extract relevant values
  type <- match.arg(type)
  Q <- object$Q[object$subset]
  Z <- is.data(Q)
  cif <- fitted(object, s = s, drop = FALSE)
  ind <- cif > 0

  # Compute residual measures
  discrete <- switch(type,
                     raw     = rep.int(1, sum(Z)),
                     inverse = 1 / cif[Z, ],
                     pearson = 1 / sqrt(cif[Z, ]))
  density  <- switch(type,
                     raw     = -cif,
                     inverse = -ind,
                     pearson = -ind * sqrt(cif))
  res <- msr(Q, discrete, density)
  res
}
