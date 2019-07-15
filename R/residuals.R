#' Residuals for a 'ppmnet' object
#'
#' Computes residual measures for a regularization path of point process models
#' fit via penalized composite likelihood.
#'
#' @param object A fitted \code{ppmnet} object.
#' @param s Value(s) of the penalty tuning parameter at which residuals are
#'        to be computed. Default is the entire sequence used to fit the
#'        regularization path.
#' @param ... Ignored
#'
#' @return If \code{s} is of length 1, a measure (of class \code{msr});
#'         otherwise, a vector-valued measure (of class \code{msr}), with each
#'         component corresponding to a model in the  regularization path.
#'
#' @examples
#' Qp <- spatstat::quadscheme(Xp)
#' fit <- ppmnet(Qp, exdata)
#' res <- residuals(fit)
#' smo <- spatstat::Smooth(unstack(res)[[40]])
#' plot(smo)
#'
#' @importFrom stats fitted
#' @export
residuals.ppmnet <- function(object, s = NULL, ...) {

  # Extract relevant values
  Q <- object$Q[object$subset]
  Z <- is.data(Q)
  cif <- fitted(object, s = s, drop = TRUE)

  # Compute residual measures
  discrete <- rep.int(1, sum(Z))
  density  <-  -cif
  res <- msr(Q, discrete, density)
  res
}
