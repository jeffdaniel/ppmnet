#' Log (composite) likelihood for a 'ppmnet' object
#'
#' Extracts the log likelihood for a regularized Poisson point process model, or
#' the log composite likelihood for a regularized Gibbs point process model.
#'
#' @param object A fitted \code{ppmnet} object.
#' @param s Value(s) of the penalty tuning parameter at which log (composite)
#'        likelihood is to be computed. Default is the entire sequence used to
#'        fit the regularization path.
#' @param ... Ignored.
#'
#' @export
logLik.ppmnet <- function(object, s = NULL, ...) {

  # Extract relevant quantities
  method <- object$method
  subset <- object$subset
  Q <- object$Q[subset]
  Z <- is.data(Q)

  # Evaluate conditional intensity at quadrature points
  cif <- fitted.ppmnet(object, s = s, drop = TRUE)

  # Compute log-likelihood
  if (method == "mpl") {
    w <- object$w[subset]
    if (dim(cif)[2] > 1) {
      ll <- colSums(log(cif[Z, ])) - colSums(w * cif)
    } else {
      ll <- sum(log(cif[Z])) - sum(w * cif)
    }
  } else if (method == "logi") {
    b <- rep(object$b, sum(subset))
    p <- cif / (b + cif)
    if (dim(cif)[2] > 1) {
      ll <- colSums(log(p/(1 - p))[Z, ]) + colSums(log(1 - p)) + sum(log(b[Z]))
    } else {
      ll <- sum(log(p/(1 - p))[Z]) + sum(log(1 - p)) + sum(log(b[Z]))
    }
  } else {
    stop("method = '", method, "' is unrecognized.", call. = FALSE)
  }
  return(ll)
}
