#' AUC values for a 'ppmnet' object
#'
#' Computes the area under the receiver operating characteristic curve (AUC) for
#' a regularization path of point process models.
#'
#' @param X A fitted \code{ppmnet} object.
#' @param data A list of pixel images (of class \code{imlist})
#'        containing the spatial covariates used to fit the model.
#' @param s Value(s) of the penalty tuning parameter at which AUC is to
#'        be calculated. Default is the entire sequence used to fit the
#'        regularization path.
#' @param ... Additional arguments passed to \code{predict.ppmnet}.
#'
#' @importFrom stats ecdf
#'
#' @examples
#' Qp <- quadscheme(Xp)
#' fitp <- ppmnet(Qp, exdata, nlambda = 10)
#' auc(fitp, exdata)
#'
#' @export
auc.ppmnet <- function(X, data, s = NULL, ...) {

  # Get point data
  Y <- X$Q$data

  # Compute predictions
  preds <- predict.ppmnet(X, data, s = s, ...)

  # Compute AUCs
  if (is.null(s) || length(s) > 1) {
    aucs <- numeric(length(preds))
    for (i in 1:length(aucs)) {
      Fl   <- ecdf(preds[[i]][])
      lamY <- preds[[i]][Y]
      aucs[i] <- mean(Fl(lamY))
    }
  } else {
    Fl <- ecdf(preds[])
    lamY <- preds[Y]
    aucs <- mean(Fl(lamY))
  }

  # If model is stationary, AUC is 0.5
  nz <- predict.glmnet(X, s = s, type = "nonzero")
  aucs[which(sapply(nz, is.null))] <- 1/2

  return(aucs)
}
