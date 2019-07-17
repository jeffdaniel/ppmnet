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

#' ROC curves for a 'ppmnet' object
#'
#' Computes receiver operating characteristic (ROC) curves for a
#' regularization path of point process models.
#'
#' @param X A fitted \code{ppmnet} object.
#' @param data A list of pixel images (of class \code{imlist})
#'        containing the spatial covariates used to fit the model.
#' @param s Value(s) of the penalty tuning parameter at which the ROC curve is
#'        to be calculated. Default is the entire sequence used to fit the
#'        regularization path.
#' @param ... Additional arguments passed to \code{predict.ppmnet}.
#'
#' @examples
#' Qp <- quadscheme(Xp)
#' fitp <- ppmnet(Qp, exdata, nlambda = 10)
#' roc(fitp, exdata)
#'
#' @export
roc.ppmnet <- function(X, data, s = NULL, ...) {

  # Compute predictions and obtain null model
  preds <- predict.ppmnet(X, data, s = s, ...)
  nullmodel <- ppm(X$Q$data)

  if (length(s) == 1) {
    preds <- list(preds)
  }

  rocs <- list()
  for (i in 1:length(preds)) {
    D <- spatialCDFframe(nullmodel, preds[[i]])
    U <- D$values$U
    E <- ecdf(1 - U)
    p <- seq(0, 1, length = 1024)
    fobs <- E(p)
    df <- data.frame(p = p, fobs = fobs, fnull = p)
    rocs[[i]] <- fv(df,
                    argu = "p",
                    ylab = "ROC(p)",
                    valu = "fobs",
                    fmla = . ~ p,
                    desc = c("fraction of area",
                             "observed fraction of points",
                             "expected fraction if no effect"),
                    fname = "ROC")
  }
  if (length(rocs) == 1) {
    rocs <- unlist(rocs)
  }
  return(rocs)
}

