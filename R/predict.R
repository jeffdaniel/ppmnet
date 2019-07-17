#' Make predictions from a 'ppmnet' object
#'
#' Makes predictions from a regularized spatial point process model fit via
#' penalized composite likelihood.
#'
#' @param object A fitted \code{ppmnet} object.
#' @param data A list of pixel images (of class \code{imlist})
#'        containing the spatial covariates used to fit the model.
#' @param window Optional. An observation window (of class \code{owin})
#'        defining the region within which predictions are to be made. Default
#'        is the window of the original data used to fit the model.
#' @param eps Optional. The height and width of pixels in the prediction
#'        image(s). A numeric value or numeric vector of length 2 specifying
#'        pixel dimensions in the \strong{x} and \strong{y} directions.
#'        Incompatible with \code{dimyx}.
#' @param dimyx Optional. The resolution of the prediction image(s). A numeric
#'        value or numeric vector of length 2 specifying the number of pixels in
#'        the \strong{y} and \strong{x} directions. Incompatible with
#'        \code{eps}.
#' @param s Value(s) of the penalty tuning parameter at which predictions are to
#'        be made. Default is the entire sequence used to fit the regularization
#'        path.
#' @param type Type of prediction required. Either "\code{trend}" for the
#'        spatial trend, "\code{intensity}" for the intensity (Poisson models
#'         only), or "\code{cif}" for the conditional intensity.
#' @param ... Additional arguments passed to \code{predict.glmnet}.
#'
#' @return A list of pixel images containing predictions, or, if \code{s} is of
#'         length 1, a single pixel image.
#'
#' @examples
#' # Predicted intensities
#' Qp <- quadscheme(Xp)
#' fitp <- ppmnet(Qp, exdata)
#' predict(fitp, exdata)
#'
#' # Predicted conditional intensities
#' Qs <- quadscheme(Xs)
#' fits <- ppmnet(Qs, exdata, Strauss(5), nlambda = 20)
#' predict(fits, exdata, type = "cif")
#'
#' @export
predict.ppmnet <- function(object, data, window = NULL,
                           eps = NULL, dimyx = NULL, s = NULL,
                           type = c("trend",  "intensity", "cif"), ...) {

  # Validate covariate data
  if (!inherits(data, "imlist")) {
    stop("The argument 'data' must be a list of pixel images.",
         call. = FALSE)
  }

  # Validate window
  if (!is.null(window) && !is.owin(window)) {
    stop("The argument 'window' must be an observation window of class 'owin'.",
         call. = FALSE)
  }

  # Get window for predictions
  if (is.null(window)) {
    masque <- as.mask(object$Q$data$window, eps, dimyx)
  } else {
    masque <- as.mask(window, eps, dimyx)
  }

  # Get prediction points
  rxy <- rasterxy.mask(masque, drop = TRUE)
  xpredict <- rxy$x
  ypredict <- rxy$y

  # Construct matrix of covariate values at prediction points
  newx <- lapply(data, lookup.im, xpredict, ypredict,
                 naok = TRUE, strict = FALSE)
  newx <- do.call(cbind, newx)

  type <- match.arg(type)
  if (!is.poisson(object$interaction)) {
    if (type == "trend") {
      # Interaction potential is ignored
      vnames <- object$vnames
      newv <- matrix(0, nrow(newx), length(vnames),
                     dimnames = list(NULL, vnames))
      newx <- cbind(newx, newv)
    } else if (type == "cif") {
      # Evaluate interaction potential at prediction points
      Q <- object$Q; X <- Q$data
      U <- ppp(x = xpredict, y = ypredict, window = masque)
      newv <- evalInteraction(X, U, equalpairs(U, X),
                              interaction = object$interaction,
                              correction = object$correction)
      dimnames(newv) <- list(NULL, object$vnames)
      newx <- cbind(newx, newv)
    }
  }

  # Make predictions
  pred <- predict.glmnet(object, newx = newx, s = s, type = "link", ...)
  pred <- exp(pred)

  # Return predictions as a list of pixel images
  output <- as.list(seq(ncol(pred)))
  names(output) <- colnames(pred)
  for (i in 1:ncol(pred)) {
    z <- as.im(masque)
    z[] <- as.vector(pred[,i])
    output[[i]] <- z
  }
  output <- as.imlist(output)
  if (length(output) == 1) {
    output <- output[[1]]
  }
  return(output)
}

#' Fitted conditional intensity from a 'ppmnet' object
#'
#' Computes the fitted conditional intensity for a regularized spatial point
#' process model at the quadrature points used to fit the model.
#'
#' @param object A fitted \code{ppmnet} object.
#' @param s Value(s) of the penalty tuning parameter at which predictions are to
#'        be made. Default is the entire sequence used to fit the regularization
#'        path.
#' @param drop Logical value. If \code{TRUE}, quadrature points that were not
#'        used to fit the model are deleted.
#' @param ... Ignored
#'
#' @return A matrix of fitted values, or, if \code{s} is of length 1, a vector
#'         of fitted values.
#'
#' @examples
#' Qp <- quadscheme(Xp)
#' fitp <- ppmnet(Qp, exdata)
#' fitted(fitp)
#'
#' @export
fitted.ppmnet <- function(object, s = NULL, drop = FALSE, ...) {
  x <- object$x
  if (drop) {
    x <- x[object$subset, ]
  }
  cif <- predict.glmnet(object, newx = x, s = s, type = "link")
  cif <- exp(cif)
  return(cif)
}
