#' Fit a spatial point process model with lasso or elastic net regularization
#'
#' Fits spatial point process models via penalized composite likelihood. The
#' regularization path is computed for the lasso or elastic net penalty along a
#' sequence of tuning parameter values. Support for Poisson point process models
#' and for Gibbs point process models.
#'
#' @param Q A quadrature scheme (of class \code{quad}) containing a
#'        point pattern.
#' @param data A list of pixel images (of class \code{imlist})
#'        representing spatial covariates.
#' @param interaction An object (of class \code{interact}) describing the
#'        point process interaction structure, or \code{NULL} indicating that a
#'        Poisson point process model is to be fit.
#' @param correction Edge correction to be used. Either "\code{border}" for
#'        border correction, or "\code{none}" for no correction.
#' @param rbord If \code{correction = "border"}, the distance by which the
#'        observation window is eroded.
#' @param method Method to be used to fit the model. Either
#'        "\code{mpl}" for penalized maximum pseudolikelihood, or "\code{logi}"
#'        for penalized maximum logistic composite likelihood.
#' @param ... Additional arguments passed to \code{glmnet} to control the
#'        fitting procedure.
#'
#' @return An object of class "\code{ppmnet}".
#'
#' @examples
#' # Poisson model fit via penalized maximum likelihood
#' Qp <- spatstat::quadscheme(Xp)
#' fitp <- ppmnet(Qp, exdata)
#'
#' # Strauss model fit via penalized maximum pseudolikelihood
#' Qs <- spatstat::quadscheme(Xs)
#' fits <- ppmnet(Qs, exdata, interaction = spatstat::Strauss(5), nlambda = 20)
#'
#' # Geyer saturation model fit via penalized logistic composite likelihood
#' Qg <- spatstat::quadscheme.logi(Xg)
#' fitg <- ppmnet(Qg, exdata, interaction = spatstat::Geyer(5, 1),
#'                method = "logi", nlambda = 20)
#'
#' @import spatstat
#' @importFrom glmnet glmnet
#' @importFrom glmnet deviance.glmnet
#' @importFrom glmnet predict.glmnet
#' @importFrom glmnet plot.glmnet
#' @export
#'
ppmnet <- function(Q, data, interaction = NULL, correction = "border",
                   rbord = NULL, method = "mpl", ...) {

  # ----

  # Validate point pattern data
  if (!inherits(Q, "quad")) {
    stop("The argument 'Q' must be a quadrature scheme.", call. = FALSE)
  }
  if (is.marked(Q) || is.multitype(Q)) {
   stop("Marked and multitype patterns are not currently supported.",
        call. = FALSE)
  }

  # Validate covariate data
  if (!inherits(data, "imlist")) {
    stop("The argument 'data' must be a list of pixel images.", call. = FALSE)
  }
  if (length(data) < 2) {
    stop("The argument 'data' must be a list of 2 or more pixel images.",
         call. = FALSE)
  }

  # Validate interaction
  if (is.null(interaction)) {
    interaction <- Poisson()
  }
  if (!inherits(interaction, "interact")) {
    stop("The argument 'interaction' must be an object of class 'interact'.",
         call. = FALSE)
  }
  if (correction == "border") {
    if (is.null(rbord)) {
      rbord <- reach(interaction)
    }
    if (is.infinite(rbord) || eroded.areas(as.owin(Q), rbord) == 0) {
      stop("The reach of this interaction is too large for this window. ",
           "Please specify a smaller value of 'rbord'.", call. = FALSE)
    }
  } else {
    if (correction != "none") {
      stop("The argument 'correction' must be one of: 'border' or 'none'.",
           call. = FALSE)
    }
  }

  # Validate fitting method
  if (method == "logi") {
    if (!inherits(Q, "logiquad")) {
      stop("method = 'logi' only makes sense when 'Q' is of class 'logiquad'",
           call. = FALSE)
    }
  } else if (method != "mpl") {
      stop("The argument 'method' must be one of: 'mpl' or 'logi'.",
         call. = FALSE)
  }

  # ----

  # Construct matrix of covariate values at quadrature points
  U <- union.quad(Q)
  x <- lapply(data, lookup.im, U$x, U$y, naok = TRUE, strict = FALSE)
  x <- do.call(cbind, x)

  # Evaluate interaction potential at quadrature points, if applicable
  if (!is.poisson(interaction)) {
    v <- evalInteraction(Q$data, U, equalpairs.quad(Q), interaction, correction)
    vnames <- dimnames(v)[[2]]
    if (is.null(vnames)) {
      nc <- ncol(v)
      if (nc == 1) {
        vnames <- "Interaction"
      } else {
        vnames <- paste("Interact.", 1:nc, sep = "")
      }
      dimnames(v) <- list(dimnames(v)[[1]], vnames)
    }
    x <- cbind(x, v)
  } else {
    vnames <- NULL
  }

  # Construct vectors of "responses" and quadrature weights
  n <- n.quad(Q)
  if (method == "mpl") {
    fam <- "poisson"
    w <- w.quad(Q)
    z <- is.data(Q)
    y <- numeric(n)
    y[z] <- 1 / w[z]
  }
  if (method == "logi") {
    fam <- "binomial"
    b <- Q$param$rho # Get offset for logistic fit
    w <- rep.int(1, n)
    y <- as.numeric(is.data(Q))
  }

  # Ignore any quadrature points with zero weight
  zeroes <- attr(w, "zeroes")
  subset <- if (is.null(zeroes)) rep.int(TRUE, n) else !zeroes

  # Ignore any quadrature points with missing covariate values
  subset <- subset & !(rowSums(is.na(x)) > 0)

  # Apply border correction, if applicable
  if (correction == "border") {
    bdp <- bdist.points(U)
    subset <- subset & (bdp >= rbord)
  }

  # Fit the regularization path
  fit <- glmnet(x = x[subset, ],
                y = y[subset],
                weights = w[subset],
                family = fam, ...)
  if (fam == "binomial") fit$a0 <- fit$a0 + log(b) # Hack to add logistic offset

  # ----

  # Results
  fit$Q <- Q
  fit$x <- x
  fit$y <- y
  fit$w <- w
  fit$subset <- subset
  fit$method <- method
  fit$vnames <- vnames
  fit$call <- match.call()
  fit$correction <- correction
  fit$interaction <- interaction
  if (fam == "binomial") fit$b <- b
  class(fit) <- c("ppmnet", class(fit))
  fit
}

#' Extract the coefficients from a 'ppmnet' object
#' @param object A fitted \code{ppmnet} object.
#' @param s Value(s) of the penalty tuning parameter at which coefficients are
#'        to be extracted. Default is the entire sequence used to fit the
#'        regularization path.
#' @param ... Additional arguments passed to \code{coef.glmnet}.
#' @export
coef.ppmnet <- function(object, s = NULL, ...) {
  predict.glmnet(object, s = s, type = "coefficients", ...)
}

#' Extract the deviance from a 'ppmnet' object
#' @param object A fitted \code{ppmnet} object.
#' @param ... Ignored.
#' @export
deviance.ppmnet <- function(object, ...) {
  deviance.glmnet(object, ...)
}

#' Plot the coefficients from a 'ppmnet' object
#' @param x A fitted \code{ppmnet} object.
#' @param ... Additional arugments passed to \code{plot.glmnet}.
#' @export
plot.ppmnet <- function(x, ...) {
  plot.glmnet(x, ...)
}

#' Print a 'ppmnet' object
#' @param x A fitted \code{ppmnet} object.
#' @param digits The number of significant digits used in the printout.
#' @param ... Ignored.
#' @export
print.ppmnet <- function(x, digits = max(3, getOption("digits") - 3),
                         ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Df = x$df, `%Dev` = signif(x$dev.ratio, digits),
        Lambda = signif(x$lambda, digits)))
}
