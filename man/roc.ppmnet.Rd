% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auc.R
\name{roc.ppmnet}
\alias{roc.ppmnet}
\title{ROC curves for a 'ppmnet' object}
\usage{
\method{roc}{ppmnet}(X, data, s = NULL, ...)
}
\arguments{
\item{X}{A fitted \code{ppmnet} object.}

\item{data}{A list of pixel images (of class \code{imlist})
containing the spatial covariates used to fit the model.}

\item{s}{Value(s) of the penalty tuning parameter at which the ROC curve is
to be calculated. Default is the entire sequence used to fit the
regularization path.}

\item{...}{Additional arguments passed to \code{predict.ppmnet}.}
}
\description{
Computes receiver operating characteristic (ROC) curves for a
regularization path of point process models.
}
\examples{
Qp <- quadscheme(Xp)
fitp <- ppmnet(Qp, exdata, nlambda = 10)
roc(fitp, exdata)

}
