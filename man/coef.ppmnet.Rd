% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppmnet.R
\name{coef.ppmnet}
\alias{coef.ppmnet}
\title{Extract the coefficients from a 'ppmnet' object}
\usage{
\method{coef}{ppmnet}(object, s = NULL, ...)
}
\arguments{
\item{object}{A fitted \code{ppmnet} object.}

\item{s}{Value(s) of the penalty tuning parameter at which coefficients are
to be extracted. Default is the entire sequence used to fit the
regularization path.}

\item{...}{Additional arguments passed to \code{coef.glmnet}.}
}
\description{
Extract the coefficients from a 'ppmnet' object
}
