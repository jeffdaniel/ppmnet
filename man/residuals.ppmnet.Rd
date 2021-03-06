% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residuals.R
\name{residuals.ppmnet}
\alias{residuals.ppmnet}
\title{Residuals for a 'ppmnet' object}
\usage{
\method{residuals}{ppmnet}(object, s = NULL, type = c("raw", "inverse",
  "pearson"), ...)
}
\arguments{
\item{object}{A fitted \code{ppmnet} object.}

\item{s}{Value(s) of the penalty tuning parameter at which residuals are
to be computed. Default is the entire sequence used to fit the
regularization path.}

\item{type}{Type of residuals to be computed. Options are \code{"raw"},
\code{"inverse"}, and \code{"pearson"}.}

\item{...}{Ignored}
}
\value{
If \code{s} is of length 1, a measure (of class \code{msr});
        otherwise, a vector-valued measure (of class \code{msr}), with each
        component corresponding to a model in the  regularization path.
}
\description{
Computes residual measures for a regularization path of point process models
fit via penalized composite likelihood.
}
\examples{
Qp <- quadscheme(Xp)
fit <- ppmnet(Qp, exdata)
res <- residuals(fit)
smo <- Smooth(unstack(res)[[40]])
plot(smo)

}
