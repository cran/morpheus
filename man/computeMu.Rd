% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeMu.R
\name{computeMu}
\alias{computeMu}
\title{Compute mu}
\usage{
computeMu(X, Y, optargs = list())
}
\arguments{
\item{X}{Matrix of input data (size nxd)}

\item{Y}{Vector of binary outputs (size n)}

\item{optargs}{List of optional argument:
\itemize{
  \item 'jd_method', joint diagonalization method from the package jointDiag:
    'uwedge' (default) or 'jedi'.
  \item 'jd_nvects', number of random vectors for joint-diagonalization
    (or 0 for p=d, canonical basis by default)
  \item 'M', moments of order 1,2,3: will be computed if not provided.
  \item 'K', number of populations (estimated with rank of M2 if not given)
}}
}
\value{
The estimated normalized parameters as columns of a matrix mu of size dxK
}
\description{
Estimate the normalized columns mu of the beta matrix parameter in a mixture of
logistic regressions models, with a spectral method described in the package vignette.
}
\examples{
io <- generateSampleIO(10000, 1/2, matrix(c(1,0,0,1),ncol=2), c(0,0), "probit")
mu <- computeMu(io$X, io$Y, list(K=2)) #or just X and Y for estimated K

}
\seealso{
\code{multiRun} to estimate statistics based on mu,
  and \code{generateSampleIO} for I/O random generation.
}
