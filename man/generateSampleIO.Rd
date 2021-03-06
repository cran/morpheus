% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleIO.R
\name{generateSampleIO}
\alias{generateSampleIO}
\title{Generate sample inputs-outputs}
\usage{
generateSampleIO(n, p, beta, b, link)
}
\arguments{
\item{n}{Number of individuals}

\item{p}{Vector of K(-1) populations relative proportions (sum (<)= 1)}

\item{beta}{Vectors of model parameters for each population, of size dxK}

\item{b}{Vector of intercept values (use rep(0,K) for no intercept)}

\item{link}{Link type; "logit" or "probit"}
}
\value{
A list with
\itemize{
  \item{X: the input matrix (size nxd)}
  \item{Y: the output vector (size n)}
  \item{index: the population index (in 1:K) for each row in X}
}
}
\description{
Generate input matrix X of size nxd and binary output of size n, where Y is subdivided
into K groups of proportions p. Inside one group, the probability law P(Y=1) is
described by the corresponding column parameter in the matrix beta + intercept b.
}
\examples{
# K = 3 so we give first two components of p: 0.3 and 0.3 (p[3] = 0.4)
io <- generateSampleIO(1000, c(.3,.3),
  matrix(c(1,3,-1,1,2,1),ncol=3), c(.5,-1,0), "logit")
io$index[1] #number of the group of X[1,] and Y[1] (in 1...K)

}
