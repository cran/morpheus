\name{morpheus-package}
\alias{morpheus-package}
\alias{morpheus}
\docType{package}

\title{
	\packageTitle{morpheus}
}

\description{
	\packageDescription{morpheus}
}

\details{
	The package devtools should be useful in development stage, since we rely on testthat for
	unit tests, and roxygen2 for documentation. knitr is used to generate the package vignette.
	Concerning the other suggested packages:
	\itemize{
		\item{tensor is used for comparing to some reference functions initially coded in R;
			it should not be required in further package versions;}
		\item{jointDiag allows to solve a joint diagonalization problem, providing a more
			robust solution compared to a single diagonalization;}
		\item{parallel (generally) permits to run the bootstrap method faster.}
	}

	The three main functions are located in R/main.R:
	\itemize{
		\item{getParamsDirs_ref: reference method to estimate parameters directions;}
		\item{getParamsDirs: method of choice to estimate parameters directions, using a
			spectral decomposition of inputs/outputs;}
		\item{getBootstrapParams: run getParamsDirs on B bootstrap replicates.}
	}
}

\author{
	\packageAuthor{morpheus}

	Maintainer: \packageMaintainer{morpheus}
}

%\references{
%	TODO: Literature or other references for background information
%}

%\examples{
%	TODO: simple examples of the most important functions
%}
