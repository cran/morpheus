# extractParam
#
# Extract successive values of a projection of the parameter(s)
#
# @inheritParams plotHist
#
extractParam <- function(mr, x=1, y=1)
{
	# Obtain L vectors where L = number of res lists in mr
	lapply( mr, function(mr_list) {
		sapply(mr_list, function(m) m[x,y])
	} )
}

#' plotHist
#'
#' Plot histogram
#'
#' @param mr Output of multiRun(), list of lists of functions results
#' @param x Row index of the element inside the aggregated parameter
#' @param y Column index of the element inside the aggregated parameter
#'
#' @examples
#' \donttest{
#' beta <- matrix(c(1,-2,3,1),ncol=2)
#' mr <- multiRun(...) #see bootstrap example in ?multiRun : return lists of mu_hat
#' mu <- normalize(beta)
#' for (i in 1:2)
#'   mr[[i]] <- alignMatrices(res[[i]], ref=mu, ls_mode="exact")
#' plotHist(mr, 2, 1) #second row, first column}
#' @export
plotHist <- function(mr, x, y)
{
	params <- extractParam(mr, x, y)
	L = length(params)
	# Plot histograms side by side
	par(mfrow=c(1,L), cex.axis=1.5, cex.lab=1.5, mar=c(4.7,5,1,1))
	for (i in 1:L)
		hist(params[[i]], breaks=40, freq=FALSE, xlab="Parameter value", ylab="Density")
}

#' plotBox
#'
#' Draw boxplot
#'
#' @inheritParams plotHist
#'
#' @examples
#' #See example in ?plotHist
#' @export
plotBox <- function(mr, x, y)
{
	params <- extractParam(mr, x, y)
	L = length(params)
	# Plot boxplots side by side
	par(mfrow=c(1,L), cex.axis=1.5, cex.lab=1.5, mar=c(4.7,5,1,1))
	for (i in 1:L)
		boxplot(params[[i]], ylab="Parameter value")
}

#' plotCoefs
#'
#' Draw coefs estimations + standard deviations
#'
#' @inheritParams plotHist
#' @param params True value of parameters matrix
#'
#' @examples
#' #See example in ?plotHist
#' @export
plotCoefs <- function(mr, params)
{
	nf <- length(mr)
	L <- nrow(mr[[1]][[1]])
	K <- ncol(mr[[1]][[1]])

	params_hat <- vector("list", nf)
	stdev <- vector("list", nf)
	for (i in 1:nf)
	{
		params_hat[[i]] <- matrix(nrow=L, ncol=K)
		stdev[[i]] <- matrix(nrow=L, ncol=K)
	}
	for (x in 1:L)
	{
		for (y in 1:K)
		{
			estims <- extractParam(mr, x, y)
			for (i in 1:nf)
			{
				params_hat[[i]][x,y] <- mean(estims[[i]])
#				stdev[[i]][x,y] <- sqrt( mean( (estims[[i]] - params[x,y])^2 ) )
				# HACK remove extreme quantile in estims[[i]] before computing sd()
				stdev[[i]][x,y] <- sd( estims[[i]] [ estims[[i]] < max(estims[[i]]) & estims[[i]] > min(estims[[i]]) ] )
			}
		}
	}

	par(mfrow=c(1,nf), cex.axis=1.5, cex.lab=1.5, mar=c(4.7,5,1,1))
	params <- as.double(params)
	o <- order(params)
	for (i in 1:nf)
	{
		avg_param <- as.double(params_hat[[i]])
		std_param <- as.double(stdev[[i]])
		matplot(cbind(params[o],avg_param[o],avg_param[o]+std_param[o],avg_param[o]-std_param[o]),
			col=c(2,1,1,1), lty=c(1,1,2,2), type="l", lwd=2, xlab="param", ylab="value")
	}

	#print(o) #not returning o to avoid weird Jupyter issue... (TODO:)
}

#' plotQn
#'
#' Draw 3D map of objective function values
#'
#' @param N Number of starting points
#' @param n Number of points in sample
#' @param p Vector of proportions
#' @param b Vector of biases
#' @param beta Regression matrix (target)
#' @param link Link function (logit or probit)
#'
#' @export
plotQn <- function(N, n, p, beta, b, link)
{
	d <- nrow(beta)
	K <- ncol(beta)
	io <- generateSampleIO(n, p, beta, b, link)
	op <- optimParams(K, link, list(X=io$X, Y=io$Y))
	# N random starting points gaussian (TODO: around true beta?)
	res <- matrix(nrow=d*K+1, ncol=N)
	for (i in seq_len(N))
	{
		beta_init <- rnorm(d*K)
		par <- op$run( c(rep(1/K,K-1), beta_init, rep(0,K)) )
		par <- op$linArgs(par)
		Qn <- op$f(par)
		res[,i] = c(Qn, par[K:(K+d*K-1)])
	}
	res #TODO: plot this, not just return it...
}
