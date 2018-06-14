#' Optimize parameters
#'
#' Optimize the parameters of a mixture of logistic regressions model, possibly using
#' \code{mu <- computeMu(...)} as a partial starting point.
#'
#' @param K Number of populations.
#' @param link The link type, 'logit' or 'probit'.
#' @param optargs a list with optional arguments:
#'   \itemize{
#'     \item 'M' : list of moments of order 1,2,3: will be computed if not provided.
#'     \item 'X,Y' : input/output, mandatory if moments not given
#'     \item 'exact': use exact formulas when available?
#'   }
#'
#' @return An object 'op' of class OptimParams, initialized so that \code{op$run(x0)}
#'   outputs the list of optimized parameters
#'   \itemize{
#'     \item p: proportions, size K
#'     \item beta: regression matrix, size dxK
#'     \item b: intercepts, size K
#'   }
#'   x0 is a vector containing respectively the K-1 first elements of p, then beta by
#'   columns, and finally b: \code{x0 = c(p[1:(K-1)],as.double(beta),b)}.
#'
#' @seealso \code{multiRun} to estimate statistics based on beta, and
#'   \code{generateSampleIO} for I/O random generation.
#'
#' @examples
#' # Optimize parameters from estimated mu
#' io = generateSampleIO(10000, 1/2, matrix(c(1,-2,3,1),ncol=2), c(0,0), "logit")
#' mu = computeMu(io$X, io$Y, list(K=2))
#' M <- computeMoments(io$X, io$Y)
#' o <- optimParams(2, "logit", list(M=M))
#' x0 <- c(1/2, as.double(mu), c(0,0))
#' par0 <- o$run(x0)
#' # Compare with another starting point
#' x1 <- c(1/2, 2*as.double(mu), c(0,0))
#' par1 <- o$run(x1)
#' o$f( o$linArgs(par0) )
#' o$f( o$linArgs(par1) )
#' @export
optimParams = function(K, link=c("logit","probit"), optargs=list())
{
	# Check arguments
	link <- match.arg(link)
	if (!is.list(optargs))
		stop("optargs: list")
	if (!is.numeric(K) || K < 2)
		stop("K: integer >= 2")

	M <- optargs$M
	if (is.null(M))
	{
		if (is.null(optargs$X) || is.null(optargs$Y))
			stop("If moments are not provided, X and Y are required")
		M <- computeMoments(optargs$X,optargs$Y)
	}

	# Build and return optimization algorithm object
	methods::new("OptimParams", "li"=link, "M1"=as.double(M[[1]]),
		"M2"=as.double(M[[2]]), "M3"=as.double(M[[3]]), "K"=as.integer(K))
}

# Encapsulated optimization for p (proportions), beta and b (regression parameters)
#
# @field li Link, 'logit' or 'probit'
# @field M1 Estimated first-order moment
# @field M2 Estimated second-order moment (flattened)
# @field M3 Estimated third-order moment (flattened)
# @field K Number of populations
# @field d Number of dimensions
#
setRefClass(
	Class = "OptimParams",

	fields = list(
		# Inputs
		li = "character", #link 'logit' or 'probit'
		M1 = "numeric", #order-1 moment (vector size d)
		M2 = "numeric", #M2 easier to process as a vector
		M3 = "numeric", #M3 easier to process as a vector
		# Dimensions
		K = "integer",
		d = "integer"
	),

	methods = list(
		initialize = function(...)
		{
			"Check args and initialize K, d"

			callSuper(...)
			if (!hasArg("li") || !hasArg("M1") || !hasArg("M2") || !hasArg("M3")
				|| !hasArg("K"))
			{
				stop("Missing arguments")
			}

			d <<- length(M1)
		},

		expArgs = function(x)
		{
			"Expand individual arguments from vector x"

			list(
				# p: dimension K-1, need to be completed
				"p" = c(x[1:(K-1)], 1-sum(x[1:(K-1)])),
				"beta" = matrix(x[K:(K+d*K-1)], ncol=K),
				"b" = x[(K+d*K):(K+(d+1)*K-1)])
		},

		linArgs = function(o)
		{
			" Linearize vectors+matrices into a vector x"

			c(o$p[1:(K-1)], as.double(o$beta), o$b)
		},

		f = function(x)
		{
			"Sum of squares (Mi - hat_Mi)^2 where Mi is obtained from formula"

			P <- expArgs(x)
			p <- P$p
			beta <- P$beta
			lambda <- sqrt(colSums(beta^2))
			b <- P$b

			# Tensorial products beta^2 = beta2 and beta^3 = beta3 must be computed from current beta1
			beta2 <- apply(beta, 2, function(col) col %o% col)
			beta3 <- apply(beta, 2, function(col) col %o% col %o% col)

			return(
				sum( ( beta  %*% (p * .G(li,1,lambda,b)) - M1 )^2 ) +
				sum( ( beta2 %*% (p * .G(li,2,lambda,b)) - M2 )^2 ) +
				sum( ( beta3 %*% (p * .G(li,3,lambda,b)) - M3 )^2 ) )
		},

		grad_f = function(x)
		{
			"Gradient of f, dimension (K-1) + d*K + K = (d+2)*K - 1"

			P <- expArgs(x)
			p <- P$p
			beta <- P$beta
			lambda <- sqrt(colSums(beta^2))
			mu <- sweep(beta, 2, lambda, '/')
			b <- P$b

			# Tensorial products beta^2 = beta2 and beta^3 = beta3 must be computed from current beta1
			beta2 <- apply(beta, 2, function(col) col %o% col)
			beta3 <- apply(beta, 2, function(col) col %o% col %o% col)

			# Some precomputations
			G1 = .G(li,1,lambda,b)
			G2 = .G(li,2,lambda,b)
			G3 = .G(li,3,lambda,b)
			G4 = .G(li,4,lambda,b)
			G5 = .G(li,5,lambda,b)

			# (Mi - hat_Mi)^2 ' == (Mi - hat_Mi)' 2(Mi - hat_Mi) = Mi' Fi
			F1 = as.double( 2 * ( beta  %*% (p * G1) - M1 ) )
			F2 = as.double( 2 * ( beta2 %*% (p * G2) - M2 ) )
			F3 = as.double( 2 * ( beta3 %*% (p * G3) - M3 ) )

			km1 = 1:(K-1)
			grad <- #gradient on p
				t( sweep(as.matrix(beta [,km1]), 2, G1[km1], '*') - G1[K] * beta [,K] ) %*% F1 +
				t( sweep(as.matrix(beta2[,km1]), 2, G2[km1], '*') - G2[K] * beta2[,K] ) %*% F2 +
				t( sweep(as.matrix(beta3[,km1]), 2, G3[km1], '*') - G3[K] * beta3[,K] ) %*% F3

			grad_beta <- matrix(nrow=d, ncol=K)
			for (i in 1:d)
			{
				# i determines the derivated matrix dbeta[2,3]

				dbeta_left <- sweep(beta, 2, p * G3 * beta[i,], '*')
				dbeta_right <- matrix(0, nrow=d, ncol=K)
				block <- i
				dbeta_right[block,] <- dbeta_right[block,] + 1
				dbeta <- dbeta_left + sweep(dbeta_right, 2,  p * G1, '*')

				dbeta2_left <- sweep(beta2, 2, p * G4 * beta[i,], '*')
				dbeta2_right <- do.call( rbind, lapply(1:d, function(j) {
					sweep(dbeta_right, 2, beta[j,], '*')
				}) )
				block <- ((i-1)*d+1):(i*d)
				dbeta2_right[block,] <- dbeta2_right[block,] + beta
				dbeta2 <- dbeta2_left + sweep(dbeta2_right, 2, p * G2, '*')

				dbeta3_left <- sweep(beta3, 2, p * G5 * beta[i,], '*')
				dbeta3_right <- do.call( rbind, lapply(1:d, function(j) {
					sweep(dbeta2_right, 2, beta[j,], '*')
				}) )
				block <- ((i-1)*d*d+1):(i*d*d)
				dbeta3_right[block,] <- dbeta3_right[block,] + beta2
				dbeta3 <- dbeta3_left + sweep(dbeta3_right, 2, p * G3, '*')

				grad_beta[i,] <- t(dbeta) %*% F1 + t(dbeta2) %*% F2 + t(dbeta3) %*% F3
			}
			grad <- c(grad, as.double(grad_beta))

			grad = c(grad, #gradient on b
				t( sweep(beta,  2, p * G2, '*') ) %*% F1 +
				t( sweep(beta2, 2, p * G3, '*') ) %*% F2 +
				t( sweep(beta3, 2, p * G4, '*') ) %*% F3 )

			grad
		},

		run = function(x0)
		{
			"Run optimization from x0 with solver..."

			if (!is.numeric(x0) || any(is.na(x0)) || length(x0) != (d+2)*K-1
				|| any(x0[1:(K-1)] < 0) || sum(x0[1:(K-1)]) > 1)
			{
				stop("x0: numeric vector, no NA, length (d+2)*K-1, sum(x0[1:(K-1) >= 0]) <= 1")
			}

			op_res = constrOptim( x0, .self$f, .self$grad_f,
				ui=cbind(
					rbind( rep(-1,K-1), diag(K-1) ),
					matrix(0, nrow=K, ncol=(d+1)*K) ),
				ci=c(-1,rep(0,K-1)) )

			expArgs(op_res$par)
		}
	)
)

# Compute vectorial E[g^{(order)}(<beta,x> + b)] with x~N(0,Id) (integral in R^d)
#                 = E[g^{(order)}(z)] with z~N(b,diag(lambda))
#
# @param link Link, 'logit' or 'probit'
# @param order Order of derivative
# @param lambda Norm of columns of beta
# @param b Intercept
#
.G <- function(link, order, lambda, b)
{
	# NOTE: weird "integral divergent" error on inputs:
	# link="probit"; order=2; lambda=c(531.8099,586.8893,523.5816); b=c(-118.512674,-3.488020,2.109969)
	# Switch to pracma package for that (but it seems slow...)

	exactComp <- FALSE #TODO: global, or argument...

	if (exactComp && link == "probit")
	{
		# Use exact computations
		sapply( seq_along(lambda), function(k) {
			.exactProbitIntegral(order, lambda[k], b[k])
		})
	}

	else
	{
		# Numerical integration
		sapply( seq_along(lambda), function(k) {
			res <- NULL
			tryCatch({
				# Fast code, may fail:
				res <- stats::integrate(
					function(z) .deriv[[link]][[order]](lambda[k]*z+b[k]) * exp(-z^2/2) / sqrt(2*pi),
					lower=-Inf, upper=Inf )$value
			}, error = function(e) {
				# Robust slow code, no fails observed:
				sink("/dev/null") #pracma package has some useless printed outputs...
				res <- pracma::integral(
					function(z) .deriv[[link]][[order]](lambda[k]*z+b[k]) * exp(-z^2/2) / sqrt(2*pi),
					xmin=-Inf, xmax=Inf, method="Kronrod")
				sink()
			})
			res
		})
	}
}

# TODO: check these computations (wrong atm)
.exactProbitIntegral <- function(order, lambda, b)
{
	c1 = (1/sqrt(2*pi)) * exp( -.5 * b/((lambda^2+1)^2) )
	if (order == 1)
		return (c1)
	c2 = b - lambda^2 / (lambda^2+1)
	if (order == 2)
		return (c1 * c2)
	if (order == 3)
		return (c1 * (lambda^2 - 1 + c2^2))
	if (order == 4)
		return ( (c1*c2/((lambda^2+1)^2)) * (-lambda^4*((b+1)^2+1) -
			2*lambda^3 + lambda^2*(2-2*b*(b-1)) + 6*lambda + 3 - b^2) )
	if (order == 5) #only remaining case...
		return ( c1 * (3*lambda^4+c2^4+6*c1^2*(lambda^2-1) - 6*lambda^2 + 6) )
}

# Derivatives list: g^(k)(x) for links 'logit' and 'probit'
#
.deriv <- list(
	"probit"=list(
		# 'probit' derivatives list;
		# TODO: exact values for the integral E[g^(k)(lambdaz+b)]
		function(x) exp(-x^2/2)/(sqrt(2*pi)),                     #g'
		function(x) exp(-x^2/2)/(sqrt(2*pi)) *  -x,               #g''
		function(x) exp(-x^2/2)/(sqrt(2*pi)) * ( x^2 - 1),        #g^(3)
		function(x) exp(-x^2/2)/(sqrt(2*pi)) * (-x^3 + 3*x),      #g^(4)
		function(x) exp(-x^2/2)/(sqrt(2*pi)) * ( x^4 - 6*x^2 + 3) #g^(5)
	),
	"logit"=list(
		# Sigmoid derivatives list, obtained with http://www.derivative-calculator.net/
		# @seealso http://www.ece.uc.edu/~aminai/papers/minai_sigmoids_NN93.pdf
		function(x) {e=exp(x); .zin(e                                    /(e+1)^2)}, #g'
		function(x) {e=exp(x); .zin(e*(-e   + 1)                         /(e+1)^3)}, #g''
		function(x) {e=exp(x); .zin(e*( e^2 - 4*e    + 1)                /(e+1)^4)}, #g^(3)
		function(x) {e=exp(x); .zin(e*(-e^3 + 11*e^2 - 11*e   + 1)       /(e+1)^5)}, #g^(4)
		function(x) {e=exp(x); .zin(e*( e^4 - 26*e^3 + 66*e^2 - 26*e + 1)/(e+1)^6)}  #g^(5)
	)
)

# Utility for integration: "[return] zero if [argument is] NaN" (Inf / Inf divs)
#
# @param x Ratio of polynoms of exponentials, as in .S[[i]]
#
.zin <- function(x)
{
	x[is.nan(x)] <- 0.
	x
}
