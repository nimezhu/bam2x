poisson_cdf <- function(x,lambda,lower.tail=TRUE)
{
	pp=ppois(x, lambda, lower.tail = lower.tail)
	return(pp)
}
