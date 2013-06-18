# plot the CDF plot <guifengwei@gmail.com>
# Cumulative Frequency Graph = cdf 

# cumulative distribution frequency

cdf <- function(data, add=FALSE, color='orange')
{
    # calculate the minimal and maximal value of the list
    data=data
    # the min and max of the data
    min_v <- min(data)
    max_v <- max(data)
    # divided into 100 bins
    # [ a,b)
    breaks = seq(min_v, max_v, length=100)
    data.cut = cut(data, breaks, right=FALSE)
    # compute the number of data in each bin
    data.freq =table(data.cut)

    cumfreq0=c(0, cumsum(data.freq))
    cumrelfreq0 = cumfreq0 / length(data)

    if(add==FALSE){

        plot(breaks, cumrelfreq0, cex=0.01, col=color)
        # join the points
        lines(breaks, cumrelfreq0, col=color, lwd=2)
    }
    
    if(add==TRUE){
        
        points(breaks, cumrelfreq0, cex=0.01, col=color)
        
        lines(breaks, cumrelfreq0, col=color, lwd=2)    
    }
}
