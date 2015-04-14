lag1_autocorrelation <- function(x) {
    temp <- acf(x, lag.max=1, "correlation", FALSE)['1']
    as.numeric(temp$acf) 
}

