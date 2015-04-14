# http://ryouready.wordpress.com/2010/01/11/progress-bars-in-r-part-ii-a-wrapper-for-apply-functions/
sapply_pb <- function(X, FUN, ...)
{
    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
    
    wrapper <- function(...){
        curVal <- get("counter", envir = env)
        assign("counter", curVal +1 ,envir=env)
        setTxtProgressBar(get("pb", envir=env), curVal +1)
        FUN(...)
    }
    res <- sapply(X, wrapper, ...)
    close(pb)
    res
}



max_abs_deviation <- function (x) {
    max(abs(x - mean(x)))
}

linear_change <- function(x) {
    Data <- data.frame(time=seq(x), y=x)
    coef(lm(y~time, Data))['time']
}


posterior_stats <- function (x) {
    # x: single dataset
    x <- x[order(x$time), ] # lag effects assume data is orded by time
    list(
        mean_of_means = mean(aggregate(y~subject, x, mean)[,2]),
        sd_of_means = sd(aggregate(y~subject, x, mean)[,2]),
        skew_of_means = skew(aggregate(y~subject, x, mean)[,2]),        
        mean_of_sds = mean(aggregate(y~subject, x, sd)[,2]),
        sd_of_sds =sd(aggregate(y~subject, x, sd)[,2]),
        mean_of_lag1autocorrelations = mean(sapply(split(x$y, x$subject), function(X) lag1_autocorrelation(X))),
        mean_of_skews = mean(aggregate(y~subject, x, skew)[,2]),
        mean_of_outliers = mean(aggregate(y~subject, x, max_abs_deviation)[,2]), 
        sd_of_outliers = sd(aggregate(y~subject, x, max_abs_deviation)[,2]), 
        sd_linearchange = sd(aggregate(y~subject, x, linear_change)[,2]),
        cor_mean_sigma = cor(aggregate(y ~ subject, x, mean)[,2],  aggregate(y ~ subject, x, sd)[,2])
    )
}

calculate_statistics <- function(simulations, dataset=Data) {
    # TODO: update to ensure that it is robust to sort order of data.frame
    
    # debugonce(stats)
    cat("\nCalculate statistics on each posterior simulated dataset\n")
    simulations_stats <- sapply_pb(simulations, posterior_stats)
    dataset_stats <- posterior_stats(dataset)
    list(simulations=simulations_stats, dataset=dataset_stats)
}



posterior_predictive_check_pwi <- function(jagsdata=data$jagsdata, rawdata=data$Data, jagsmodel="jags/m1", transformed = FALSE, inits=NULL) {
    # 1. get posterior samples
    # have a good burnin: 3000
    # have many iterations: 10000 
    # have thinning: 10
    
    pc_model <- run_jags(jagsmodel, jagsdata,  'yhat', dic.run=FALSE, showplots=FALSE, showsamples=FALSE, inits=inits)
    
    
    # 2. create simulated datasets
    cat("Create simulated datasets\n")
    pc_samples <- do.call(rbind, pc_model$samples)
    simulations <- vector("list", length = nrow(pc_samples))
    for (i in seq(nrow(pc_samples)) ) {
        simulations[[i]] <- data.frame(subject=jagsdata$subject,  time=jagsdata$time, y=pc_samples[i,])
    }
    
    if(transformed) rawdata$y <-  (pnorm(rawdata$y) * 11) - 0.5
    
    # 3. calculate statistics on simulated data and data
    # debugonce(calculate_statistics)
    x <- calculate_statistics(simulations, rawdata)
    summaryx <- summarise_statistics(x)
    list(summary_statistics=summaryx, raw_statistics=x)
}

