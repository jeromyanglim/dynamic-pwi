analysis_settings <- function(analysis='standard') {
    # settings <- analysis_settings('quick') 
    # i.e., always assign to variable called settings
    if (analysis == "ultraquick")
        settings <- list(n.chains= 2,  n.adapt= 50, burnin = 5, n.iter= 10, thin=1)
    if (analysis == "quick")
        settings <- list(n.chains= 4,  n.adapt= 100, burnin = 50, n.iter= 100, thin=1)
    if (analysis == "publication")
        settings <- list(n.chains= 4,  n.adapt= 1000, burnin = 2000, n.iter= 25000, thin=1)
    if (analysis == "standard")
        settings <- list(n.chains= 4,  n.adapt= 200, burnin = 1000, n.iter= 1000, thin=1)
    if (analysis == "medium")
        settings <- list(n.chains= 2,  n.adapt= 200, burnin = 500, n.iter= 1000, thin=1)
    settings
}

run_jags <- function(script, data, variable.names, dic.run=FALSE,
                     n.chains=settings$n.chains, n.adapt=settings$n.adapt, burnin=settings$burnin, 
                     n.iter=settings$n.iter, thin=settings$thin, dic.type='pD', showplots=TRUE, showsamples=TRUE,
                     plots.file=NULL, inits = NULL) {  
    # if script is file name, then import as file name
    if(file.exists(script))  script <-  paste( readLines(script, warn=FALSE) , collapse="\n")
    cat('jags.model: \n')
    tcs <- textConnection(script)
    mod <- jags.model(tcs, data=data,  n.chains=n.chains, n.adapt=n.adapt, inits=inits)
    close(tcs)
    cat('update: \n')
    update(mod, n.iter=burnin) # burn in
    cat('coda.samples: \n')
    samples <- coda.samples(model=mod, n.iter=n.iter, thin=thin, variable.names=variable.names)
    
    if (showplots) {
        plot(samples, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    }
    
    if (!is.null(plots.file)) {
        pdf(file = plots.file, onefile=TRUE)
        plot(samples, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
        dev.off()
    }
    
    if (showsamples) {
        print(summary(samples))
    }
    
    cat('dic.samples: \n')
    dic <- NULL
    if (dic.run) {
        dic <- dic.samples(mod, n.iter=n.iter, thin=thin, type=dic.type) # Deviance Information Criterion    
    }
    list(mod=mod, samples=samples, dic=dic, data=data, script=script, variable.names=variable.names)
}

summary_table_posterior <- function(mcmc_object, parameters=row.names(summary(mcmc_object)$statistics),
                        statistics=c('Mean', '2.5%', '97.5%'), digits=2,
                                    collapse=FALSE) {
    s <- summary(mcmc_object)
    s <- cbind(rbind(s$statistics), rbind(s$quantiles))
    snames <-  dimnames(mcmc_object[[1]])[[2]]
    row.names(s) <- snames
    if (length(setdiff(parameters, snames)) > 0) {
        stop(paste("these parameters not in mcmc_object:", setdiff(parameters, row.names(s))))
    }
    sb <- s[parameters,statistics, drop=FALSE]
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
    x <- trim(format(round(sb, digits), nsmall = digits))
    if (collapse) {
        x <- cbind(  apply(x, 1, function(X) paste0(X['Mean'], " (", X['2.5%'], ", ", X["97.5%"], ")")))
    }
    x
}

plot_posterior_density <- function(mcmc_object, 
    variable, ...) {
    x <- as.vector(sapply(mcmc_object, function(X) X[, variable]))
    plot(density(x), main="", ...) 
}

simulate_many <- function(thetas, FUN) {
    lapply(seq(nrow(thetas)), function(X)  FUN(thetas[X, ]))
}

# calculate_statistics <- function(simulations, dataset=Data) {
#     stats <- function (x) {
#         x$w <- as.logical(x$w)
#         list(cor = cor(x$y_d, x$y_a), 
#              mean_y_d = mean(x$y_d), 
#              mean_y_a = mean(x$y_a), 
#              sd_y_d = sd(x$y_d),
#              sd_y_a = sd(x$y_a),
#              y_a_gt_y_d = mean(x$y_a > x$y_d),
#              y_a_lt0.5sd_y_d = mean(x$y_a < (x$y_d - 0.5 * sd(x$y_d))),
#              y_a_gt0.5sd_y_d = mean(x$y_a > (x$y_d + 0.5 * sd(x$y_d))),
#              y_a_gt1sd_y_d = mean(x$y_a > (x$y_d + 1 * sd(x$y_d))),
#              u_diff_minus_w_diff =  mean(x$y_a[!x$w] - x$y_d[!x$w]) -
#                  mean(x$y_a[x$w] - x$y_d[x$w]) 
#         )
#     }
#     simulations_stats <- sapply(simulations, stats)
#     dataset_stats <- stats(Data)
#     list(simulations=simulations_stats, dataset=dataset_stats)
# }

summarise_statistics <- function(x) {
    # x: object returned from calculate statistics
    for_one_theta <- function(theta_i) {
        d <- x$dataset[[theta_i]]
        s <- unlist(x$simulations[theta_i, ])
        pvalue <- min(mean(d < s), 
                      mean(d > s)) * 2
        meanvalue <- mean(s)
        sdvalue <- sd(s)
        low95value <- as.numeric(quantile(s, .025))
        high95value <- as.numeric(quantile(s, .975))
        c(dataset=d, pvalue=pvalue, 
             mean=meanvalue, sd=sdvalue, 
             low95 =low95value,
             high95 = high95value)
    }
    
    thetas <-  names(x$dataset)
    sapply(thetas, for_one_theta)
}

get_dic <- function(x) {
    # x is dic object
    # returns deviance, penalty, and penalised deviance
    deviance <- sum(x$deviance)
    psum <- sum(x[[2]])
    penalised_deviance <- deviance + psum
    c(deviance=deviance, penalty=psum, penalised_deviance=penalised_deviance)
}

create_parameter_table <- function(samples, dics, parameter_order) {
    # samples: list of MCMC samples of parameters (length = to number of models)
    # dics: list of DIC samples from MCMC (length = to number of models)
    # parameter_order: vector of named parameter in desired order
    parameters <- lapply(samples, function(X) data.frame(summary_table_posterior(X, collapse=TRUE)))
    # naming variables is necessary for the subsequent merging
    for(i in seq(parameters)) {
        names(parameters[[i]]) <- names(parameters)[i]
    }

    parameter_table <- merge(parameters[[1]], parameters[[2]], by=0, all=TRUE)
    if (length(parameters) > 2) {
        for (i in seq(3, length(parameters))) {
            parameter_table <- merge(parameter_table, parameters[[i]], 
                                     by.x='Row.names', by.y='row.names', all=TRUE)
        }
    }
    names(parameter_table) <- c('parameter', names(samples))
    
    row.names(parameter_table) <- parameter_table$parameter
    
    parameter_table <- parameter_table[parameter_order, ]
    
    dic_table <- sapply(dics, function(X) unclass(get_dic(X)))
    dic_table <- cbind(parameter=row.names(dic_table), dic_table)
    
    rbind(parameter_table, dic_table)
}

line_num_cat <- function(x){
    tmp <- unlist(strsplit(x, "\n"))
    cat(paste0(seq_len(length(tmp)), ": ", tmp, collapse = "\n"), "\n")
}
