extract_parameter <- function(samples=cs, parameter="sigma", FUN = mean) {
    cols <- grep(paste0(parameter, "\\[.*\\]$"), colnames(samples), value=TRUE)
    df <- samples[,cols]
    params <- apply(df, 2, FUN)
    dat <- data.frame(id=seq(params), params)    
    names(dat) <- c("id", parameter)
    dat
}

merge_parameters <- function(samples=cs, parameters=coef_param) {
    x <- extract_parameter(samples, parameters[1])    
    for (i in seq(2, length(parameters))) {
        y <- extract_parameter(samples, parameters[i])    
        x <- merge(x, y)
    }
    x
}


