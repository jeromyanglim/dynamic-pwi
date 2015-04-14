group.center <- function(var,grp) {
    return(var-tapply(var,grp,mean,na.rm=T)[as.character(grp)])
}

within_person_correlation <- function(variables, id, data, ...) {
    # variables: character vector of variable names in data
    # id: characheter scalar of variable name representing grouping id variable in data
    # center all vaiables by id
    # correlate centered variables
    split_data <-   split(data[,variables], data[,id])
    csplit_data <- lapply(split_data, function(X) as.data.frame(scale(X, scale=FALSE)))
    cdata <- unsplit(csplit_data, data[,id])
    cor(cdata, ...)
}


between_person_correlation <- function(variables, id, data, na.rm=TRUE, ...) {
    gdata <- aggregate(data[,variables], list(data[,id]), mean, na.rm=na.rm)
    cor(gdata[,-1], ...)
}

between_within_descriptives <- function(variables, id, data, na.rm=TRUE) {
    gdata <- aggregate(data[,variables], list(data[,id]), mean, na.rm=na.rm)
    gdata <- gdata[,-1]
    mean_of_means <-  apply(gdata, 2, mean, na.rm=TRUE)
    sd_of_means <- apply(gdata, 2, sd, na.rm=TRUE)
    
    
    split_data <-  split(data[,variables], data[,id])
    sds <- data.frame(t(sapply(split_data, function(X) sapply(X, function(Y) sd(Y, na.rm=TRUE)))))
    mean_of_sds <- sapply(sds, mean)
    sd_of_sds <- sapply(sds, sd)
    data.frame(mean_of_means, sd_of_means, mean_of_sds, sd_of_sds)
}

reshape_simple <- function(long_data, id, variable_name, response) {
    # long_data: data.frame that contains the three variables, id, variable_name, and response
    # id: name of the id variable in long_data (e.g., "id")
    # variable_name: name of variable in long_data that represents the variable when in wide format
    # response: the name of the 
    # return: returns a data frame in wide format
    tempRawData <- long_data[, c(id, variable_name, response)]
    names(tempRawData) <- c(id, variable_name, 'response') 
    tempWideData <- reshape(tempRawData, idvar=id, direction='wide', timevar=variable_name)  
    names(tempWideData) <-  sub('response.', '', names(tempWideData))
    tempWideData
}


wavetime_by_id <- function(x, id, wave) {
    # Purpose: convert wave numbers e.g., (8, 10, 12) into time numbers (0,1,2)
    # x: data frame
    # id: name of id variable
    # wave: name of wave number of other numeric value indicating how to sort the data
    # return: vector of integeer ranked positions
    wavebyid <- split(x[,wave], x[,id])
    wavebyidseq <- lapply(wavebyid, function(X) rank(X) -1)
    unsplit(wavebyidseq, x[,id])
}
