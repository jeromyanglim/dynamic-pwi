valid_age_sequence <- function(x) {
    # x is a sequence of ages ordered by wave
    # determine if valid
    
    test <- list()
    xn <- na.omit(x)
    
    # non-decreasing age
    test$nondecreasing <- all(diff(xn) >= 0)
    
    # no consecutive increase greater than 5
    test$noincreasegt5 <- all(diff(xn) <= 5)
    
    all(unlist(test))
}

valid_gender_sequence <- function(x) {
    # at least three responses with a given gender
    # no more than one other gender
    test <- list()
    xn <- na.omit(x)
    
    # at least three responses with a given gender
    max_response <- max(table(xn))
    test$atleastthree <-  max_response >= 3
    
    # no more than one response for other gender
    other_response <- length(xn) - max_response
    test$otherlt1 <- other_response <=1 

    all(unlist(test))
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