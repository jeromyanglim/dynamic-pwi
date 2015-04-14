# source code fromt the runjags package on CRAN
# It was added here because the cluster had an issue installing runjags, and this was the only function
# I needed
combine.mcmc <- function (mcmc.objects = list(), thin = 1, return.samples = NA, 
          collapse.chains = if (length(mcmc.objects) == 1) TRUE else FALSE, 
          vars = NA) 
{
    if (class(mcmc.objects) != "list") {
        if (any(class(mcmc.objects) == c("mcmc.list", "mcmc", 
                                         "runjags"))) {
            mcmc.objects <- list(mcmc.objects)
        }
        else {
            stop("Data must be provided as a list of or single mcmc object(s), or a list of or single mcmc.list(s) (for multiple chains)")
        }
    }
    mcmc.objects <- lapply(mcmc.objects, function(x) {
        if (class(x) == "runjags") 
            return(x$mcmc)
        else return(x)
    })
    returnmcmcs <- all(sapply(mcmc.objects, class) == "mcmc")
    if (returnmcmcs && collapse.chains) {
        warning("Can't collapse chains of single-chain mcmc objects")
        collapse.chains <- FALSE
    }
    mcmc.objects <- lapply(mcmc.objects, function(x) {
        if (class(x) == "mcmc") 
            return(as.mcmc.list(x))
        else return(x)
    })
    no.objects <- length(mcmc.objects)
    if (length(mcmc.objects) == 0) 
        stop("The list provided cannot be empty")
    n.chains <- integer(length = no.objects)
    n.params = (rowlengths <- vector("list", length = no.objects))
    vnames <- lapply(mcmc.objects, varnames)
    if (length(mcmc.objects) > 1) {
        allsame <- sapply(vnames, function(x) return(all(x == 
                                                             vnames[[1]])))
        if (!all(allsame)) 
            stop("Non matching variable names for supplied MCMC objects")
    }
    fornames <- mcmc.objects[[1]][[1]]
    selected <- matchvars(vars, varnames(fornames))
    if (is.null(dimnames(fornames)) || is.null(dimnames(fornames)[[1]])) 
        iterstart <- 1
    else iterstart <- as.numeric(dimnames(fornames)[[1]])[1]
    iterthin <- thin(fornames)
    if (no.objects > 1) {
        for (i in 1:no.objects) {
            if (class(mcmc.objects[[i]]) == "mcmc.list") {
                n.chains[i] <- length(mcmc.objects[[i]])
                returnlist <- TRUE
            }
            else {
                if (class(mcmc.objects[[i]]) == "mcmc") {
                    n.chains[i] <- 1
                    mcmc.objects[[i]] <- mcmc.list(mcmc.objects[[i]])
                    returnlist <- FALSE
                }
                else {
                    stop("Data must be provided as a list of mcmc objects, or a list of mcmc.lists (for multiple chains).")
                }
            }
            n.params[[i]] <- integer(length = n.chains[i])
            rowlengths[[i]] <- integer(length = n.chains[i])
            for (j in 1:n.chains[i]) {
                n.params[[i]][j] <- nvar(mcmc.objects[[i]][[j]])
                rowlengths[[i]][j] <- niter(mcmc.objects[[i]][[j]])
            }
            if (!all(rowlengths[[i]] == rowlengths[[i]][1])) 
                stop(paste("The chain lengths were not equal for object ", 
                           i, ".", sep = ""))
            rowlengths[[i]] <- rowlengths[[i]][1]
        }
        rowlengths <- unlist(rowlengths)
        paramsequal <- all(unlist(lapply(n.params, function(x) if (all(x == 
                                                                       n.params[[1]][1])) return(TRUE) else return(FALSE))))
        if (!(all(n.chains == n.chains[1]))) 
            ("There was an unequal number of chains between mcmc objects")
        if (!paramsequal) 
            stop("There was an unequal number of monitored variables (columns) between chains / mcmc objects")
        n.chains <- n.chains[1]
        n.params <- n.params[[1]][1]
        newobjects <- vector("list", length = n.chains)
        for (i in 1:n.chains) {
            newobjects[[i]] <- matrix(NA, nrow = 0, ncol = n.params, 
                                      dimnames = list(NULL, dimnames(mcmc.objects[[1]][[1]])[[2]]))
            for (j in 1:no.objects) {
                if (is.null(dim(mcmc.objects[[j]][[i]]))) 
                    dim(mcmc.objects[[j]][[i]]) <- c(niter(mcmc.objects[[j]][[i]]), 
                                                     1)
                newobjects[[i]] <- rbind(newobjects[[i]], mcmc.objects[[j]][[i]])
            }
            rowlengths <- nrow(newobjects[[i]])
            newiternames <- (((iterstart - 1)/iterthin):((((iterstart - 
                                                                1)/iterthin) + rowlengths) - 1) * iterthin) + 
                1
            dimnames(newobjects[[i]])[[1]] <- newiternames
            newobjects[[i]] <- mcmc(newobjects[[i]], start = iterstart, 
                                    thin = iterthin)
        }
        if (returnlist) 
            newobjects <- as.mcmc.list(newobjects)
        else newobjects <- newobjects[[1]]
    }
    else {
        newobjects <- mcmc.objects[[1]]
    }
    rowlengths <- niter(newobjects)
    startretsamples <- return.samples
    startthin <- thin
    if (!is.na(return.samples)) {
        if (return.samples > rowlengths) {
            thin <- 1
            if (return.samples != Inf) 
                warning("Specified return.samples was longer than the chains provided - returning shorter MCMC object length")
        }
        else {
            thin <- rowlengths/return.samples
        }
    }
    else {
        return.samples <- Inf
    }
    currentthin <- thin(newobjects)
    thin <- floor(thin) * currentthin
    endat <- (start(newobjects) + (thin * return.samples)) - 
        1
    suppressWarnings(newobjects <- window(newobjects, thin = thin, 
                                          end = endat))
    thevarnames <- dimnames(newobjects[[1]])
    newobjects <- newobjects[, selected, drop = FALSE]
    for (i in 1:length(newobjects)) {
        dimnames(newobjects[[i]]) <- list(thevarnames[[1]], thevarnames[[2]][selected])
    }
    if (is.null(dimnames(newobjects[[1]])) && is.null(dimnames(newobjects[[1]])[[1]])) {
        warning("NULL iteration names produced")
    }
    if (collapse.chains) {
        class(newobjects) <- "list"
        newobjects <- combine.mcmc(newobjects, collapse.chains = FALSE, 
                                   return.samples = startretsamples, thin = 1, vars = NA)
    }
    else {
        if (returnmcmcs) 
            newobjects <- newobjects[[1]]
    }
    return(newobjects)
}

matchvars <- function(vars, names){
    
    vars <- as.character(na.omit(vars))
    
    if(length(vars)>0){
        #	matched <- vapply(vars, function(m) return(grepl(paste("^",m,sep=""),names)), logical(length(names)))
        matched <- vapply(vars, function(m) return(grepl(m,paste("^",names,sep=""),fixed=TRUE)), logical(length(names)))
        
        exact <- vapply(vars, function(m) return(tolower(gsub("'","",gsub('"','',m,fixed=TRUE),fixed=TRUE)) == names), logical(length(names)))
        
        exactneeded <- t(matrix((grepl("'",vars,fixed=TRUE) | grepl('"',vars,fixed=TRUE)), ncol=length(names), nrow=length(vars)))
        
        selected <- which(apply( (matched & !exactneeded) | (exact & exactneeded) , 1, any))
        
        if(length(selected)==0) stop("No matches for the variable names supplied", call.=FALSE)
        
    }else{
        
        selected <- 1:length(names)
    }
    
    return(selected)	
    
}
