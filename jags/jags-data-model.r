set_analysis <- function(analysis='standard') {
    if (analysis == "quick")
        settings <- list(n.chains= 4,  n.adapt= 100, burnin = 50, n.iter= 100, thin=1)
    if (analysis == "publication")
        settings <- list(n.chains= 4,  n.adapt= 1000, burnin = 2000, n.iter= 25000, thin=1)
    if (analysis == "standard")
        settings <- list(n.chains= 4,  n.adapt= 200, burnin = 1000, n.iter= 2000, thin=1)
    settings
}

get_data <- function(dv='pwi_mean', timevar='time0', subjectid='user_id', observationid="wave_id", include_lag=FALSE,
                     data=cclong, jitter_data=FALSE, standardise_dv=FALSE, agevar = "age_est", agecenter = TRUE, jitter_sd=1000) { 
    Data <- data[ ,c(observationid, subjectid, timevar, agevar, dv)]
    Data <- na.omit(Data)
    names(Data) <- c('observationid', 'subjectid', 'time', "age", 'y')
    Data$subject<- as.numeric(as.factor(Data$subjectid))
    
    
    # remove any missing data (essentially on the DV)
    Data <- Data[order(Data$time), ]
#     if (include_lag) {
#         lag1 <- function(x) {
#             c(NA, x[-length(x)])
#         }
#         
#         Data$ylag1 <- unsplit(lapply(split(Data$y, Data$subject), function(X)  lag1(X)), Data$subject)
#     }
    
    # removes 1st observation (i.e., missing lag data)
    Data <- na.omit(Data)
    
    N <- length(unique(Data$subject))
    if (standardise_dv) {
        Data$y <- as.numeric(scale(Data$y))
    }    
    # age
    print(round(mean(Data$age)))
    Data$age <- Data$age - round(mean(Data$age))
    
    # center time
    Data$time <- unsplit(lapply(split(Data$time, Data$subject), function(X) X - mean(X)), Data$subject)
    
    # j . i.e., observation number for participant starting with 1
    Data$time1 <-  unsplit(lapply(split(Data$time, Data$subjectid), rank), Data$subjectid)
    
    # number of observations for the ith participant
    ni <- sapply(sort(unique(Data$subject)), function(X) sum(Data$subject == X))
    
    # Convert data to input format required by JAGS
    Data <- Data[order(Data$subject), ]
    
    # index matrix
    Data$row_number <- seq(nrow(Data))
    index <- matrix(0, nrow = N, ncol = max(ni))
    for (i in seq(nrow(Data))) {
        index[Data[i, "subject"], Data[i, "time1"]] <- Data[i, "row_number"]
    }
    Data
    
    jagsdata <- list(subject=Data$subject, time=Data$time, time1 = Data$time1, 
                     index = index, age = Data$age, 
                     y=Data$y, N=N, ni = ni, observationid=Data$observationid)
    
    # add a little noise to DV
    # need to think about jittering lag effect
    if (jitter_data) {
        jagsdata$y <- jagsdata$y + rnorm(length(jagsdata$y), 0, sd(jagsdata$y)/jitter_sd)
      }
    
    list(Data=Data, jagsdata=jagsdata)
}

jags_model <- function (random_error = FALSE, linear_quadratic = FALSE, yhat = FALSE, transformed = FALSE, lag_effect=FALSE,
                        deviance = FALSE, age_effect = TRUE, ypred = FALSE ) {
    script <- 
        "
model {
for (i in 1:N) {
#$(YLAG)    mu_prime_y0[i] <- beta_intercept[i] 
#$(YLAG)        + theta_age_linear * (age[index[i,1]] - 1) + theta_age_quadratic * (age[index[i,1]] - 1)^2 
#$(YLAG)            + theta_age_cubic * (age[index[i,1]] - 1)^3
#$(YLAG)        #$(LINEAR_QUADRATIC)  + beta_linear[i] * (time[index[i,1]] - 1) + beta_quadratic[i] * (time[index[i,1]] - 1)^2
#$(YLAG)    #$(YHAT) mu_prime_y02[i] <- beta_intercept2[i] 
#$(YLAG)    #$(YHAT)     + theta_age_linear * (age[index[i,1]] - 1) + theta_age_quadratic * (age[index[i,1]] - 1)^2 
#$(YLAG)    #$(YHAT) + theta_age_cubic * (age[index[i,1]] - 1)^3
#$(YLAG)    #$(YHAT) #$(LINEAR_QUADRATIC)  + beta_linear2[i] * (time[index[i,1]] - 1) + beta_quadratic2[i] * (time[index[i,1]] - 1)^2

#$(YLAG)    mu_prime_y[index[i,1]] <- beta_intercept[i] 
#$(YLAG)        + theta_age_linear * (age[index[i,1]]) + theta_age_quadratic * (age[index[i,1]])^2 
#$(YLAG)        + theta_age_cubic * (age[index[i,1]])^3
#$(YLAG)        #$(LINEAR_QUADRATIC) + beta_linear[i] * time[index[i,1]] + beta_quadratic[i] * time[index[i,1]]^2
#$(YLAG)    mu_y[index[i,1]] <- mu_prime_y[index[i,1]] + mu_beta_lag * (y0[i] - mu_prime_y0[i])

#$(YLAG)    #$(YHAT) mu_prime_y2[index[i,1]] <- beta_intercept2[i] 
#$(YLAG)    #$(YHAT)    + theta_age_linear * (age[index[i,1]]) + theta_age_quadratic * (age[index[i,1]])^2 
#$(YLAG)    #$(YHAT)    + theta_age_cubic * (age[index[i,1]])^3
#$(YLAG)    #$(YHAT)    #$(LINEAR_QUADRATIC) + beta_linear2[i] * time[index[i,1]] + beta_quadratic2[i] * time[index[i,1]]^2
#$(YLAG)    #$(YHAT) mu_y2[index[i,1]] <- mu_prime_y2[index[i,1]] + mu_beta_lag * (y02[i] - mu_prime_y02[i])
#$(YLAG)    #$(YHAT) yhat_raw[index[i,1]]  ~ dnorm(mu_y2[index[i,1]], 
#$(YLAG)    #$(YHAT)    #$(RANDOM_ERROR)   tausq_y2[i])
#$(YLAG)    #$(YHAT)    #$(FIXED_ERROR) tausq_fixed_y)
    
#$(YNOLAG)  for (j in 1:ni[i]) {      
#$(YLAG)    for (j in 2:ni[i]) {
        mu_prime_y[index[i,j]] <- beta_intercept[i] 
            + theta_age_linear * (age[index[i,j]]) + theta_age_quadratic * (age[index[i,j]])^2 
            + theta_age_cubic * (age[index[i,j]])^3
            #$(LINEAR_QUADRATIC) + beta_linear[i] * time[index[i,j]] + beta_quadratic[i] * time[index[i,j]]^2
        mu_y[index[i,j]] <- mu_prime_y[index[i,j]] #$(YLAG) + mu_beta_lag * ( y[index[i,j-1]] - mu_prime_y[index[i,j-1]])
        

        #$(YHAT) mu_prime_y2[index[i,j]] <- beta_intercept2[i] 
        #$(YHAT)    + theta_age_linear * (age[index[i,j]]) + theta_age_quadratic * (age[index[i,j]])^2 
        #$(YHAT)    + theta_age_cubic * (age[index[i,j]])^3
        #$(YHAT)    #$(LINEAR_QUADRATIC) + beta_linear2[i] * time[index[i,j]] + beta_quadratic2[i] * time[index[i,j]]^2
        #$(YHAT) mu_y2[index[i,j]] <- mu_prime_y2[index[i,j]] #$(YLAG) + mu_beta_lag * ( yhat_raw[index[i,j-1]] - mu_prime_y2[index[i,j-1]])
        #$(YHAT) yhat_raw[index[i,j]]  ~ dnorm(mu_y2[index[i,j]], 
        #$(YHAT)    #$(RANDOM_ERROR)   tausq_y2[i])
        #$(YHAT)    #$(FIXED_ERROR) tausq_fixed_y)

    }
    
    for (j in 1:ni[i]) {
        #$(RANDOM_ERROR) #$(DEVIANCE) y_log[index[i,j]] <- log(1/(sigma_y[i] * 2.506628275) * exp(0-(y[index[i,j]] - mu_y[index[i,j]])^2 / (2 * sigma_y[i]^2)))
        #$(FIXED_ERROR) #$(DEVIANCE) y_log[index[i,j]] <- log(1/(sigma_fixed_y * 2.506628275) * exp(0-(y[index[i,j]] - mu_y[index[i,j]])^2 / (2 * sigma_fixed_y^2)))
        
        y[index[i,j]]  ~ dnorm(mu_y[index[i,j]], 
        #$(RANDOM_ERROR)   tausq_y[i])
        #$(FIXED_ERROR) tausq_fixed_y)

        #$(YPRED) ypred_raw[index[i,j]] ~ dnorm(mu_y[index[i, j]], 
        #$(YPRED)    #$(RANDOM_ERROR)   tausq_y[i])
        #$(YPRED)    #$(FIXED_ERROR) tausq_fixed_y)
        #$(YPRED) #$(UNTRANSFORMED)  ypred[index[i,j]]  <- ypred_raw[index[i,j]]
        #$(YPRED) #$(TRANSFORMED) ypred[index[i,j]]  <- ( phi(ypred_raw[index[i,j]]) * 11 ) - 0.5

        #$(YHAT) #$(UNTRANSFORMED)  yhat[index[i,j]]  <- yhat_raw[index[i,j]]
        #$(YHAT) #$(TRANSFORMED) yhat[index[i,j]]  <- ( phi(yhat_raw[index[i,j]]) * 11 ) - 0.5
    }
}


    

# Random coefficients
for (i in 1:N) {    
    beta[i, 1:4] ~ dmnorm(mu_beta[], prec_beta[,])    
    beta_intercept[i] <- beta[i,1]
    beta_linear[i] <- beta[i, 2]
    beta_quadratic[i] <- beta[i, 3]
    sigma_y[i] <-  exp(beta[i,4])
    tausq_y[i] <-  pow(sigma_y[i], -2)

    #$(YHAT) beta2[i, 1:4] ~ dmnorm(mu_beta[], prec_beta[,])    
    #$(YHAT) beta_intercept2[i] <- beta2[i,1]
    #$(YHAT) beta_linear2[i] <- beta2[i, 2]
    #$(YHAT) beta_quadratic2[i] <- beta2[i, 3]
    #$(YHAT) sigma_y2[i] <-  exp(beta2[i,4])
    #$(YHAT) tausq_y2[i] <-  pow(sigma_y2[i], -2)
    
#$(YLAG)    y0[i] ~ dnorm(mu_prime_y0[i],  
#$(YLAG)        #$(RANDOM_ERROR) tausq_y[i])
#$(YLAG)        #$(FIXED_ERROR) tausq_fixed_y)
#$(YLAG)    #$(YHAT) #$(RANDOM_ERROR) y02[i] ~ dnorm(beta_intercept2[i], tausq_y[i])
#$(YLAG)    #$(YHAT) #$(FIXED_ERROR) y02[i] ~ dnorm(beta_intercept2[i], tausq_fixed_y)
}

# Priors
#$(YLAG) mu_beta_lag ~ dunif(-.6,.6)

#$(TRANSFORMED) mu_beta_intercept  ~ dunif(-2,2)
#$(UNTRANSFORMED) mu_beta_intercept  ~ dunif(0,10)
mu_beta_linear <- 0 # ~ dunif(-1, 1)
mu_beta_quadratic <- 0 # ~ dunif(-1, 1)
mu_beta_sigma ~ dunif(-4, 1)    
mu_beta[1] <- mu_beta_intercept
mu_beta[2] <- mu_beta_linear
mu_beta[3] <-  mu_beta_quadratic
mu_beta[4] <-  mu_beta_sigma

#$(AGE_EFFECT) theta_age_linear ~ dnorm(0, .066666667) # tausq = (2/30) 
#$(AGE_EFFECT) theta_age_quadratic ~ dnorm(0, 202500) # tausq = (2 / 30^2 ) ^-2
#$(AGE_EFFECT) theta_age_cubic ~ dnorm(0, 182250000) # tausq = (2 / 30^3) ^-2
#$(NOAGE_EFFECT) theta_age_linear <- 0
#$(NOAGE_EFFECT) theta_age_quadratic <- 0 
#$(NOAGE_EFFECT) theta_age_cubic <- 0


for (j in 1:4){ for (k in 1:4){  Omega[j,k] <-equals(j,k)*.1 } }
prec_beta[1:4,1:4] ~ dwish(Omega[,],4)

Sigma_mu[1:4, 1:4] <- inverse(prec_beta[,])
rho_12 <- Sigma_mu[1,2] / sqrt(Sigma_mu[1,1]* Sigma_mu[2,2])
rho_13 <- Sigma_mu[1,3] / sqrt(Sigma_mu[1,1]* Sigma_mu[3,3])
rho_14 <- Sigma_mu[1,4] / sqrt(Sigma_mu[1,1]* Sigma_mu[4,4])
rho_23 <- Sigma_mu[2,3] / sqrt(Sigma_mu[2,2]* Sigma_mu[3,3])
rho_24 <- Sigma_mu[2,4] / sqrt(Sigma_mu[2,2]* Sigma_mu[4,4])
rho_34 <- Sigma_mu[3,4] / sqrt(Sigma_mu[3,3]* Sigma_mu[4,4])

sigma_beta_intercept <- sqrt(Sigma_mu[1,1])
sigma_beta_linear <- sqrt(Sigma_mu[2,2])
sigma_beta_quadratic <- sqrt(Sigma_mu[3,3])
sigma_beta_sigma <- sqrt(Sigma_mu[4,4])

#$(FIXED_ERROR) sigma_fixed_y ~ dunif(0.01, 3)

# Transformations
#$(DEVIANCE) deviance_y <- -2 * sum(y_log)
mean_sigma_y <- mean(sigma_y)
sd_sigma_y <- sd(sigma_y)


#$(FIXED_ERROR) tausq_fixed_y <- 1/(sigma_fixed_y^2)
}"

    macros <- list(
    list("#$(YHAT)",
             ifelse(yhat, "", "#")),
    list("#$(YPRED)",
         ifelse(ypred, "", "#")),
    list("#$(LINEAR_QUADRATIC)",
         ifelse(linear_quadratic, "", "#")),
    list("#$(RANDOM_ERROR)",
         ifelse(random_error, "", "#")),
    list("#$(FIXED_ERROR)",
         ifelse(!random_error, "", "#")),
    list("#$(TRANSFORMED)",
         ifelse(transformed, "", "#")),
    list("#$(UNTRANSFORMED)",
         ifelse(!transformed, "", "#")),
    list("#$(YLAG)",
         ifelse(lag_effect, "", "#")),
    list("#$(YNOLAG)",
         ifelse(!lag_effect, "", "#")),
    list("#$(AGE_EFFECT)",
         ifelse(age_effect, "", "#")),
    list("#$(NOAGE_EFFECT)",
         ifelse(!age_effect, "", "#")),
    list("#$(DEVIANCE)",
         ifelse(deviance, "", "#"))
    
    )
    
    # apply macros
    for (m in seq(macros)) {
        script <- gsub(macros[[m]][1], macros[[m]][2], script, fixed=TRUE)
    }
    
    parameters <- c("mu_beta_intercept", "sigma_beta_intercept", 
                    "theta_age_linear", "theta_age_quadratic", "theta_age_cubic")
    if (yhat) {
        parameters <- c(parameters, "yhat")
    }

    if (random_error) {
        parameters <- c(parameters, "mu_beta_sigma", "sigma_beta_sigma", "mean_sigma_y", "sd_sigma_y", "rho_14")
    } else {
        parameters <- c(parameters, "sigma_fixed_y")
    }
    if (lag_effect) {
        parameters <- c(parameters, "mu_beta_lag")
    }
    if (deviance) {
        parameters <- c(parameters, "deviance_y")
    }

    if (linear_quadratic) {
        parameters <- c(parameters, # "mu_beta_linear",  "mu_beta_quadratic", 
                        "sigma_beta_linear","sigma_beta_quadratic",
                        "rho_12", "rho_13", "rho_23")
    }

    if (linear_quadratic & random_error) {
        parameters <- c(parameters, "rho_24", "rho_34")
    }
    
    
    list(script=script, parameters = parameters)
}


