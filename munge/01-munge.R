# see export in raw-data folder for more options.
# library(ProjectTemplate); load.project(list(munging=FALSE)) # use to debug munging file
# create variable name lists
v <- list()
v$main_scales <- c("hpmood_mean", "anxiety_mean", "depression_mean", 
                   "stress_mean", "pwi_mean", "optimism_mean", "selfesteem_mean", 
                   "extraversion_mean", "neuroticism_mean", "nwi_mean")

v$mood_scales <- c("hpmood_mean", "anxiety_mean", "depression_mean", 
                   "stress_mean", "pwi_mean",  "selfesteem_mean", 
                   "nwi_mean")
v$pwi_items <- c("s1mat", "s2hea", "s3pro", "s4int", "s5saf", 
                 "s6com", "s7sec")
v$neo_items <- paste0("neo", 1:60)
v$gls <- c("lifesat")
v$pwiall <- c(v$gls, "pwi_mean", v$pwi_items)

# length(unique(clong$user_id))
# 
# cbind(table(clong$survey))


# add wave information 
clong <- merge(clong, meta.waves[,c('wave', 'year', 'month')])
clong$wave_id <- clong$user_id * 1000 + clong$wave # fix original wave id

# merge content and contented
clong$content <- ifelse(!is.na(clong$content), clong$content, clong$contented)

# score scales
score_tests <- list()
keys <- list()

# hpmood
keys$hpmood <- data.frame(hpmood=meta.hpmood$reverse)
score_tests$hpmood <- scoreItems(keys$hpmood, items = clong[, meta.hpmood$item], impute = "none")

# dass
keys$dass <- make.keys(meta.dass$item, split(meta.dass$item, meta.dass$scale))
score_tests$dass <- scoreItems(keys$dass, items = clong[, meta.dass$item], impute = "none")

# pwi
keys$pwi <- data.frame(pwi=meta.pwi$reverse)
score_tests$pwi <- scoreItems(keys$pwi, items = clong[, meta.pwi$item], impute = "none")

# optimism
keys$optimism <- data.frame(optimism=meta.optimism$reverse)
score_tests$optimism <- scoreItems(keys$optimism, items = clong[, meta.optimism$item], impute = "none")


# selfesteem
keys$selfesteem <- data.frame(selfesteem=meta.selfesteem$reverse)
score_tests$selfesteem <- scoreItems(keys$selfesteem, items = clong[, meta.selfesteem$item], impute = "none")

# tipi
meta.tipi$item_reversed <- ifelse(meta.tipi$reverse == 1, meta.tipi$item, paste0("-", meta.tipi$item))
keys$tipi <- make.keys(meta.tipi$item, split(meta.tipi$item_reversed, meta.tipi$scale))
score_tests$tipi <- scoreItems(keys$tipi, items = clong[, meta.tipi$item], impute = "none")

# nwi
keys$nwi <- data.frame(nwi=meta.nwi$reverse)
score_tests$nwi <- scoreItems(keys$nwi, items = clong[, meta.nwi$item], impute = "none")


# Create and merge all scale scores
scale_scores <- sapply(score_tests, function(X) X$scores)
scale_scores <-  do.call(cbind, scale_scores)
scale_scores <- data.frame(scale_scores)
names(scale_scores) <- paste0(names(scale_scores), "_mean")
clong<- cbind(clong, scale_scores)

# get an sd measure for pwi
clong$pwi_sd <- apply(clong[, v$pwi_items], 1, sd)

# gender
clong$gender_text <- car::recode(clong$gender, "1='female'; 0='male'")
clong$male <- clong$gender_text == "male"


# date time
clong$year_month <- clong$year + clong$month / 12

clong$pwi_missing_count <-  apply(clong[,v$pwi_items], 1, function(X) sum(is.na(X)))
clong <- clong[clong$pwi_missing_count == 0, ]
clong <- clong[!is.na(clong$pwi_mean), ]

# create a time0 version (i.e., first wave is 0)
clong$time0 <- wavetime_by_id(clong, "user_id", "wave")



# checks
# tempids <- sample(clong$user_id, 5)
# temp <- clong[ clong$user_id %in% tempids ,c('user_id', 'wave', 'time0')]
# temp[order(temp$user_id), ]

# Create cases file

rcases <- aggregate(clong$user_id, list(user_id=clong$user_id), length)
names(rcases) <- c('user_id', 'waves_completed')

# number of waves to include in main analysis
rcases$valid_wavecount <- rcases$waves_completed >= 5
tempids <- rcases[ rcases$valid_wavecount, 'user_id']
cclong <-  clong[ clong$user_id %in% tempids, ]


# problematic cases
val_gen <- sapply(split(cclong$gender, cclong$user_id), valid_gender_sequence)
val_gen <- data.frame(user_id=as.numeric(names(val_gen)), valid_gender=val_gen)
rcases <- merge(rcases, val_gen, all.x=TRUE)

# val_age <- sapply(split(cclong$age, cclong$user_id), valid_age_sequence)
# val_age <- data.frame(user_id=as.numeric(names(val_age)), valid_age=val_age)
# rcases <- merge(rcases, val_age, all.x=TRUE)

# rcases$retain <-  (rcases$valid_age & rcases$valid_gender & rcases$valid_wavecount) %in% TRUE

# get age estimate
cclong$est_dob <- cclong$year_month - cclong$age + .5
temp <- aggregate(cclong$est_dob, list(user_id = cclong$user_id), median, na.rm=TRUE)
names(temp) <- c("user_id", "dob")
rcases <- merge(rcases, temp)



temp <- aggregate(cclong$gender_text, list(user_id=cclong$user_id), function(X) names(sort(table(X), decreasing=TRUE))[1] )
names(temp) <- c("user_id", "gender_text")
rcases <- merge(rcases, temp)

cclong$gender_text <- NULL
cclong <- merge(rcases[,c('user_id', 'waves_completed', 'gender_text', 'dob')], cclong, by="user_id")
cclong$ageraw <- cclong$age
cclong$age <- NULL
cclong$age_est <- cclong$year_month - cclong$dob

# add first wave date to rcases 
temp <- cclong[cclong$time0 == 0, c('user_id', 'year_month')]
names(temp) <- c('user_id', 'date_first_wave')
rcases <- merge(rcases, temp)

rcases$age_first_wave <- rcases$date_first_wave - rcases$dob


temp <- aggregate(pwi_mean ~ user_id, cclong, mean)
names(temp) <- c("user_id", "mean_pwi_mean")
rcases <- merge(rcases, temp)

temp <- aggregate(pwi_mean ~ user_id, cclong, sd)
names(temp) <- c("user_id", "sd_pwi_mean")
rcases <- merge(rcases, temp)


# remove cases where age is far from prediction
cclong <- cclong[!(abs(cclong$ageraw - cclong$age_est) > 3) %in% TRUE, ]

# remove invalid gender waves
cclong$invalid_gender <- (cclong$gender == 0 & cclong$gender_text == "female") %in% TRUE |  (cclong$gender == 1 & cclong$gender_text == "male") %in% TRUE
cclong <- cclong[!cclong$invalid_gender, ]


# re-remove less than 5 waves
validids <- as.numeric(names(table(cclong$user_id))[table(cclong$user_id) >= 5])
rcases <- rcases[ rcases$user_id %in% validids, ]
cclong <- cclong[ cclong$user_id %in% validids, ]

# remove neo items
cclong <- cclong[,setdiff(names(cclong), v$neo_items)] 

cclong$pwi_mean_probit <- qnorm((cclong$pwi_mean + .5) /11)
#  x <- c(0, 3, 5, 7, 10) ;  y <- qnorm((x + .5) /11);  cbind(x, y=round(y, 2))

temp <- aggregate(pwi_mean_probit ~ user_id, cclong, mean)
names(temp) <- c("user_id", "mean_pwi_mean_probit")
rcases <- merge(rcases, temp)

temp <- aggregate(pwi_mean_probit ~ user_id, cclong, sd)
names(temp) <- c("user_id", "sd_pwi_mean_probit")
rcases <- merge(rcases, temp)


# recreate a time0 version (i.e., first wave is 0)
cclong$time0 <- wavetime_by_id(cclong, "user_id", "wave")



temp <- aggregate(cclong$user_id, list(user_id=cclong$user_id), length)
names(temp) <- c('user_id', 'waves_completed2')            
rcases <- merge(rcases, temp)
rcases$waves_completed <- rcases$waves_completed2
rcases$waves_completed2 <- NULL


cclong <- merge(cclong, rcases[,c("user_id", "mean_pwi_mean", "mean_pwi_mean_probit", 
                                  "sd_pwi_mean", "sd_pwi_mean_probit")], all.x = TRUE)
cclong$residual_pwi_mean <- cclong$pwi_mean - cclong$mean_pwi_mean
cclong$stdresidual_pwi_mean <- cclong$residual_pwi_mean / cclong$sd_pwi_mean
cclong$residual_pwi_mean_probit <- cclong$pwi_mean_probit - cclong$mean_pwi_mean_probit
cclong$stdresidual_pwi_mean_probit <- cclong$residual_pwi_mean_probit / cclong$sd_pwi_mean_probit









