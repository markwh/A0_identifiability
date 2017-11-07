# prior.R
# Mark Hagemann
# 11/06/2017


# Need to get prior on 3 parameters: logA0 (median), sigma_logA, and 
# sigma_logA0, where the last one is standard deviation of median from mean.

pairs(train_smry %>% select(-xs, -A0))

plot(logA0 ~ logW_sd, train_smry)

A0lm <- lm(logA0 ~ logW_mean + logW_sd, data = train_smry)

summary(A0lm)
anova(A0lm)

# Prior is normal(-1.28197 + 1.47542 * logW_mean - 0.27367 * logW_sd, sigma = 0.5125)


# sigma_logA --------------------------------------------------------------

siglogA_lm <- lm(logA_sd ~ logW_sd + logW_mean, data = train_smry)
summary(siglogA_lm)
anova(siglogA_lm)

# prior is truncated normal(0.305138 + 1.123502 * logW_sd - 0.017514 * logW_mean, sigma = 0.2029)



# sigma_logA0 -------------------------------------------------------------

meanmed_dev <- with(train_smry, logA0_mean - logA0)
plot(meanmed_dev^2 ~ train_smry$n, log = "y") # I see no substantial dependency here.
hist(meanmed_dev)
sd(meanmed_dev)

# parameter value is 0.1520559



# Predicting using priors only --------------------------------------------

logA0_pred <- predict(A0lm, newdata = test_smry)

plot(logA0_pred ~ test_smry$logA0)
sd(logA0_pred - test_smry$logA0)
