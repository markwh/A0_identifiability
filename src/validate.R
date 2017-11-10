# validate.R
# Mark Hagemann
# 11/10/2017

test_smry %>% dim()
sigma1 %>% dim()


val_smry <- test_smry %>% 
  cbind(sigma1)




plot(test_smry$A0, A02$mean)
abline(0, 1)
# plot(test_smry$A0, A03$mean, log = "xy")
plot(test_smry$A0, A03$mean)
abline(0, 1)
# plot(test_smry$A0, A04$mean, log = "xy")
plot(test_smry$A0, A04$mean)
abline(0, 1)

plot(test_smry$logA_sd, sigmalogA_pred)
abline(0, 1)
plot(test_smry$logA_sd, sigma2$mean)
abline(0, 1) 
plot(test_smry$logA_sd, sigma3$mean)
abline(0, 1)
plot(test_smry$logA_sd, sigma4$mean)
abline(0, 1)



resid_prior <- logA0_pred - test_smry$logA0
# resid_post2 <- log(A0) - test_smry$logA0
resid_post3 <- log(A03$mean) - test_smry$logA0
resid_post4 <- log(A04$mean) - test_smry$logA0

sum(abs(resid_post3) < abs(resid_prior))
sum(abs(resid_post3) > abs(resid_prior))
sum(abs(resid_post4) < abs(resid_prior))
sum(abs(resid_post4) > abs(resid_prior))

sd(resid_prior)
sd(resid_post3)
sd(resid_post4)


nn <- test_full %>% 
  group_by(xs) %>% 
  summarize(n = n()) %>% 
  `[[`("n")


sd(resid_prior[nn > 30])
sd(resid_post3[nn > 30])
sd(resid_post4[nn > 30])
