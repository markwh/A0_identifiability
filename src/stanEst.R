# stanEst.R
# Mark Hagemann
# 11/08/2017



# 0. Prep data ------------------------------------------------------------

makeStanPriors <- function(smry) {
  
  logW_sd <- smry$logW_sd
  logW_mean <- smry$logW_mean
  
  
  sigma_logA_hat <- 0.305138 + 1.123502 * logW_sd - 0.017514 * logW_mean
  sigma_logA_sd <- 0.2029
  logA0_hat <- -1.28197 + 1.47542 * logW_mean - 0.27367 * logW_sd
  logA0_sd <- 0.5125
  nt <- smry$n
  sigma_logA0 <- 0.1520559
  
  out <- list(sigmalogA_hat = sigma_logA_hat,
              sigmalogA_sd = sigma_logA_sd, 
              logA0_hat = logA0_hat,
              logA0_sd = logA0_sd, 
              nt = nt,
              sigmalogA0 = sigma_logA0)
}

dA_test <- test_full %>% 
  split(.$xs) %>% 
  map(~.[["dA"]])

stanPriorList <- test_smry %>% 
  split(.$xs) %>% 
  map(makeStanPriors)
cache("stanPriorList")

sum(!(names(stanPriorList) == names(dA_test)))

standatalist <- stanPriorList %>% 
  map2(dA_test, function(x, y) {x$dA = y; x})

cache("standatalist")

maxest <- 958 # set to optionally do smaller set of estimates

# 1. Flat, unbounded (above) priors on both params ------------------------


mod1 <- stan_model("src/A0_flatBoth.stan")

optims1 <- standatalist[1:maxest] %>% 
  map(safely(optimizing), object = mod1)

errs1 <- map_lgl(optims1, ~!is.null(.[[2]]))
sum(errs1)
vals1 <- optims1[!errs1] %>% 
  map(~.[["result"]][["par"]][1:2]) %>% 
  map(~data.frame(A0 = .[1], sigma_A = .[2])) %>% 
  bind_rows(.id = "xs")

# For batchtools / HPC:
reg1 <- makeRegistry(file.dir = "reg_flatBoth", packages = "A0est", 
                     make.default = FALSE, 
                     seed = 1586558)

batchMap(A0_estimate, A0data = standatalist, more.args = list(variant = "flatBoth"),
         method = "sampling", reg = reg1)
testJob(1, reg = reg1)
submitJobs(resources = list(walltime = 10), reg = reg1)




# 2. Prior on sigsq_logA only ---------------------------------------------

# mod2 <- stan_model("src/A0_flatA0.stan")

# optims2 <- standatalist[1:maxest] %>% 
  # map(safely(optimizing), object = mod2)
# cache("optims2")

errs2 <- map_lgl(optims2, ~!is.null(.[[2]]))
sum(errs2)
vals2 <- optims2[!errs2] %>% 
  map(~.[["result"]][["par"]][1:2]) %>% 
  map(~data.frame(A0 = .[1], sigma_A = .[2])) %>% 
  bind_rows(.id = "xs")

# For batchtools / HPC:
reg2 <- makeRegistry(file.dir = "reg_flatA0", packages = "A0est", 
                     make.default = FALSE, 
                     seed = 1586558)

batchMap(A0_estimate, A0data = standatalist, more.args = list(variant = "flatA0"),
         method = "sampling", reg = reg2)
testJob(1, reg = reg2)
submitJobs(resources = list(walltime = 10), reg = reg2)


# 3. Prior on both --------------------------------------------------------

# mod3 <- stan_model("src/A0.stan")
# 
# optims3 <- standatalist[1:maxest] %>% 
#   map(safely(optimizing), object = mod3)
# cache("optims3")

errs3 <- map_lgl(optims3, ~!is.null(.[[2]]))
sum(errs3)
vals3 <- optims3[!errs3] %>% 
  map(~.[["result"]][["par"]][1:2]) %>% 
  map(~data.frame(A0 = .[1], sigma_A = .[2])) %>% 
  bind_rows(.id = "xs")

# For batchtools / HPC:
reg3 <- makeRegistry(file.dir = "reg_noLatent", packages = "A0est", 
                     make.default = FALSE, 
                     seed = 1586558)

batchMap(A0_estimate, A0data = standatalist, more.args = list(variant = "noLatent"),
         method = "sampling", reg = reg3)
testJob(1, reg = reg3)
submitJobs(resources = list(walltime = 10), reg = reg3)

# 4. Prior on both, with latent param for asymmetry -----------------------

mod4 <- stan_model("src/A0_withLatent.stan")

optims4 <- standatalist[1:maxest] %>%
  map(safely(optimizing), object = mod4)
cache("optims4")

errs4 <- map_lgl(optims4, ~!is.null(.[[2]]))
sum(errs4)
vals4 <- optims4[!errs4] %>% 
  map(~.[["result"]][["par"]][1:2]) %>% 
  map(~data.frame(A0 = .[1], sigma_A = .[2])) %>% 
  bind_rows(.id = "xs")

# For batchtools / HPC:
reg4 <- makeRegistry(file.dir = "reg_withLatent", packages = "A0est", 
                     make.default = FALSE, 
                     seed = 1586558)

batchMap(A0_estimate, A0data = standatalist, more.args = list(variant = "withLatent"),
         method = "sampling", reg = reg4)
testJob(1, reg = reg4)
submitJobs(resources = list(walltime = 10), reg = reg4)


# Validation --------------------------------------------------------------




plot(test_smry$A0[1:maxest], vals2$A0)
abline(0, 1)
plot(test_smry$A0[1:maxest], vals3$A0, log = "xy")
abline(0, 1)
plot(test_smry$A0[1:maxest], vals4$A0, log = "xy")
abline(0, 1)

plot(test_smry$logA_sd[1:maxest], vals2$sigma_A)
abline(0, 1) 
plot(test_smry$logA_sd[1:maxest], vals3$sigma_A)
abline(0, 1)


# What fraction of posterior A0 were closer to actual than prior was?

resid_prior <- logA0_pred - test_smry$logA0
resid_post2 <- log(vals2$A0) - test_smry$logA0
resid_post3 <- log(vals3$A0) - test_smry$logA0
resid_post4 <- log(vals4$A0) - test_smry$logA0

sum(abs(resid_post3) < abs(resid_prior))
sum(abs(resid_post3) > abs(resid_prior))
sum(abs(resid_post4) < abs(resid_prior))
sum(abs(resid_post4) > abs(resid_prior))

sd(resid_prior)
sd(resid_post3)
sd(resid_post4)
