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
              sigma_logA0 = sigma_logA0)
}

dA_test <- test_full %>% 
  split(.$xs) %>% 
  map(~.[["dA"]])

stanPriorList <- test_smry %>% 
  split(.$xs) %>% 
  map(makeStanPriors)

sum(!(names(stanPriorList) == names(dA_test)))

standatalist <- stanPriorList %>% 
  map2(dA_test, function(x, y) {x$dA = y; x})

# 1. Flat, unbounded (above) priors on both params ------------------------


mod1 <- stan_model("src/A0_flatBoth.stan")

optims1 <- standatalist[1:15] %>% 
  map(safely(optimizing), object = mod1)

errs1 <- map_lgl(optims1, ~!is.null(.[[2]]))
sum(errs1)
vals1 <- optims1[!errs1] %>% 
  map(~.[["result"]][["par"]][1:2]) %>% 
  map(~data.frame(A0 = .[1], sigma_A = .[2])) %>% 
  bind_rows(.id = "xs")


# 2. Prior on sigsq_logA only ---------------------------------------------


# 3. Prior on both --------------------------------------------------------


# 4. Prior on both, with latent param for asymmetry -----------------------


