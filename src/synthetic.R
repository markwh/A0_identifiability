# Synthetic data that obey lognormal assumption, for comparison.

synthmu <- log(test_smry$A0)
synthsigma <- test_smry$logA_sd
synthn <- test_smry$n

synthdA <- pmap(.l = list(n = synthn, meanlog = synthmu, sdlog = synthsigma),
                .f = rlnorm) %>% 
  map(function(x) x - median(x))



stanPriorList <- test_smry %>% 
  split(.$xs) %>% 
  map(makeStanPriors) # from stanEst.R


synthdatalist <- stanPriorList %>% 
  map2(synthdA, function(x, y) {x$dA = y; x})


cache("synthdatalist")



# runs (on HPC) -----------------------------------------------------------

resubJobs <- function(reg) {
  jstat <- getJobStatus(reg = reg)
  needResub <- which(is.na(jstat$done))
  print(needResub)
  print(length(needResub))
  resetJobs(needResub, reg = reg)
  submitJobs(reg = reg, resources = list(walltime = 20))
}

reg1 <- makeRegistry(file.dir = "synth_flatBoth", packages = "A0est", 
                     make.default = FALSE, 
                     seed = 1586558)

batchMap(A0_estimate, A0data = synthdatalist, more.args = list(variant = "flatBoth"),
         method = "sampling", reg = reg1)
testJob(1, reg = reg1)
submitJobs(resources = list(walltime = 10), reg = reg1)


reg2 <- makeRegistry(file.dir = "synth_flatA0", packages = "A0est", 
                     make.default = FALSE, 
                     seed = 1586558)

batchMap(A0_estimate, A0data = synthdatalist, more.args = list(variant = "flatA0"),
         method = "sampling", reg = reg2)
testJob(1, reg = reg2)
submitJobs(resources = list(walltime = 10), reg = reg2)

reg3 <- makeRegistry(file.dir = "synth_noLatent", packages = "A0est", 
                     make.default = FALSE, 
                     seed = 1586558)

batchMap(A0_estimate, A0data = synthdatalist, more.args = list(variant = "noLatent"),
         method = "sampling", reg = reg3)
testJob(1, reg = reg3)
submitJobs(resources = list(walltime = 10), reg = reg3)


reg4 <- makeRegistry(file.dir = "synth_withLatent", packages = "A0est", 
                     make.default = FALSE, 
                     seed = 1586558)

batchMap(A0_estimate, A0data = synthdatalist, more.args = list(variant = "withLatent"),
         method = "sampling", reg = reg4)
testJob(1, reg = reg4)
submitJobs(resources = list(walltime = 10), reg = reg4)



# Eval --------------------------------------------------------------------

jobsmry <- function(jobid, reg) {
  print(jobid)
  allsmry <- loadResult(jobid, reg = reg) %>% 
    summary %>% 
    `[[`("summary")# %>% 
  out <- data.frame(median = allsmry[1:2, "50%"],
                    mean = allsmry[1:2, "mean"],
                    variable = c("A0", "sigma_logA"), row.names = NULL)
  out
}
# Flat priors -------------------------------------------------------------

reg1 <- loadRegistry("src/A0_experiments/synth_flatBoth/", make.default = FALSE)

smry1 <- 1:nrow(reg1$defs) %>% 
  map(jobsmry, reg = reg1) %>% 
  bind_rows() %>% 
  split(.$variable)

synth_A01 <- smry1$A0
synth_sigma1 <- smry1$sigma_logA

cache("synth_A01")
cache("synth_sigma1")

# Flat A0 prior -----------------------------------------------------------

reg2 <- loadRegistry("src/A0_experiments/synth_flatA0/", make.default = FALSE)

smry2 <- 1:nrow(reg2$defs) %>% 
  map(jobsmry, reg = reg2) %>% 
  bind_rows() %>% 
  split(.$variable)

synth_A02 <- smry2$A0
synth_sigma2 <- smry2$sigma_logA

cache("synth_A02")
cache("synth_sigma2")


# Priors on both ----------------------------------------------------------

reg3 <- loadRegistry("src/A0_experiments/synth_noLatent/", make.default = FALSE)

smry3 <- 1:nrow(reg3$defs) %>% 
  map(jobsmry, reg = reg3) %>% 
  bind_rows() %>% 
  split(.$variable)

synth_A03 <- smry3$A0
synth_sigma3 <- smry3$sigma_logA

cache("synth_A03")
cache("synth_sigma3")

# With latent variable ----------------------------------------------------

reg4 <- loadRegistry("src/A0_experiments/synth_withLatent/", make.default = FALSE)

smry4 <- 1:nrow(reg4$defs) %>% 
  map(jobsmry, reg = reg4) %>% 
  bind_rows() %>% 
  split(.$variable)

synth_A04 <- smry4$A0
synth_sigma4 <- smry4$sigma_logA

cache("synth_synth_A04")
cache("synth_sigma4")



# Adding autocorrelation --------------------------------------------------

rho <- 0.95
varcoef <- sum(rho^(1:1000 * 2))
sigsq_w <- synthsigma^2 / varcoef
sig_w <- sqrt(sigsq_w)

calcBigSig <- function(rho, dayvec) {
  mat1 <- abs(outer(dayvec, dayvec, FUN = `-`))
  bigsig <- rho^(dayvec)
  bigsig
  
}
