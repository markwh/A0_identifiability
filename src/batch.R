# batch.R
# Mark Hagemann
# 11/10/2017
# For running and analyzing models on HPC


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

reg1 <- loadRegistry("src/A0_experiments/reg_flatBoth/", make.default = FALSE)

smry1 <- 1:nrow(reg1$defs) %>% 
  map(jobsmry, reg = reg1) %>% 
  bind_rows() %>% 
  split(.$variable)

A01 <- smry1$A0
sigma1 <- smry1$sigma_logA

cache("A01")
cache("sigma1")

# Flat A0 prior -----------------------------------------------------------

reg2 <- loadRegistry("src/A0_experiments/reg_flatA0/", make.default = FALSE)

smry2 <- 1:nrow(reg2$defs) %>% 
  map(jobsmry, reg = reg2) %>% 
  bind_rows() %>% 
  split(.$variable)

A02 <- smry2$A0
sigma2 <- smry2$sigma_logA

cache("A02")
cache("sigma2")


# Priors on both ----------------------------------------------------------

reg3 <- loadRegistry("src/A0_experiments/reg_noLatent/", make.default = FALSE)

smry3 <- 1:nrow(reg3$defs) %>% 
  map(jobsmry, reg = reg3) %>% 
  bind_rows() %>% 
  split(.$variable)

A03 <- smry3$A0
sigma3 <- smry3$sigma_logA

cache("A03")
cache("sigma3")

# With latent variable ----------------------------------------------------

reg4 <- loadRegistry("src/A0_experiments/reg_withLatent/", make.default = FALSE)

smry4 <- 1:nrow(reg4$defs) %>% 
  map(jobsmry, reg = reg4) %>% 
  bind_rows() %>% 
  split(.$variable)

A04 <- smry4$A0
sigma4 <- smry4$sigma_logA

cache("A04")
cache("sigma4")
