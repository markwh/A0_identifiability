# likvis examples and tests. Will require some modification to
# use as tests in eventual package. 
# Migrated from ad-hoc testing when writing lib/likvis.R (formerly src/likvis.R)
# May not actually be useful, but then again may be, and I needed it out of likvis.R.


# this bit from likvis.Rmd
testrow <- test_smry %>% 
  filter(n == max(n)) %>% 
  `[`(1, )
testdata <- test_full %>% 
  filter(xs == testrow$xs) %>% 
  mutate(Date = as.Date(datetime)) %>% 
  group_by(Date) %>% 
  filter(c(TRUE, rep(FALSE, n() - 1))) %>% 
  ungroup()

real_A0 <- testrow$A0
real_sigmalogA <- testrow$logA_sd

obs <- testdata$dA
obsdates <- as.Date(testdata$datetime)
priors <- stanPriorList[[as.character(testrow$xs)]]
# (end thisbit)


lpar1 <- logPost_ar1(datadf = data.frame(obs = obs, dates = obsdates),
                     muhat = priors$logA0_hat, 
                     sigmahat = priors$sigmalogA_hat, 
                     sigmasd = priors$sigmalogA_sd)


# test out contour function
xv1 <- seq(from = log(-min(obs - 1)), to = 10, length.out = 100)
yv1 <- 10^seq(from = -1.5, to = 1, length.out = 100)
gridar1 <- llikGrid(xvec = xv1, yvec = yv1, nllfun = lpar1)

llikContour(gridar1, realx = log(real_A0), realy = real_sigmalogA)


df_dens <- tibble(dA = obs[1:ndat])
df_hats <- mles[nrow(mles),]
densfun <- function(x, logmean, logsd) {
  xadj <- x + exp(logmean)
  out0 <- dlnorm(x = xadj, meanlog = logmean, sdlog = logsd)
  out <- out0 / max(out0)
  out
}

# test out denshist function
densfun <- function(x, params) {
  # browser()
  logmean <- params[1]
  logsd <- params[2]
  xadj <- x + exp(logmean)
  out0 <- dlnorm(x = xadj, meanlog = logmean, sdlog = logsd)
  out <- out0 / max(out0)
  out
}

llik_denshist(obs, logPost_ar1(datadf = data.frame(obs, obsdates), 
                               priors = list(priors$logA0_hat, 
                                             priors$sigmalogA_hat,
                                             priors$sigmalogA_sd)),
              densfun = densfun, p = c(7, 2))



# Things get messy below...

foo <- mledf(logPost_ar1, 
             datadf = data.frame(obs = obs, dates = obsdates), 
             p = c(8, 2),
             muhat = priors$logA0_hat, 
             sigmahat = priors$sigmalogA_hat, 
             sigmasd = priors$sigmalogA_sd)

datavecs <- 1:length(obs) %>% 
  map(~1:.) %>% 
  map(~obs[.])

datevecs <- 1:length(obs) %>% 
  map(~1:.) %>% 
  map(~obsdates[.])

logPosts <- map2(datavecs, datevecs, logPost_ar1, muhat = priors[["logA0_hat"]],
                 sigmahat = priors$sigmalogA_hat, sigmasd = priors$sigmalogA_sd)

attr(logPosts[[1]], "gradient")(c(7, 0.4))
attr(logPosts[[13]], "gradient")(c(7, 0.4))
attr(logPosts[[1]], "gradient")(c(8, 1))


mle1 <- nlm(logPosts[[1]], p = c(priors$logA0_hat, priors$sigmalogA_hat))
mle1

mynlm <- function(...) suppressWarnings(nlm(...))

mles <- logPosts %>%
  map(mynlm, p = c(priors$logA0_hat + 1, priors$sigmalogA_hat / 5)) %>% 
  map(~data.frame(logA0 = .$estimate[1], sigmalogA = .$estimate[2])) %>% 
  setNames(1:length(obs)) %>% 
  bind_rows(.id = "nobs") %>% 
  mutate(nobs = as.numeric(nobs))

head(mles, 10)

# further down
area1est <- obs[1] + exp(mles$logA0)

sdfun <- function(x) {
  # browser()
  if (length(x) < 3) return(NA)
  sd(x)
}
sdts <- 1:nrow(testdata) %>% 
  map_dbl(~sdfun(log(testdata$area_m2[seq(1, .)])))

plot(area1est, ylim = c(0, 1000))
abline(h = obs[1] + real_A0)

obs[1] - min(obs)

plot(mles$sigmalogA - sdts)
abline(h = 0)


