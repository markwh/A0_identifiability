# A0 likelihood


makeGrad <- function(dA) {
  gr <- function(pars) {
    A0 <- pars[1]; sigsq_logA <- pars[2]
    aptl <- - sum(-1 / (dA + A0) - 1 / sigsq_logA * (log(dA + A0) - log(A0)) *
                    (1 / (dA + A0) - 1 / A0))
    sigsqptl <- - sum(-1 / (2 * sigsq_logA) + 1 / (2 * sigsq_logA^2) * 
                        (log(dA + A0) - log(A0))^2)
    
    out <- c(aptl, sigsqptl)
    out
  }
  gr
}


makeNLL <- function(dA, trace = FALSE, gr_att = TRUE) {
  
  gradfun <- makeGrad(dA)
  
  out <- function(pars) {
    A0 <- pars[1]
    sigsq_logA <- pars[2]
    if (A0 <= 0 | sigsq_logA <= 0 | any(dA < -A0)) {
      nll <- Inf
    } else {
      logquan <- suppressWarnings(log(dA + A0))
      loglik <- sum(-logquan - 1/2 * log(sigsq_logA) - 
                      1 / (2 * sigsq_logA) * (logquan - log(A0))^2)
      nll <- -loglik
      
    }
    if (trace) {
      cat(A0, " ", sigsq_logA, " ", nll, "\n")
    }
    
    if (gr_att) {
      # browser()
      attr(nll, "gradient") <- gradfun(c(A0, sigsq_logA))
    }
    
    nll
  }
}

# Actual maximum likelihood optimization
dA_test <- test_full %>% 
  split(.$xs) %>% 
  map(~.[["dA"]])

plist <- dA_test %>% 
  map(~c(-min(.) + 1, 0.1))

dA_ests <- dA_test %>% 
  map(makeNLL) %>% 
  map2(plist, function(...) suppressWarnings(nlm(...))) %>% 
  map(~data.frame(A0 = .[["estimate"]][1], sigsq_logA = .[["estimate"]][2])) %>% 
  bind_rows(.id = "xs") %>% 
  mutate(xs = as.numeric(xs))


# Compare to the real thing.
mle_compare <- test_smry %>% 
  left_join(rename(dA_ests, A0_est = A0, sigsq_logA_est = sigsq_logA), 
            by = "xs")
# plot(A0 ~ A0_est, mle_compare, log = "xy")
plot(A0 ~ A0_est, mle_compare)
abline(0, 1)

plot(logA_sd ~ sqrt(sigsq_logA_est), mle_compare)
# plot(logA_sd ~ sqrt(sigsq_logA_est), mle_compare, log = "xy")
abline(0, 1)


# That's a wrap on this section! On to stan!
