# A0 likelihood

ao_loglik <- function(dA_vec, A0, sigsq_logA) {
  logquan <- log(dA_vec + A0)
  loglik <- sum(-logquan - 1/2 * log(sigsq_logA) - 
                  1 / (2 * sigsq_logA) * (logquan - log(A0))^2)
  loglik
}

# try using mle function in stats4.
makeNLL <- function(dA) {
  out <- function(A0, sigsq_logA) {
    logquan <- log(dA + A0)
    loglik <- sum(-logquan - 1/2 * log(sigsq_logA) - 
                    1 / (2 * sigsq_logA) * (logquan - log(A0))^2)
    nll <- -loglik
    nll
  }
}

dA_test <- test_full %>% 
  split(.$xs) %>% 
  map(~.[["dA"]])

dA_mle <- dA_test %>% 
  map(makeNLL) %>% 
  map(mle, start = list(A0 = exp(9.1)x))

test1 <- mle(makeNLL(dA = test_full$dA), start = list(A0 = 500, sigsq_logA = 0.005))

test2 <- mle(makeNLL(dA = testao$dA_med), start = list(A0 = 500), 
             fixed = list(sigsq_logA = 0.005))

