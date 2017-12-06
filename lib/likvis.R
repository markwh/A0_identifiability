# visualization and animation functions


# Negative log likelihood -------------------------------------------------

logPost <- function(datadf, priors = NULL) {
  
  stopifnot(inherits(datadf, "data.frame"))
  datavec <- datadf[[1]]

  if (!is.null(priors)) {
    stopifnot(is.list(priors) && length(priors) == 3)
    # browser()
    muhat <- priors[[1]]
    sigmahat <- priors[[2]]
    sigmasd <- priors[[3]]
    
    logPrior <- function(logA0, sigma_logA) {
      # Change this for different applications
      
      out_mu <- dnorm(logA0, muhat, sigma_logA, log = TRUE)
      out_sigma <- dnorm(sigma_logA, sigmahat, sigmasd, log = TRUE)
      
      out_mu + out_sigma
      # end changes
    }
    
    # Prior gradient
    gradPrior <- function(logA0, sigma_logA) {
      dlogA0 <- - sigma_logA^(-2) * (logA0 - muhat)
      dsigma <- - sigma_logA^(-1) + sigma_logA^(-3) * (logA0 - muhat)^2 - 
        sigmasd^(-2) * (sigma_logA - sigmahat)
      out <- c(dlogA0, dsigma)
      out
    }
    
    if (length(data) == 0)
      return(logPrior)
  }

  
  loglikvec <- function(logA0, sigma_logA) {
    # Need to redo likelihood with logA0 instead of A0. 
    # And with sigma instead of sigmasquared.
    if (sigma_logA <= 0 || min(datavec + exp(logA0)) <= 0) return(datavec - Inf)
    
    - log(datavec + exp(logA0)) -  # This bit because of jacobian for dA
      log(sigma_logA) -
      1 / (2 * sigma_logA^2) * (log(datavec + exp(logA0)) - logA0)^2
  }
  
  out <- function(params, cumulative = FALSE) {
    logA0 <- params[1]
    sigma_logA <- params[2]
    llvec <- loglikvec(logA0, sigma_logA)
    if (cumulative) {
      llpiece <- cumsum(llvec)
    } else {
      llpiece <- sum(llvec)
    }
    
    target <- -llpiece
    
    if (!is.null(priors)) {
      target <- target - logPrior(logA0, sigma_logA)
    }
    target
  }
  
  attr(out, "gradient") <- function(params, cumulative = FALSE) {
    logA0 <- params[1]
    sigma_logA <- params[2]
    
    
    if (!is.null(prior)) {
      priorgrad <- gradPrior(logA0 = logA0, sigma_logA = sigma_logA)
    } else {
      priorgrad <- c(0, 0)
    }
    
    dlogA0vec <- - sigma_logA^(-2) * (log(datavec + exp(logA0)) - logA0) *
      (exp(logA0) / (datavec + exp(logA0)) - 1) - 
      (exp(logA0) / (datavec + exp(logA0)))
    dsigmavec <- - sigma_logA^(-1) + 
      sigma_logA^(-3) * (log(datavec + exp(logA0)) - logA0)^2
    
    out <- - c(priorgrad + c(sum(dlogA0vec), sum(dsigmavec)))
    out
  }
  
  out
}

# AR(1) version, phi is known

#' priors should be a list of length 3: muhat, sigmahat, sigmasd
logPost_ar1 <- function(datadf, priors = NULL, phi = 0.96) {

  datavec <- datadf[[1]]
  datevec <- datadf[[2]]
  
  if (!is.null(priors)) {
    stopifnot(is.list(priors) && length(priors) == 3)
    # browser()
    muhat <- priors[[1]]
    sigmahat <- priors[[2]]
    sigmasd <- priors[[3]]
    
    logPrior <- function(logA0, sigma_logA) {
      # Change this for different applications
      
      out_mu <- dnorm(logA0, muhat, sigma_logA, log = TRUE)
      out_sigma <- dnorm(sigma_logA, sigmahat, sigmasd, log = TRUE)
      
      out_mu + out_sigma
      # end changes
    }
    
    # Prior gradient
    gradPrior <- function(logA0, sigma_logA) {
      dlogA0 <- - sigma_logA^(-2) * (logA0 - muhat)
      dsigma <- - sigma_logA^(-1) + sigma_logA^(-3) * (logA0 - muhat)^2 - 
        sigmasd^(-2) * (sigma_logA - sigmahat)
      out <- c(dlogA0, dsigma)
      out
    }
    if (length(datavec) == 0)
      return(logPrior)
  }
  
  # AR(1) mean is \mu_t = \phi^{h_t} \mu_{t - 1}
  ndata <- length(datavec)
  diffdays <- diff(as.numeric(as.Date(datevec)))
  
  loglikvec <- function(logA0, sigma_logA) {
    
    if (sigma_logA <= 0 || min(datavec + exp(logA0)) <= 0) return(datavec - Inf)
    
    lagdata <- datavec[1:(ndata - 1)]
    
    devs <- log(lagdata + exp(logA0)) - logA0
    muadj <- c(0, phi^(diffdays) * devs)
    muvec <- logA0 + muadj
    sigmaW <- sigma_logA * sqrt(1 - phi^2)
    
    sigmafun <- function (h) 
      vapply(h, function(x) sum(phi^(0:(x - 1) * 2)), numeric(1)) 
    sigmavec <- c(sigma_logA, sigmaW * sqrt(sigmafun(diffdays)))
    
    - log(datavec + exp(logA0)) -  # This bit because of jacobian for dA
      log(sigmavec) -
      1 / (2 * sigmavec^2) * (log(datavec + exp(logA0)) - muvec)^2
  }
  
  out <- function(params, cumulative = FALSE) {
    logA0 <- params[1]
    sigma_logA <- params[2]
    llvec <- loglikvec(logA0, sigma_logA)
    if (cumulative) {
      llpiece <- cumsum(llvec)
    } else {
      llpiece <- sum(llvec)
    }
    target <- -llpiece
    
    if (!is.null(priors)) {
      target <- target - logPrior(logA0, sigma_logA)
    }
    target
  }
  
  attr(out, "gradient") <- function(params, cumulative = FALSE) {
    logA0 <- params[1]
    sigma_logA <- params[2]
    sigmaW <- sigma_logA * sqrt(1 - phi^2)
    
    muadj <- c(0, phi^(diffdays) * log(datavec[1:(ndata - 1)] + exp(logA0)))
    muvec <- logA0 + muadj

    sigmafun <- function (h) 
      vapply(h, function(x) sum(phi^(0:(x - 1) * 2)), numeric(1)) 
    sigmavec <- c(sigma_logA, sigmaW * sqrt(sigmafun(diffdays)))
    
    if (!is.null(prior)) {
      priorgrad <- gradPrior(logA0 = logA0, sigma_logA = sigma_logA)
    } else {
      priorgrad <- c(0, 0)
    }
    
    dmuvec <- c(1, 1 + phi^diffdays / (datavec[1:(ndata - 1)] + exp(logA0)) * exp(logA0))
    dlogA0vec <- - sigmavec^(-2) * (log(datavec + exp(logA0)) - muvec) *
      (exp(logA0) / (datavec + exp(logA0)) - dmuvec) - 
      (exp(logA0) / (datavec + exp(logA0)))
    
    dsigifun <- function(h) {
      sums <- vapply(h, function(x) sum(phi^(0:(x - 1) * 2)), numeric(1))
      sqrt((1 + phi^2) * sums)
    }
    
    dsigivec <- c(1, dsigifun(diffdays))
    dsigmavec <- - sigmavec^(-1) * dsigivec +
      sigmavec^(-3) * (log(datavec + exp(logA0)) - muvec)^2 * dsigivec
    
    out <- - c(priorgrad + c(sum(dlogA0vec), sum(dsigmavec)))
    out
  }
  
  out
}


# Likelihood with covariance matrix

logPost_mvn <- function(datavec, muhat, Sigma) {
  vecQuadForm <- function(mat, vec) {
    if (all(dim(as.matrix(mat)) == c(1, 1)))
      return(as.numeric(vec %*% mat %*% vec))
    
    nr <- nrow(mat)
    stopifnot(nr == ncol(mat) && nr == length(vec))
    
    oldbits <- vecQuadForm(mat[1:(nr - 1), 1:(nr - 1)], vec[1:(nr - 1)])
    
    newmat <- mat[-nr, nr]
    newbit <- oldbits[nr - 1] + 2 * vec[nr] * (vec[1:(nr - 1)] %*% newmat) +
      vec[nr]^2 * mat[nr, nr]
    
    out <- c(oldbits, newbit)
    out
  }
  
  
  
  loglikvec <- function(logA0, sigma_logA, distmat, phi = 0.95) {
    if (sigma_logA <= 0 || min(datavec + exp(logA0)) <= 0) return(datavec - Inf)
    
    R <- phi ^ distmat
    # Rinv <- solve(R)
    S <- R * sigma_logA^2
    Sinv <- solve(S)
    
    nn <- nrow(distmat)
    dets <- vapply(1:nn, function(x) det(S[1:x, 1:x]))
    
    devs <- log(datavec + exp(logA0)) - logA0
    
    - log(datavec + exp(logA0)) -  # This bit because of jacobian for dA
      1 / 2 * log(dets) -
      1 / 2 / sigma_logA^2 * vecQuadForm(devs)
  }
  
}

# Cumulative quadratic form using recursive matrix multiplication


# Construct gridded log-posterior for visualization -----------------------

# 2-D (static) version

llikGrid <- function(xvec, yvec, nllfun, nllData, ..., accum = FALSE,
                     logshift = 1) {
  
  nx <- length(xvec)
  ny <- length(yvec)
  
  xshift <- xvec[c(2:nx, nx)]
  yshift <- yvec[c(2:ny, ny)]
  
  plotparams <- data.frame(x = rep(xvec, ny), 
                           xmax = rep(xshift, ny),
                           y = rep(yvec, each = nx),
                           ymax = rep(yshift, each = nx))
  
  # nllfun <- nllFactory(nllData, ...)
  # 
  out <- plotparams %>% 
    mutate(logpost = map_dbl(map2(x, y, c), nllfun),
           lp_adj = log(logpost - min(logpost) + logshift))
  
  # if (accum) {
  #   out <- 
  # }
  
  out
}

llikContour <- function(griddf, realx, realy) {
  
  df_hats <- filter(griddf, logpost == min(logpost))
  
  ggplot(griddf, aes(x = x, y = y)) +
    geom_rect(aes(fill = lp_adj, xmin = x, xmax = xmax,
                  ymin = y, ymax = ymax)) +
    geom_contour(aes(z = lp_adj), color = "gray80") +
    scale_y_log10() +
    scale_fill_continuous(type = "viridis") + 
    annotation_logticks(sides = "l") +
    geom_point(data = data.frame(x = realx, y = realy),
               aes(x = x, y = y), size = 5) +
    geom_point(data = df_hats, size = 3, shape = 21, color = "gray")
}

#' requires a user-specified function giving the density as a function of parameters.
llik_denshist <- function(obs, nllfun, densfun, p) {
  # browser()  
  df_dens <- data.frame(dA = obs)
  
  estparams <- nlm(nllfun, p = p)$estimate
  
  out <- ggplot(df_dens) + 
    geom_histogram(aes(x = dA, y = ..ncount..), bins = 30) +
    # geom_density(aes(x = dA, y = ..scaled..)) +
    xlim(min(obs), max(obs)) +
    stat_function(fun = densfun, 
                  args = list(params = estparams))
  out
}


# MLE estimates, formatted as data.frame ----------------------------------

mledf <- function(factory, datadf, p, ..., cumulative = FALSE) {
  
  if (! cumulative) {
    logPost <- factory(datadf, ...)
    est <- nlm(logPost, p = p)[["estimate"]]
    return(as.data.frame(matrix(est, nrow = 1)))
  }
  
  dfslist <- 1:nrow(datadf) %>% 
    map(~1:.) %>% 
    map(~datadf[., ])
  # browser()
  
  logPosts <- map(dfslist, .f = factory, ...)
  
  mle1 <- nlm(logPosts[[1]], p = p)
  
  out <- matrix(NA_real_, nrow = length(logPosts), ncol = length(mle1$estimate))
  out[1, ] <- mle1$estimate
  
  lastest <- mle1$estimate
  for (i in 2:length(logPosts)) {
    mlei <- nlm(logPosts[[i]], p = lastest)
    out[i, ] <- mlei$estimate
    lastest <- mlei$estimate
  }
  
  as.data.frame(out)
}
