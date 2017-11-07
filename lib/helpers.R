# 


# Fix timezones
# Not exactly precise--I'm making some strong assumptions regarding data input errors.
fixtzs <- function(tz) {
  olsons <- c(ADT = "America/Blanc-Sablon", 
              AKDT = "America/Anchorage", 
              AKST = "America/Anchorage", 
              AST = "America/Blanc-Sablon",
              CAST = "America/Chicago", 
              CDT = "America/Chicago", 
              CETDST = "America/Chicago", 
              CST = "America/Chicago", 
              DST = "America/New_York", 
              EDT = "America/New_York",
              HDT = "America/New_York", 
              MDT = "America/Denver", 
              PDT = "America/Los_Angeles",
              PST = "America/Los_Angeles")
  out <- plyr::revalue(tz, replace = olsons, warn_missing = FALSE)
  out
}

#' Least-squares slope coefficient for y ~ x
#' @param x independent variable
#' @param y dependent variable
slope <- function(x, y)
  cor(x, y) * sd(y) / sd(x)
senslope <- function(x, y) 
  zyp.sen(y ~ x, data.frame(x = x, y = y))[["coefficients"]][[2]]

#' calculate dH/dW for various width quantiles
#' 
#' Returns a single number that is the average slope for a given range of width percentiles
#' 
#' @param H stage height or eleveation
#' @param W stream width
#' @param lwr_p Fraction of widths smaller than range considered
#' @param upr_p Fraction of widths larger than range considered
#' 
dHdW <- function(H, W, lwr_p, upr_p) {
  ords <- order(W)
  
  stopifnot(all(ords >= 0))
  H <- H[ords]
  W <- W[ords]
  
  minw <- quantile(W, lwr_p)
  maxw <- quantile(W, upr_p)
  
  inds <- W <= maxw & W >= minw
  stopifnot(sum(inds) > 1)
  
  out <- slope(x = W[inds], y = H[inds])
  out
}



# Plotting station variables ----------------------------------------------


ploth <- function(siteno, ...) {
  hswot %>% 
    filter(site_no == siteno) %>% 
    # `[[`("stage_va") %>% 
    plot(stage_va ~ datetime, ., main = station_nm[1], 
         xlab = "Date", ylab = "Stage (ft)", ...)
}

# uses xsdat from xsarea.R
plothw <- function(siteno, ...) {
  xsdat %>% 
    filter(xs == siteno) %>% 
    # `[[`("stage_va") %>% 
    plot(h_m ~ w_m, ., main = xsname[1], 
         xlab = "Width (m)", ylab = "Stage (m)", ...)
}


plotqw <- function(siteno, log = "xy", ...) {
  siteno <- as.character(siteno)
  
  hswot %>% 
    filter(site_no == siteno) %>% 
    plot(q_va ~ stream_wdth_va, ., main = station_nm[1],
         xlab = "width (feet)", ylab = "flow (CFS)", ...)
}


orthog <- function(x, y) {
  b <- y / sqrt(y %*% y)
  a <- x %*% b
  out <- x - a * b
  out
}

## Reference for hydroswot columns
fieldref <- read.csv("data/FieldDefinitions.txt", sep = "\t")
fr <- function(field)
  dplyr::filter(fieldref, grepl(field, Field))


## changepoint statistic
cpstat <- function(x, ...) {
  cp <- cpt.mean(x, method = "BinSeg", Q = 1, ...)
  out <- cp@pen.value.full
  out
}

## AHG b parameter
getb <- function(logW, logQ)
  cor(logW, logQ) * sd(logW) / sd(logQ)
