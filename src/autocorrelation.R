# autocorrelation.R
# Mark Hagemann
# 11/14/2017

# Assume logA AR1 coefficient the same as logQ AR1 coefficient.

# Get Q for a bunch of locations.


train_smry %>% glimpse()

train_full %>% 
  arrange(xs, datetime) %>% 
  glimpse() %>% 
  group_by(xs) %>% 
  mutate(difftime = c(NA, diff(as.numeric(datetime))) / 3600 / 24) %>% 
  glimpse() %>% 
  summary()
 
# Median difftime is 42 days. Will be hard to get AR1 coefficient here.
# So get Q from NWIS.

trainxs_nwis <- sprintf("%08d", train_smry$xs)
trainxs_nwis <- trainxs_nwis[!grepl("NA", trainxs_nwis)] # get rid of NAs

set.seed(86142)
randorder <- sample(1:length(trainxs_nwis), replace = FALSE)


library(dataRetrieval)

flow10 <- readNWISdv(siteNumbers = trainxs_nwis[randorder][1:10], parameterCd = "00060",
                     startDate = "2000-01-01")

flow100 <- readNWISdv(siteNumbers = trainxs_nwis[randorder][1:100], parameterCd = "00060",
                      startDate = "2000-01-01")

flow100 %>% 
  group_by(site_no) %>% 
  summarize(n = n())

flowAC1 <- flow100 %>% 
  transmute(site_no, Date, logQ = log(X_00060_00003)) %>% 
  filter(is.finite(logQ)) %>% 
  split(.$site_no) %>% 
  map(~arima(.$logQ, order = c(1, 0, 0))) 

flowAC3 <- flow100 %>% 
  transmute(site_no, Date, logQ = log(X_00060_00003)) %>% 
  filter(is.finite(logQ)) %>% 
  split(.$site_no) %>% 
  map(~arima(.$logQ, order = c(3, 0, 0))) 


fac1 <- flowAC1 %>% 
  map(~as.data.frame(c(as.list(.[["coef"]]), 
                       list(sigma = sqrt(.[["sigma2"]]))))) %>%
  bind_rows()
summary(fac1)  

fac3 <- flowAC3 %>% 
  map(~as.data.frame(c(as.list(.[["coef"]]), 
                       list(sigma = sqrt(.[["sigma2"]]))))) %>%
  bind_rows()

fooacf <- pacf(foo$logQ)

foodat <- whatNWISdata(siteNumber = trainxs_nwis[1])

# Simulation
sim3 <- arima.sim(model = list(ar = c(1.22, -0.368, 0.142)), n = 1000, sd = 0.215) + 5
plot(exp(sim3))


sim1 <- arima.sim(model = list(ar = c(0.96)), n = 1000, sd = 0.24) + 5
plot(exp(sim1))

# Cool. That looks convincing.
