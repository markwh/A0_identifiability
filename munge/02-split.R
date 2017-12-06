
# To create test set, need HUC info. 

library(sf)
hucInfo <- st_read("data/USGS_Streamgages-NHD_Locations.shp") %>% 
  as.data.frame() %>% 
  transmute(xs = as.numeric(SITE_NO), HUC)

hswot1 <- hswot %>% 
  left_join(hucInfo, by = "xs")

glimpse(hswot1)
summary(as.factor(hswot1$HUC))

sum(is.na(hswot1$HUC)) / nrow(hswot1)

nohuc <- is.na(hswot1$HUC)

glimpse(hswot1[nohuc,])

#Need to get HUCs from dataRetrieval where the above joining doesn't work.
library(dataRetrieval)

xs_nohuc <- with(hswot1, xs[nohuc]) %>% 
  sprintf("%08d", .)

# Some don't follow 8-digit format and return NA's:
longxs_nohuc <- with(hswot1, xs[nohuc])[grepl("NA", xs_nohuc)] %>% 
  as.character()

xs_nohuc[grepl("NA", xs_nohuc)] <- longxs_nohuc


# Split retrievals into site strings of length 1000
retinds <- rep(1:72, length.out = length(xs_nohuc)) %>% sort()

siteinfo <- xs_nohuc %>% 
  split(retinds) %>% 
  map(safely(readNWISsite))

sum(sapply(siteinfo, function(x) !is.null(x$error))) # If zero, then I can bind them.

sitedf <- siteinfo %>% 
  map(~.[[1]]) %>% 
  markstats::bind_rows2() # had data type mismatches, so use my own function.


# Put the pieces together
sitemap <- sitedf %>% 
  transmute(xs = as.numeric(site_no), HUC_retr = huc_cd)

hswot2 <- hswot1 %>% 
  left_join(sitemap, by = "xs") %>% 
  mutate(HUC = ifelse(is.na(HUC), HUC_retr, HUC), 
         HUC2 = substr(HUC, 1, 2), Date = as.Date(datetime)) %>% 
  select(-HUC_retr) %>% 
  filter(!is.na(area_m2), w_m > 0, area_m2 > 0) %>% 
  group_by(Date, xs, xsname, lat, lon, HUC, HUC2) %>% 
  summarize_all(mean) %>% # Ensure there is only 1 observation per day (added 20171203)
  group_by(xs) %>% 
  mutate(n = n(), A0 = median(area_m2), dA = area_m2 - A0) %>% 
  ungroup() %>% 
  filter(n > 30) %>% 
  arrange(xs, datetime) %>% 
  unique() %>% 
  group_by(xs) %>% 
  mutate(n = n()) %>% # Do this again after comp-intensive unique() call
  ungroup() %>% 
  filter(n > 30) %>% 
  arrange(xs, datetime)


# Subset outliers based on studentized deleted residuals. First done in notebook20171202.Rmd.
refitRemove <- function(df, exclude = rep(FALSE, nrow(df)), thresh = 8) {
  
  lm1 <- lm(logQ ~ logA, df[!exclude, ])
  
  nas <- is.na(df$logQ | is.na(df$logA))
  
  dsrs <- delResid(lm1) / sqrt(delMSE(lm1))
  dsrbad <- abs(dsrs) > thresh
  dsrbad[is.na(dsrbad)] <- FALSE
  
  exclude[!exclude & !nas] <- dsrbad
  
  out <- df
  out$exclude <- exclude
  
  if (sum(dsrbad) > 0)
    out <- refitRemove(out, exclude = exclude, thresh = thresh)
  out
}

cleanse_A_outliers <- function(hswot2, thresh = 8) {
  lms <- hswot2 %>% 
    mutate(logA = log(area_m2), logQ = ifelse(is.infinite(logQ), NA, logQ)) %>% 
    split(.$xs) %>% 
    map(refitRemove, thresh = thresh) %>% 
    bind_rows()
  
  print(sprintf("Removed %s area outliers based on logQ ~ logA sdr", 
                sum(lms$exclude)))
  
  out <- lms %>% 
    filter(!exclude) %>% 
    select(-exclude, -logA) %>% 
    group_by(xs) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    filter(n > 30)
  out
}

hswot3 <- cleanse_A_outliers(hswot2, thresh = 8)

cache("hswot3")

# End outlier filtering.



sum(is.na(hswot3$HUC)) # less than 1/10 of original missing values.

# Make train and test sets
train_full <- hswot3 %>% 
  filter(as.numeric(HUC2) %% 2 == 1)
train_smry <- train_full %>% 
  group_by(xs, lat, lon, HUC, HUC2) %>% 
  summarize(logW_mean = mean(logW), logW_sd = sd(logW),
         A0 = median(area_m2), logA0 = log(A0),
         logA_sd = sd(log(area_m2)), logA0_mean = mean(log(area_m2)),
         n = n()) %>% 
  ungroup()


test_full <- hswot3 %>% 
  filter(as.numeric(HUC2) %% 2 == 0)
test_smry <- test_full %>% 
  group_by(xs, lat, lon, HUC, HUC2) %>% 
  summarize(logW_mean = mean(logW), logW_sd = sd(logW),
            A0 = median(area_m2), logA0 = log(A0),
            logA_sd = sd(log(area_m2)), logA0_mean = mean(log(area_m2)),
            n = n()) %>% 
  ungroup()

cache("train_full")
cache("train_smry")
cache("test_full")
cache("test_smry")
