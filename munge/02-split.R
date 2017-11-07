
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
         HUC2 = substr(HUC, 1, 2)) %>% 
  select(-HUC_retr) %>% 
  filter(!is.na(area_m2), w_m > 0, area_m2 > 0) %>% 
  group_by(xs) %>% 
  mutate(n = n(), A0 = median(area_m2), dA = area_m2 - A0) %>% 
  ungroup() %>% 
  filter(n > 30) %>% 
  arrange(xs, datetime)

sum(is.na(hswot2$HUC)) # less than 1/10 of original missing values.

# Make train and test sets
train_full <- hswot2 %>% 
  filter(as.numeric(HUC2) %% 2 == 1)
train_smry <- train_full %>% 
  group_by(xs, lat, lon, HUC, HUC2) %>% 
  summarize(logW_mean = mean(logW), logW_sd = sd(logW),
         A0 = median(area_m2), logA0 = log(A0),
         logA_sd = sd(log(area_m2)), logA0_mean = mean(log(area_m2)),
         n = n())


test_full <- hswot2 %>% 
  filter(as.numeric(HUC2) %% 2 == 0)
test_smry <- test_full %>% 
  group_by(xs, lat, lon, HUC, HUC2) %>% 
  summarize(logW_mean = mean(logW), logW_sd = sd(logW),
            A0 = median(area_m2), logA0 = log(A0),
            logA_sd = sd(log(area_m2)), logA0_mean = mean(log(area_m2)),
            n = n())
