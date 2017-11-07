# xsarea.R
# Mark Hagemann
# 10/31/2016
# cross-sectional area computation

# Area as a function of W
plot(area_m2 ~ logW, sample_n(xsdat, 10000), log = "y")

# Distribution of Ao
aos <- unique(xsdat$Ao)
hist(log10(aos))


pairs(aodat[-1:-2])
aodat %>% 
  filter(n > 50) %>% 
  select(lwbar : ha75, n) %>% 
  pairs()

# model lAo ~ lwbar + lwsd
alm1 <- lm(lAo ~ lwbar + lwsd, data = aodat)
summary(alm1)
AIC(alm1)
alm2 <- lm(lAo ~ lwbar, data = aodat)
summary(alm2)
AIC(alm2)

## Save model 

alm <- alm1
alm$data <- aodat
cache("alm")

## Save plot
png(filename = "Ao_fit.png", width = 600, height = 300)
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
visreg::visreg(alm1)
par(op)
dev.off()

# Look at a single xs data

# xs1 <- sample(xsdat$xs, 1)
xs1 <- filter(aodat, n > 50)$xs[which.min(filter(aodat, n > 50)$ha75)]
ad1 <- filter(xsdat, xs == xs1)
with(ad1, dHdW(h_m, w_m, 0, 1))
plot(h_m ~ w_m, ad1)
ad1$h_m


# incorporating work done in reports/heightAnomalies.Rmd

summary(aodat_nobad)
pairs(aodat_nobad[-1:-2])
head(arrange(aodat_nobad, ha25))
aodat_nobad[which.min(aodat_nobad$ha25),]
plothw(2237700)

# doesn't look like dh/dw tells us anything.


# Checking distribution of logA -- does median equal mean? ----------------

foo <- xsdat %>% glimpse() %>% 
  group_by(xs) %>% 
  mutate(n = n(), logA = log(area_m2)) %>% 
  filter(n > 20) %>% glimpse() %>% 
  group_by(xs) %>% 
  summarize(min = min(logA), median = median(logA), 
            mean = mean(logA), n = n())
  

# png(filename = "graphs/area_mean_median.png", width = 4, height = 3, 
#     units = "in", res = 100)
plot(mean ~ median, foo, xlab = "median log-cross-sectional area",
     ylab = "mean log-cross-sectional area")
abline(0, 1)
# dev.off()

hist(foo$n)
quantile(foo$n, 0:10 / 10)


foo %>% 
  filter(n > 50) %>%
  ggplot(aes(x = (median - mean) / mean)) +
  geom_histogram() +
  xlim(-.25, 0.25) + 
  theme_bw()
ggsave("graphs/area_relative_difference_mean_median.png", 
       width = 4, height = 3)



# Models for A0 as median -------------------------------------------------


# model lAo ~ lwbar + lwsd
alm3 <- lm(logAmed ~ lwbar + lwsd, data = aodat)
summary(alm3)
AIC(alm3)
BIC(alm3)
alm4 <- lm(logAmed ~ lwbar, data = aodat)
summary(alm4)
AIC(alm4)
BIC(alm4)

amedlm <- alm3
cache("amedlm")
