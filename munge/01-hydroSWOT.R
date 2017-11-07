# Reading in data

hswot0 <- read.csv("data/SWOT_ADCP_Dataset.txt", sep = "\t")

# make a data frame with sensible column names and formats
hswot <- hswot0 %>% 
  mutate(q_meas_td = fixtzs(q_meas_td)) %>% 
  
  # (attempt to) deal with messy timezones
  group_by(q_meas_td) %>% 
  mutate(datetime = mdy_hms(q_meas_dt, tz = q_meas_td[1])) %>% 
  ungroup() %>% 
  
  # arrange by date per site
  group_by(site_no) %>% 
  arrange(datetime) %>% 
  ungroup() %>% 
  
  # subset, rename columns
  transmute(datetime = datetime,
            q_m3s = q_va * 0.028317,
            w_m = stream_wdth_va / 3.28084,
            logQ = log(q_m3s),
            logW = log(w_m),
            h_m = stage_va / 3.28084,
            xs = site_no, 
            xsname = station_nm, 
            area_m2 = xsec_area_va * 0.092903,
            lat = dec_lat_va,
            lon = dec_long_va,
            xs = site_no)


### Data for inspecting cross-sections

xsdat <- hswot %>%  
  filter(!is.na(area_m2),
         area_m2 > 0,
         !is.na(w_m),
         w_m > 0, 
         !is.na(h_m)) %>%
  tbl_df() %>% 
  group_by(xs) %>% 
  mutate(Ao = min(area_m2),
         Amed = median(area_m2),
         Abar = mean(area_m2),
         dA = area_m2 - Ao, 
         n = n()) %>% 
  ungroup()

# Ao as a function of W, dA statistics
aodat <- xsdat %>% 
  filter(n > 20) %>% 
  group_by(xs, xsname) %>% 
  dplyr::summarize(lwbar = mean(logW),
            lwsd = sd(logW),
            sdabar = mean(sqrt(dA)),
            sdasd = sd(sqrt(dA)),
            lAo = log(Ao[1]),
            logAmed = log(Amed[1]),
            logAbar = log(Abar[1]),
            ha25 = dHdW(h_m, W = w_m, lwr_p = 0, upr_p = 0.25),
            ha50 = dHdW(h_m, W = w_m, lwr_p = .25, upr_p = .50),
            ha75 = dHdW(h_m, W = w_m, lwr_p = .5, upr_p = .75),
            n = n()) %>% 
  ungroup()

## height anomalies, computed in reports/heightAnomalies.Rmd

cstats <- xsdat %>%
  filter(n > 20) %>% 
  group_by(xs, xsname) %>% 
  dplyr::summarize(cpstat = cpstat(h_m)) %>% 
  ungroup() %>% 
  mutate(xs = as.character(xs)) %>% 
  arrange(desc(cpstat))

# Removing anomalous heights using work in reports/heightAnomalies.Rmd
badHstas <- cstats$xs[1:25]
aodat_nobad <- aodat %>% 
  filter(!(xs %in% badHstas))

# data for computing AHG parameters
bdata <- hswot %>% 
  filter(w_m > 0, q_m3s > 0) %>% 
  tbl_df()

# data for computing sdlogQ hyperprior

varqdat <- xsdat %>% 
  filter(!(xs %in% badHstas),
         n > 50) %>% 
  group_by(xs, xsname) %>% 
  dplyr::summarize(lwbar = mean(logW),
            lwsd = sd(logW),
            lqbar = mean(logQ),
            lqsd = sd(logQ),
            # lqcv = sd(logQ) / mean(logQ),
            sdabar = mean(sqrt(dA)),
            sdasd = sd(sqrt(dA)),
            hsd = sd(h_m)) %>% 
  ungroup()
