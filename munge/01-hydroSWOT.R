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

cache("hswot")