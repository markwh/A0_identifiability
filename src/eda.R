# eda.R
# Mark Hagemann
# 11/7/2017



# Train / test maps -------------------------------------------------------

ggmap(mapdata)


ggplot(map_data("state")) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  geom_point(aes(x = lon, y = lat), data = train_smry, color = "gray30") +
  geom_point(aes(x = lon, y = lat), data = test_smry, color = "gray70") +
  coord_map("bonne", parameters = 35)
  