# Background art was taken under a CC-BY license from:
# https://www.vecteezy.com/vector-art/66620012-scenic-mountain-trail-landscape-hiking-path-through-appalachian-mountains
library(hexSticker)

# remove white borders to set to transparency
sticker(
  "analysis/background_track.png",
  package="homeranger",
  p_size = 20,
  s_x=1,
  s_y=0.85,
  s_width = 0.8,
  p_y = 0.6,
  filename="logo.png",
  h_color = "#253338",
  h_size = 2.1,
  u_family = "Aller_Rg",
  white_around_sticker = TRUE
)

