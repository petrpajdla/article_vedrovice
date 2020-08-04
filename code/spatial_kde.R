# Project "Vedrovice"
# Script nr. 3.2
# SPATIAL ARRANGEMENT: Kernel density estimations
# author: Petr Pajdla
# KDE

library(here)
library(sf)
# library(raster)
library(tidyverse)

# data --------------------------------------------------------------------

ved_sf <- read_sf(here("data/temp/", "layout.shp")) %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")))
ved_exc <- read_sf(here("data/temp/", "window.shp"))

# kernel density est. -----------------------------------------------------

get_kde <- function(x, win = NA, variable = NA, value = NA) {
  if (any(is.na(win))) {
    win <- x
  }
  ncells <- nrow(x) * 2
  if (!is.na(variable) & !is.na(value)) {
    x <- filter(x, !!as.symbol(variable) == value)
  }
  coord <- st_coordinates(x)
  bbox <- st_bbox(win)
  bbox <- c(bbox[1], bbox[3], bbox[2], bbox[4])
  kde <- MASS::kde2d(coord[, 1], coord[, 2], lims = bbox, n = ncells)
  raster::raster(kde)
}

# plot(get_kde(ved_sf, ved_exc, "sex", "M"))
# plot(get_kde(ved_sf, ved_exc, "sex", "F"))
# plot(get_kde(ved_sf, ved_exc, "sex", "n. a."))
# plot(get_kde(ved_sf, ved_exc, "sex", "ind."))
# 
# plot(cut(get_kde(ved_sf, ved_exc, "sex", "M"), breaks = 4))
# 
# contour(get_kde(ved_sf, ved_exc, "sex", "M"))
# contour(get_kde(ved_sf, ved_exc, "sex", "F"))

# raster to tibble for computations...
# as_tibble(rasterToPoints(get_kde(ved_sf, ved_exc, "sex", "M")))

# plot kde ----------------------------------------------------------------

ved <- bind_cols(ved_sf, 
               x = st_coordinates(ved_sf)[, 1],
               y = st_coordinates(ved_sf)[, 2])

ved_bg <- ved %>% transmute(s = sex, x, y)

ved_fg <- ved %>% filter(sex != "ind.")

ggplot(data = ved_fg) +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  stat_density_2d(aes(x, y), color = "black", contour = TRUE, alpha = 0.4) +
  geom_point(data = ved_bg, aes(x, y, shape = s), 
             color = "gray", fill = "white", size = 1.4, show.legend = FALSE) +
  geom_sf(aes(shape = sex), fill = "white", size = 2) + 
  scale_shape_manual(values = c(22, 21, 24, 4)) +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br",
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m",
                              location = "br",
                              pad_y = unit(1, "cm")) +
  labs(shape = "sex") +
  facet_wrap(vars(sex)) +
  theme_void() +
  theme(strip.text = element_text(face = "italic"),
        legend.position = c(0.96, 0.8))

ggsave(here("plots", "plan_kde.pdf"), width = 16, height = 5)

# voronoi polygons --------------------------------------------------------

ved_voro <- st_voronoi(st_union(ved_sf), envelope = st_geometry(ved_exc))
ved_voro <- st_intersection(st_cast(ved_voro), ved_exc)
ved_voro <- bind_cols(st_as_sf(ved_voro), ved_sf)

ggplot(data = ved_sf) +
  geom_sf(data = ved_exc, color = "gray90", fill = NA, size = 4) +
  geom_sf(aes(shape = sex), fill = "white") + 
  geom_sf(data = ved_voro, aes(fill = sex), alpha = 0.5) +
  scale_shape_manual(values = c(22, 21, 24, 4)) +
  # ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
  #                                   location = "br", 
  #                                   pad_y = unit(2, "cm")) +
  # ggspatial::annotation_scale(plot_unit = "m", 
  #                             location = "br",
  #                             pad_y = unit(1, "cm")) +
  # theme(legend.position = c(0.9, 0.8)) +
  theme_void()
