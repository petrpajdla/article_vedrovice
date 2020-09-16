# Project "Vedrovice"
# Script nr. 3.1
# SPATIAL ARRANGEMENT
# author: Petr Pajdla
# Spatial organization within the cemetery

library(here)
# library(spatstat)
library(sf)
library(tidyverse)

# read data
ved <- read_rds(here("data/temp", "vedrovice_dataset.RDS"))
origin <- read_rds(here("data/temp/local.rds")) %>% 
  select(id_burial, age_mean, origin)
ei <- read_csv(here("data/temp", "exceptionality.csv")) %>% 
  mutate(ei_clust = factor(fct),
         id_burial = as.character(burial)) %>% 
  select(-burial, -ei_cluster, -fct)

# data prep ====================================================================
scale_m <- 4.762 # scale by a factor so 1 unit is 1 meter
angle <- 40 / 360 * 100 * 1/2 # rotate by an angle of 40° so north faces upwards

ved_layout <- DescTools::Rotate(ved$layout * scale_m, 
                                theta = angle, mx = 0, my = 0)
ved_layout <- tibble(layout_x = ved_layout$x, layout_y = ved_layout$y)

# create window ----------------------------------------------------------------
# not a rectangle, but a polygon closely bounding the burials using convex hull
ved_sf <- st_as_sf(bind_cols(ved_layout, ved$metadata), 
                   coords = c("layout_x", "layout_y")) %>% 
  left_join(origin, by = "id_burial") %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind."))) %>% 
  left_join(ei, by = "id_burial")
# ved_sf_buffer <- st_buffer(ved_sf, dist = 1.6)
# ved_sf_buffer <- st_simplify(ved_sf_buffer, dTolerance = 1)
ved_sf_hull <- st_convex_hull(st_union(st_geometry(ved_sf)))

# excavation polygon
ved_exc <- matrix(c(c(0.5, 0.5, 0, 0, -1, -1, 1.6, 1.6, 7.6, 7.6, 9.1, 9.1, 
                      12.3, 12.3, 13.5, 13.5, 12.3, 12.3, 8.6, 8.6),
                    c(-1, 2.5, 2.5, 10.3, 10.3, 14.3, 14.3, 13.6, 13.6, 
                      12.8, 12.8, 13.4, 13.4, 12.4, 12.4, 4.5, 4.5, 1.5, 1.5, -1)), 
                  ncol = 2, byrow = FALSE)
ved_exc <- DescTools::Rotate(ved_exc * scale_m, theta = angle, mx = 0, my = 0)
ved_exc <- st_cast(
  st_as_sf(
    st_geometry(
      st_multipoint(
        matrix(c(ved_exc$x, ved_exc$y), ncol = 2)))), "POLYGON")

# window_rect <- owin(xrange = c(10, 22), yrange = c(10, 22.7))
window_poly <- spatstat::as.owin(st_cast(ved_sf_hull, "POLYGON"))
window_exc <- spatstat::as.owin(ved_exc)

# plot
gg_baseplan <- ggplot(data = ved_sf) +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  # geom_sf(data = ved_exc, fill = "gray90", color = NA) +
  # geom_sf(data = ved_sf_hull, fill = NA, color = "gray40", linetype = 3) +
  geom_sf(aes(shape = sex), fill = "white") + 
  scale_shape_manual(values = c(22, 21, 24, 4)) +
  theme_void() +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br", 
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m", 
                              location = "br",
                              pad_y = unit(1, "cm")) +
  theme(legend.position = c(0.9, 0.8))

gg_baseplan

gg_ids <- gg_baseplan + ggsflabel::geom_sf_text_repel(aes(label = id_burial))

ggsave(here("plots", "plan_base.pdf"), gg_baseplan, scale = 2)
ggsave(here("plots", "plan_ids.pdf"), gg_ids, scale = 2)

# plan local vs non-local
ggplot(data = filter(ved_sf, !is.na(origin))) +
  geom_sf(data = select(ved_sf, -origin), color = "gray90") +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  stat_density2d(aes(st_coordinates(filter(ved_sf, !is.na(origin)))[, 1], 
                     st_coordinates(filter(ved_sf, !is.na(origin)))[, 2]),
                 color = "black", alpha = 0.4) +
  geom_sf(shape = 21, fill = "white") + 
  facet_wrap(~origin) +
  theme_void() +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br",
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m",
                              location = "br",
                              pad_y = unit(1, "cm")) +
  theme(legend.position = c(0.9, 0.8))

ggsave(here("plots", "plan_origin.pdf"), width = 10.5, height = 5)

#plan age
ggplot(data = filter(ved_sf, !is.na(age_sim))) +
  geom_sf(data = select(ved_sf, -age_sim), color = "gray90") +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  stat_density2d(aes(st_coordinates(filter(ved_sf, !is.na(age_sim)))[, 1], 
                     st_coordinates(filter(ved_sf, !is.na(age_sim)))[, 2]),
                 color = "black", alpha = 0.4) +
  geom_sf(shape = 21, fill = "white") + 
  facet_wrap(~age_sim) +
  theme_void() +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br",
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m",
                              location = "br",
                              pad_y = unit(1, "cm")) +
  theme(legend.position = c(0.9, 0.8))

ggsave(here("plots", "plan_age.pdf"), width = 15.5, height = 5)

#plan sex
ggplot(data = filter(ved_sf, sex != "ind.")) +
  geom_sf(data = select(ved_sf, -sex), color = "gray90") +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  stat_density2d(aes(st_coordinates(filter(ved_sf, sex != "ind."))[, 1], 
                     st_coordinates(filter(ved_sf, sex != "ind."))[, 2]),
                 color = "black", alpha = 0.4) +
  geom_sf(shape = 21, fill = "white") + 
  facet_wrap(~sex) +
  theme_void() +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br",
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m",
                              location = "br",
                              pad_y = unit(1, "cm")) +
  theme(legend.position = c(0.9, 0.8))

ggsave(here("plots", "plan_sex.pdf"), width = 15.5, height = 5)

#plan ei
ggplot(data = filter(ved_sf, !is.na(ei_clust))) +
  geom_sf(data = select(ved_sf, -ei_clust), color = "gray90") +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  # stat_density2d(aes(st_coordinates(filter(ved_sf, !is.na(ei_clust)))[, 1], 
  #                    st_coordinates(filter(ved_sf, !is.na(ei_clust)))[, 2]),
  #                color = "black", alpha = 0.4) +
  geom_sf(shape = 21, fill = "white") + 
  facet_wrap(~ei_clust) +
  theme_void() +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br",
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m",
                              location = "br",
                              pad_y = unit(1, "cm")) +
  theme(legend.position = c(0.9, 0.8))

ggsave(here("plots", "plan_ei.pdf"), width = 15.5, height = 10)

# write layouts -----------------------------------------------------------

write_sf(ved_sf, here("data/temp", "layout.shp"))
write_sf(ved_exc, here("data/temp", "window.shp"))

