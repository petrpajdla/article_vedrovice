# Project "Vedrovice"
# Script nr. 1.2
# SPATIAL ARRANGEMENT
# author: Petr Pajdla
# Data prep for spatial analysis and basic plans

library(here)
# library(spatstat)
library(sf)
library(tidyverse)

# read data
ved <- read_rds(here("data/temp", "vedrovice_dataset.RDS"))

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
  # left_join(ei, by = "id_burial") %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")),
         # age_sim = if_else(is.na(age_sim), "ind.", as.character(age_sim)),
         age_sim = fct_relevel(age_sim, c("juv.", "ad.", "mat.", "ind.")),
         # origin = if_else(is.na(origin), "ind.", origin),
         origin = fct_relevel(origin, c("local", "non-local", "ind."))
  )

# ved_sf_buffer <- st_buffer(ved_sf, dist = 1.6)
# ved_sf_buffer <- st_simplify(ved_sf_buffer, dTolerance = 1)
# ved_sf_hull <- st_convex_hull(st_union(st_geometry(ved_sf)))

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
# window_poly <- spatstat::as.owin(st_cast(ved_sf_hull, "POLYGON"))
window_exc <- spatstat.geom::as.owin(ved_exc)

# plot
# g_exc <- ggplot() + 
#   geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
#   ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
#                                     location = "br",
#                                     pad_y = unit(2, "cm")) +
#   ggspatial::annotation_scale(plot_unit = "m",
#                               location = "br",
#                               pad_y = unit(1, "cm")) +
#   theme_void() +
#   theme(legend.position = c(0.9, 0.8))

g_exc <- ggplot() + 
  geom_sf(data = ved_exc, fill = "gray90", color = NA) +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br",
                                    pad_y = unit(0.8, "cm"),
                                    pad_x = unit(-0.4, "cm"), 
                                    height = unit(1, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m",
                              location = "br", 
                              height = unit(0.2, "cm")) +
  theme_void() +
  theme(legend.position = c(0.9, 0.8))

gg_baseplan <- g_exc +
  geom_sf(data = ved_exc, fill = NA, color = "gray80", size = 4) +
  # geom_sf(data = ved_exc, fill = "gray90", color = NA) +
  # geom_sf(data = ved_sf_hull, fill = NA, color = "gray40", linetype = 3) +
  geom_sf(data = ved_sf, aes(shape = pres), fill = "white") + 
  scale_shape_manual(values = c(15, 1)) +
  labs(shape = "preservation")

gg_baseplan

gg_ids <- gg_baseplan + 
  ggsflabel::geom_sf_text_repel(data = ved_sf, aes(label = id_burial))

# ggsave(here("plots", "plan_base.pdf"), gg_baseplan, scale = 2)
# ggsave(here("plots", "plan_base_ids.pdf"), gg_ids, scale = 2)

# function to get x and y coords for density estimation
get_coords <- function(sf, var) {
  var <- sym(var)
  res <- sf %>% 
    bind_cols(as_tibble(st_coordinates(sf))) %>% 
    filter(pres != "dist", !!var != "ind.") %>% 
    select(id_burial, !!var, X, Y) %>% 
    st_drop_geometry()
  return(res)
}

# plan local vs non-local
g_exc + 
  geom_sf(data = select(ved_sf, -origin), color = "gray") +
  stat_density2d(data = get_coords(ved_sf, "origin"), aes(X, Y),
                 color = "black") +
  geom_sf(data = filter(ved_sf, pres != "dist.", origin != "ind."), 
          shape = 21, fill = "white") + 
  facet_wrap(~origin)

ggsave(here("plots", "plan_origin.pdf"), width = 12, height = 6, units = "cm")
ggsave(here("plots", "plan_origin.eps"), width = 12, height = 6, units = "cm")


# plan age
g_exc + 
  geom_sf(data = select(ved_sf, -age_sim), color = "gray") +
  stat_density2d(data = get_coords(ved_sf, "age_sim"), aes(X, Y),
                 color = "black") +
  geom_sf(data = filter(ved_sf, pres != "dist.", age_sim != "ind."), 
          shape = 21, fill = "white") + 
  facet_wrap(~age_sim)

ggsave(here("plots", "plan_age.pdf"), width = 19, height = 6, units = "cm")
ggsave(here("plots", "plan_age.eps"), width = 19, height = 6, units = "cm")

# plan sex
g_exc +
  geom_sf(data = select(ved_sf, -sex), color = "gray") +
  stat_density2d(data = get_coords(ved_sf, "sex"), aes(X, Y),
                 color = "black") +
  geom_sf(data = filter(ved_sf, pres != "dist.", sex != "ind."),
          shape = 21, fill = "white") + 
  facet_wrap(~sex)

ggsave(here("plots", "plan_sex.pdf"), width = 19, height = 6, units = "cm")
ggsave(here("plots", "plan_sex.eps"), width = 19, height = 6, units = "cm")

# write layouts -----------------------------------------------------------

write_sf(ved_sf, here("data/temp", "layout.shp"))
write_sf(ved_exc, here("data/temp", "window.shp"))

