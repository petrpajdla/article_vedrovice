# Project "Vedrovice"
# Script nr. 
# SPATIAL ARRANGEMENT - Percolation analysis
# author: Petr Pajdla
# Spatial organization within the cemetery

# package percopackage
# devtools::install_github("SCSchmidt/percopackage")
# percolation analysis
# uses in archaeology: Sophie C. Schmid 
# https://osf.io/7extc/
# https://github.com/SCSchmidt/percopackage

# libs --------------------------------------------------------------------

# library(percopackage)
library(igraph)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggforce)
library(ggspatial)
library(concaveman)
library(sf)
library(forcats)

# theme for ggplot graphics
theme_universe <- theme(panel.border = element_rect(colour = "black", 
                                                    fill = NA, 
                                                    size = 0.8),
                        panel.background = element_blank(),
                        line = element_blank(),
                        strip.background = element_blank(), 
                        strip.text = element_text(face = "italic"),
                        # axis.text.y = element_blank(), 
                        # axis.title.y = element_blank(), 
                        panel.spacing = unit(1, "lines"))

# data --------------------------------------------------------------------

ved_sf <- sf::read_sf(here::here("data/temp", "layout.shp")) %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")),
         age_sim = fct_relevel(age_sim, c("juv.", "ad.", "mat.", "ind.")),
         origin = fct_relevel(origin, c("local", "non-local", "ind."))) %>% 
  filter(pres != "dist.")
ved_coord <- bind_cols(id = as.character(ved_sf$id_burial), 
                       sf::st_coordinates(ved_sf)) %>% 
  rename(x = X, y = Y)
ved_exc <- sf::read_sf(here::here("data/temp", "window.shp"))

# experiments -------------------------------------------------------------

# ved <- readr::read_rds(here::here("data/temp", "vedrovice_dataset.RDS"))
# 
# ved_coord <- sf::st_coordinates(ved_sf)
# ved_df <- as.data.frame(cbind(ved_sf$id_burial, ved_coord))
# colnames(ved_df) <- c("PlcIndex", "Easting", "Northing")
# 
# percolate(ved_df, upper_radius = 20, lower_radius = 1, step_value = 2, limit = 25, radius_unit = 1)
# 
# plot(analysis_by_radius$radius, analysis_by_radius$num_clust)
# plot(analysis_by_radius$radius, analysis_by_radius$max_clust_size)
# 
# plot(ved_df)
# y <- percolate(data = x[1:600, ], limit = 1000, radius_unit = 1000, upper_radius = 1000, lower_radius = 200, step_value = 200)
# analysis_by_radius
# 
# sf::st_coordinates(ved_sf)


# percolation -------------------------------------------------------------

percolation <- function(x, lower = NULL, upper = NULL, step = NULL) {
  # stopifnot(is.matrix(x))
  m <- matrix(c(x$x, x$y), ncol = 2)
  rownames(m) <- x$id
  
  dist <- as.matrix(stats::dist(m, method = "euclidean"))
  
  # radius
  # triangle <- dist[lower.tri(dist)]
  # if (is.null(lower)) {
  #   lower <- min(triangle)
  # }
  # if (is.null(upper)) {
  #   upper <- max(triangle) / 4
  # }
  # if (is.null(step)) {
  #   step <- (upper - lower) / 10
  # }
  radii <- seq(lower, upper, by = step)
  
  res_lst <- vector("list", length(radii))
  names(res_lst) <- radii
  
  # distances
  dist_df <- reshape2::melt(dist)
  colnames(dist_df) <- c("from", "to", "value")
  dist_df <- dist_df[!(dist_df$from == dist_df$to), ]
  
  for (i in seq_along(radii)) {
    index <- dist_df$value <= radii[i]
    nodes <- dist_df[index, c("from", "to")]
    nodes <- sapply(nodes, as.character)
    g <- graph_from_edgelist(nodes)
    
    # get clusters
    clust_memb <- clusters(g)$membership
    # clust_memb <- clusters(g, mode = "weak")$membership
    memb <- matrix(ncol = 2, nrow = length(clust_memb))
    colnames(memb) <- c("id", "clust")
    memb[, "clust"]  <- unname(clust_memb)
    memb[, "id"] <- names(clust_memb)
    
    # result - list
    res_lst[[as.character(radii[i])]] <- memb
  }
  
  res_tbl <- lapply(res_lst, as_tibble) %>% 
    bind_rows(.id = "radius") %>% 
    mutate(radius = as.numeric(radius))
  
  # membership for different radii
  res <- vector("list", 3)
  names(res) <- c("stats", "membership_wide", "membership_long")
  
  res[["membership_long"]] <- res_tbl %>% 
    mutate(id = as.character(id)) %>% 
    arrange(id)
  
  # res[["membership_wide"]] <- res_tbl %>%
  #   pivot_wider(values_from = clust, 
  #               names_from = radius, 
  #               names_prefix = "radius") %>% 
  #   mutate(id = as.numeric(id)) %>% 
  #   arrange(id)
  
  # number of nodes
  nodes_num <- res_tbl %>% 
    group_by(radius) %>% 
    summarise(clust = unique(as.numeric(clust)), .groups = "keep") %>% 
    count() %>% 
    arrange(radius)
  
  # maximum cluster size and normalized cluster size etc.
  clust_stats <- res_tbl %>% 
    group_by(radius, clust) %>% 
    count() %>% 
    ungroup(clust) %>% 
    summarise(max_nodes = max(n),
              mean_nodes = mean(n),
              median_nodes = median(n), .groups = "drop") %>% 
    arrange(radius) %>% 
    mutate(max_nodes_norm = max_nodes / nrow(x))
  
  res$stats <- left_join(nodes_num, clust_stats, by = "radius")
  
  return(res)
}


# percolation on ved ------------------------------------------------------

p <- percolation(ved_coord, lower = 1, upper = 11, step = 0.1)


# percolation stats -------------------------------------------------------

p_breaks <- c(2.8, 3.9, 4.6, 5.5, 6.5, 7.8, 8.8, 10.4)

# number of clusters
# p$stats %>% 
#   ggplot(aes(x = radius, y = n)) +
#   geom_point() +
#   geom_path() + 
#   theme_minimal() +
#   labs(y = "number of clusters")

# maximum number of nodes per cluster
# p$stats %>%
#   ggplot(aes(x = radius)) +
#   geom_path(aes(y = max_nodes))

# normalized maximum number of nodes per cluster
p$stats %>% 
  mutate(label_radius = if_else(radius %in% p_breaks, 
                                radius, NA_real_)) %>% 
  ggplot(aes(x = radius, y = max_nodes_norm)) +
  geom_point(aes(x = label_radius), color = "gray") +
  geom_path() +
  geom_text(aes(x = radius - 0.5, y = max_nodes_norm + 0.06, label = label_radius)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(y = "max. nr of nodes (normalized)", x = "distance (m)") +
  theme_universe

ggsave(here::here("plots", "perc_distance.pdf"))

# mean and median nr of nodes
# p$stats %>% 
#   ggplot(aes(x = radius)) +
#   geom_path(aes(y = mean_nodes)) +
#   geom_path(aes(y = median_nodes), color = "blue")


# membership in clusters --------------------------------------------------

# membership
# init_coords <- bind_cols(id = 1:nrow(ved_coord), x = ved_coord[, "X"], y = ved_coord[, "Y"])

# wide tab solution
# g_wide <- left_join(p$membership_wide, init_coords, by = "id") %>% 
#   mutate(across(starts_with("radius"), factor))
# 
# g_wide %>% 
#   ggplot(aes(x, y)) +
#   geom_mark_hull(aes(filter = !is.na(radius5), group = radius5), fill = "blue", color = NA) +
#   geom_mark_hull(aes(filter = !is.na(radius4), group = radius4), fill = "green", color = NA) +
#   geom_mark_hull(aes(filter = !is.na(radius3), group = radius3), fill = "red", color = NA) +
#   # geom_mark_hull(aes(filter = !is.na(radius1.4), group = radius1.4), fill = "yellow") +
#   # geom_mark_hull(aes(filter = !is.na(radius1.2), group = radius1.2), fill = "green") +
#   # geom_mark_hull(aes(filter = !is.na(radius1), group = radius1), fill = "blue") +
#   geom_point(shape = 4) +
#   coord_fixed() +
#   theme_void()

# long tab solution
# g_long <- left_join(p$membership_long, init_coords, by = "id") %>% 
#   filter(radius %in% seq(1, 10, by = 0.4),
#          radius <= 9.4)

g_long <- left_join(p$membership_long, ved_coord, by = "id") %>% 
  filter(radius %in% c(1.6, p_breaks[1:7])) %>% # 1.6, 3.6,
  mutate(label_radius = paste("dist.", radius, "m"),
         label_radius = forcats::fct_reorder(label_radius, radius),
         id_burial = as.character(id))

# clust_col <- rep("gray", 19)

scale_facet <- tibble::tibble(
  label_radius = factor(c("dist. 7.8 m")),
  location = c("br"),
)

# using ggforce - wrongly displays some of the groups
# ggplot() +
#   geom_sf(data = ved_exc, fill = NA, color = "black", size = 0.1) +
#   ggforce::geom_mark_hull(data = g_long, aes(x, y, fill = clust), color = NA, expand = 0.025,
#                           show.legend = FALSE) +
#   # geom_point(data = g_long, aes(x, y, color = clust), show.legend = FALSE) +
#   scale_fill_manual(values = clust_col) +
#   scale_color_manual(values = clust_col) +
#   geom_sf(data = ved_sf, aes(shape = sex), size = 0.6) + 
#   scale_shape_manual(values = c(22, 21, 24, 4)) +
#   facet_wrap(vars(label_radius), ncol = 2) +
#   coord_sf() +
#   annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
#                          aes(location = location), 
#                          data = scale_facet,
#                          pad_y = unit(0.8, "cm")) +
#   annotation_scale(aes(location = location),
#                    data = scale_facet,
#                    plot_unit = "m") +
#   theme_void() +
#   theme(strip.text = element_text(face = "italic"),
#         panel.spacing = unit(1, "lines"), 
#         legend.position = "bottom", 
#         legend.direction = "horizontal")
# 
# ggsave(here::here("plots", "perc_plan.pdf"), width = 7, height = 14)

# using concaveman algorithm
get_concaveman <- function(pts) {
  map(pts, concaveman)
}

polygs <- ved_sf %>% 
  left_join(p$membership_long, by = c("id_burial" = "id")) %>% 
  select(id_burial, radius, clust) %>% 
  filter(radius %in% c(1.6, p_breaks[1:7])) %>% 
  group_by(radius) %>% 
  nest() %>% 
  mutate(data = map(data, group_by, clust),
         data = map(data, nest),
         data = map(data, rename, pts = data),
         data = map(data, ~map(.x$pts, concaveman))) %>% 
  unnest(cols = data) %>% 
  unnest(cols = data) %>% 
  st_as_sf() %>% 
  st_make_valid() %>% 
  st_buffer(1.6) %>% 
  mutate(label_radius = paste("dist.", radius, "m"),
         label_radius = forcats::fct_reorder(label_radius, radius))

disp_radii <- p_breaks[1:6]

ggplot() +
  geom_sf(data = ved_exc, fill = NA, color = "black", size = 0.2) +
  geom_sf(data = filter(polygs, radius %in% disp_radii), 
          fill = "gray", color = NA, alpha = 0.4) +
  geom_sf(data = ved_sf, 
          aes(shape = sex), 
          size = 0.6) +
  scale_shape_manual(values = c(22, 21, 24, 4)) +
  facet_wrap(vars(label_radius), ncol = 2) +
  coord_sf() +
  annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                         aes(location = location),
                         data = scale_facet,
                         pad_y = unit(0.8, "cm")) +
  annotation_scale(aes(location = location),
                   data = scale_facet,
                   plot_unit = "m") +
  theme_void() +
  theme(strip.text = element_text(face = "italic"),
        panel.spacing = unit(1, "lines"),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave(here::here("plots", "perc_plan.pdf"), width = 7, height = 10)

# numbers of clusters
polygs_ids <- filter(polygs, radius %in% disp_radii) %>% 
  group_by(radius) %>%
  mutate(cluster = 1:n())

ggplot() +
  geom_sf(data = ved_exc, fill = NA, color = "black", size = 0.1) +
  geom_sf(data = ved_sf, aes(shape = sex), size = 0.4, alpha = 0.2) +
  geom_sf(data = filter(polygs, radius %in% disp_radii), 
          fill = "gray", color = NA, alpha = 0.4) +
  # ggspatial::geom_spatial_label_repel(data = polygs_ids, aes(label = cluster)) +
  ggrepel::geom_text_repel(data = polygs_ids, 
                           aes(label = cluster, geometry = polygons),
                           stat = "sf_coordinates", 
                           size = 1.8, nudge_x = 1, nudge_y = 1) +
  scale_shape_manual(values = c(22, 21, 24, 4)) +
  facet_wrap(vars(label_radius), ncol = 2) +
  coord_sf() +
  annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                         aes(location = location),
                         data = scale_facet,
                         pad_y = unit(0.8, "cm")) +
  annotation_scale(aes(location = location),
                   data = scale_facet,
                   plot_unit = "m") +
  theme_void() +
  theme(strip.text = element_text(face = "italic"),
        panel.spacing = unit(1, "lines"),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave(here::here("plots", "perc_plan_ids.pdf"), width = 7, height = 10)

# cluster assignment ------------------------------------------------------
g_long %>% select(id, radius, cluster = clust) %>% 
  readr::write_csv(here::here("data/temp", "perc_clusters.csv"))
