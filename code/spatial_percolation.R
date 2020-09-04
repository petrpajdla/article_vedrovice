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
library(ggplot2)
library(ggforce)

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

ved_sf <- sf::read_sf(here::here("data/temp", "layout.shp"))
ved_coord <- sf::st_coordinates(ved_sf)
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
  stopifnot(is.matrix(x))
  
  dist <- as.matrix(stats::dist(x, method = "euclidean"))
  
  # radius
  triangle <- dist[lower.tri(dist)]
  if (is.null(lower)) {
    lower <- min(triangle)
  }
  if (is.null(upper)) {
    upper <- max(triangle) / 4
  }
  if (is.null(step)) {
    step <- (upper - lower) / 10
  }
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
    mutate(id = as.numeric(id)) %>% 
    arrange(id)
  
  res[["membership_wide"]] <- res_tbl %>%
    pivot_wider(values_from = clust, 
                names_from = radius, 
                names_prefix = "radius") %>% 
    mutate(id = as.numeric(id)) %>% 
    arrange(id)
  
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

p <- percolation(ved_coord, lower = 1, upper = 11, step = 0.2)


# percolation stats -------------------------------------------------------

p_breaks <- c(2.8, 4, 4.6, 5.6, 6.4, 8.8, 10.2)

# number of clusters
p$stats %>% 
  ggplot(aes(x = radius, y = n)) +
  geom_point() +
  geom_path() + 
  theme_minimal() +
  labs(y = "number of clusters")

# maximum number of nodes per cluster
# p$stats %>% 
#   ggplot(aes(x = radius)) +
#   geom_path(aes(y = max_nodes))

# normalized maximum number of nodes per cluster
p$stats %>% 
  mutate(label_radius = if_else(radius %in% p_breaks, 
                                radius, NA_real_)) %>% 
  ggplot(aes(x = radius, y = max_nodes_norm)) +
  geom_path() +
  geom_point(aes(x = label_radius)) +
  geom_text(aes(x = radius - 0.5, y = max_nodes_norm + 0.06, label = label_radius)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(y = "max. nr of nodes (normalized)", x = "distance (m)") +
  theme_universe

ggsave(here::here("plots", "perc_distance.pdf"))

# mean and median nr of nodes
p$stats %>% 
  ggplot(aes(x = radius)) +
  geom_path(aes(y = mean_nodes)) +
  geom_path(aes(y = median_nodes), color = "blue")


# membership in clusters --------------------------------------------------

# membership
init_coords <- bind_cols(id = 1:nrow(ved_coord), x = ved_coord[, "X"], y = ved_coord[, "Y"])

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

g_long <- left_join(p$membership_long, init_coords, by = "id") %>% 
  filter(radius %in% c(1.6, 3.6, p_breaks)) %>% 
  mutate(label_radius = paste("dist.", radius, "m"),
         label_radius = forcats::fct_reorder(label_radius, radius))

ggplot() +
  geom_sf(data = ved_exc, fill = NA, color = "black", size = 0.1) +
  geom_mark_hull(data = g_long, aes(x, y, fill = clust), color = NA, expand = 0.02,
                 show.legend = FALSE) +
  scale_fill_manual(values = rep("gray60", 19)) +
  geom_point(data = init_coords, aes(x, y), shape = 4, color = "gray40", size = 0.8) + 
  facet_wrap(vars(label_radius)) +
  coord_sf() +
  theme_void() +
  theme(strip.text = element_text(face = "italic"),
        panel.spacing = unit(1, "lines"))

ggsave(here::here("plots", "perc_plan.pdf"), width = 9, height = 9)
