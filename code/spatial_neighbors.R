# Project "Vedrovice"
# Script nr. ?
# NEIGHBOURS
# author: Petr Pajdla
# Mean number of neighbours with a given sex

library(here)
library(tidyverse)
library(broom)
library(sf)
library(spdep)
library(igraph)

set.seed(42)

# theme for ggplot graphics
theme_universe <- theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
                        panel.background = element_blank(),
                        line = element_blank(),
                        strip.background = element_blank(), 
                        axis.text.y = element_blank(), 
                        axis.title.y = element_blank(), 
                        panel.spacing = unit(1, "lines"))

# read data
ved <- read_rds(here("data/temp", "vedrovice_dataset.RDS"))

# data prep ====================================================================
# solving problems with negative coordinaes of polygon window
ved_layout <- ved$layout + 10

# create window ----------------------------------------------------------------
# not a rectangle, but a polygon closely bounding the burials using convex hull
ved_sf <- st_as_sf(bind_cols(as_tibble(ved_layout), ved$metadata), 
                   coords = c("layout_x", "layout_y")) %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")))

# excavation polygon
ved_exc <- matrix(c(c(0.5, 0.5, 0, 0, -1, -1, 1.6, 1.6, 7.6, 7.6, 9.1, 9.1, 
                      12.3, 12.3, 13.5, 13.5, 12.3, 12.3, 8.6, 8.6),
                    c(-1, 2.5, 2.5, 10.3, 10.3, 14.3, 14.3, 13.6, 13.6, 
                      12.8, 12.8, 13.4, 13.4, 12.4, 12.4, 4.5, 4.5, 1.5, 1.5, -1)), 
                  ncol = 2, byrow = FALSE)
ved_exc <- st_cast(st_as_sf(st_geometry(st_multipoint(ved_exc + 10))), "POLYGON")

# check neighbours =============================================================
ved_gabriel <- gabrielneigh(ved_layout)

ved_gabriel_lines <- ved_gabriel %>% 
  graph2nb(row.names = rownames(ved_layout)) %>% 
  nb2lines(coords = ved_layout) %>% 
  st_as_sf()

ggplot() +
  geom_sf(data = ved_exc, fill = "gray90", color = NA) +
  geom_sf(data = ved_gabriel_lines, linetype = 3) +
  geom_sf(data = ved_sf, shape = 21, fill = "white") +
  theme_void()

# ggsave(here("plots", "plan_vedrovice_gabriel.pdf"), scale = 2)

# get sex of neighbors
ved_g <- graph_from_edgelist(cbind(from = ved_gabriel$from, to = ved_gabriel$to), 
                             directed = FALSE)
V(ved_g)$name <- ved$metadata$id_burial

ved_sex <- as.character(ved$metadata$sex)
names(ved_sex) <- ved$metadata$id_burial

neigh_sex <- function(g, sex_vector) {
  res <- vector(mode = "list", length = length(V(g)))
  names(res) <- V(g)$name
  for (i in seq_along(V(g))) {
    x <- attr(neighbors(g, v = i), "names")
    res[[i]] <- unname(sex_vector[x])
  }
  return(res)
}

nb_sex <- neigh_sex(ved_g, ved_sex) %>% 
  map(tibble) %>% 
  map(set_names, "to_sex") %>% 
  bind_rows(.id = "id") %>% 
  mutate(from_sex = unname(ved_sex[id])) %>% 
  filter(from_sex != "ind.", to_sex != "ind.") %>% 
  group_by(id) %>%
  count(to_sex) %>% 
  mutate(from_sex = unname(ved_sex[id])) %>% 
  ungroup() %>% 
  group_by(from_sex, to_sex) %>% 
  summarise(mean = mean(n), sum = sum(n), .groups = "drop")

# randomization of sex for neighbors
randomize_neigh_sex <- function(g, n_sim, metadata) {
  n_sim <- n_sim
  n_bur <- nrow(metadata)
  prob <- metadata %>%
    select(id_burial, sex) %>%
    group_by(sex) %>%
    count() %>%
    mutate(prob = n / n_bur)
  res <- vector("list", n_sim)
  for (i in 1:n_sim) {
    rand_sex_vec <- randomizr::simple_ra(n_bur, prob_each = prob$prob, 
                                         conditions = prob$sex) %>% 
      as.character()
    names(rand_sex_vec) <- unlist(metadata[, "id_burial"])
    res[[i]] <- neigh_sex(g, rand_sex_vec) %>% 
      map(tibble) %>% 
      map(set_names, "to_sex") %>% 
      bind_rows(.id = "id") %>% 
      mutate(from_sex = unname(rand_sex_vec[id])) %>% 
      filter(from_sex != "ind.", to_sex != "ind.") %>% 
      group_by(id) %>%
      count(to_sex) %>% 
      mutate(from_sex = unname(rand_sex_vec[id])) %>% 
      ungroup() %>% 
      group_by(from_sex, to_sex) %>% 
      summarise(mean = mean(n), sum = sum(n), .groups = "drop")
    print(i)
  }
  return(res)
}

ved_rand_sex <- randomize_neigh_sex(ved_g, n_sim = 99, metadata = ved$metadata)

ved_rand_sex %>% 
  bind_rows() %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray80", color = NA, alpha = 0.6) +
  geom_vline(data = nb_sex, aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(to_sex ~ from_sex, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("mean neighbours") +
  theme_universe

ggsave(here("plots", "prob_neighbour_sex.pdf"), width = 12, height = 6)
