# Project "Vedrovice"
# Script nr. 3.2
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

# read data
# ved <- read_rds(here("data/temp", "vedrovice_dataset.RDS"))

# data prep ====================================================================
ved_sf <- read_sf(here("data/temp", "layout.shp")) %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")))
ved_layout <- st_coordinates(ved_sf)
ved_exc <- read_sf(here("data/temp", "window.shp"))

# check neighbours =============================================================
ved_gabriel <- gabrielneigh(ved_layout)

ved_gabriel_lines <- ved_gabriel %>% 
  graph2nb(row.names = rownames(ved_sf)) %>% 
  nb2lines(coords = ved_layout) %>% 
  st_as_sf()

ggplot() +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  geom_sf(data = ved_gabriel_lines, linetype = 3) +
  geom_sf(data = ved_sf, aes(shape = sex), fill = "white") +
  scale_shape_manual(values = c(22, 21, 24, 4)) +
  theme_void() + 
  theme(legend.position = c(0.9, 0.8)) +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br", 
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m", 
                              location = "br",
                              pad_y = unit(1, "cm"))

ggsave(here("plots", "plan_gabriel.pdf"), scale = 2)

# get sex of neighbors
ved_g <- graph_from_edgelist(cbind(from = ved_gabriel$from, to = ved_gabriel$to), 
                             directed = FALSE)
V(ved_g)$name <- ved_sf$id_burial

# variable vector
ved_sex <- ved_sf$sex
names(ved_sex) <- ved_sf$id_burial

mean_neighbors <- function(g, variable_vector) {
  res <- vector("list", length(V(g)))
  names(res) <- V(g)$name
  for (i in seq_along(V(g))) {
    x <- attr(neighbors(g, v = i), "names")
    res[[i]] <- unname(variable_vector[x])
  }
  res %>% map(tibble) %>% 
    map(set_names, "to") %>% 
    bind_rows(.id = "id") %>% 
    # mutate(from = unname(variable_vector[id])) %>% 
    group_by(id) %>%
    count(to) %>% 
    mutate(from = unname(variable_vector[id]),
           across(c(from, to), fct_relevel, levels(variable_vector))) %>% 
    ungroup() %>% 
    group_by(from, to) %>% 
    summarise(mean = mean(n), sum = sum(n), .groups = "drop")
}

nb_sex <- mean_neighbors(ved_g, ved_sex)

# randomization of sex for neighbors --------------------------------------
# sf - the attribute table - first column is ID
randomize_neighbors <- function(g, sf, variable, n_sim = 99) {
  n_bur <- nrow(sf)
  variable <- as.symbol(variable)
  
  prob <- sf %>%
    st_drop_geometry() %>% 
    select(!!variable) %>%
    group_by(!!variable) %>%
    count() %>%
    mutate(prob = n / n_bur)
  
  lvls <- prob %>% pull(!!variable) %>% levels()
  
  res <- vector("list", n_sim)
  
  for (i in 1:n_sim) {
    rand_vect <- randomizr::simple_ra(n_bur, prob_each = prob$prob, 
                                      conditions = pull(prob, !!variable)) %>% 
      as.character()
    names(rand_vect) <- sf %>% st_drop_geometry() %>% pull(1)
    res[[i]] <- mean_neighbors(g, rand_vect)
  }
  res %>% bind_rows %>% 
    mutate(across(c(from, to), fct_relevel, lvls))
}

# simulation (99 iterations) ----------------------------------------------

ved_rand_sex <- randomize_neighbors(ved_g, ved_sf, "sex", n_sim = 200)

# plot
ved_rand_sex %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = nb_sex, aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  xlab("mean neighbours") +
  theme_universe

ggsave(here("plots", "nb_gabriel_sex.pdf"), width = 12, height = 6)
