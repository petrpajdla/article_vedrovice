# Project "Vedrovice"
# Script nr. 3.3
# SPATIAL NEIGHBORS: Gabriel graph
# author: Petr Pajdla
# Mean number of neighbours with a given sex based on Gabriel graph

library(here)
library(tidyverse)
library(broom)
library(sf)
library(spdep)
library(igraph)

set.seed(42)

# theme -------------------------------------------------------------------

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

# ved <- read_rds(here("data/temp", "vedrovice_dataset.RDS"))

ved_sf <- read_sf(here("data/temp", "layout.shp")) %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")))
ved_layout <- st_coordinates(ved_sf)
ved_exc <- read_sf(here("data/temp", "window.shp"))


# create delaunay + gabriel neighborhoods ---------------------------------
ved_gabriel <- gabrielneigh(ved_layout)
ved_delaunay <- tri2nb(ved_layout)

# plot delaunay triangulation and gabriel graph ---------------------------
ved_gabriel_lines <- ved_gabriel %>% 
  graph2nb(row.names = rownames(ved_sf)) %>% 
  nb2lines(coords = ved_layout) %>% 
  st_as_sf()

ved_delaunay_lines <- ved_delaunay %>% 
  nb2lines(coords = ved_layout) %>% 
  st_as_sf()

g_gabriel <- ggplot() +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  # geom_sf(data = ved_delaunay_lines, linetype = 1, size = 0.2) +
  geom_sf(data = ved_gabriel_lines, linetype = 1, size = 0.2, color = "gray40") +
  geom_sf(data = ved_sf, aes(shape = sex), fill = "white") +
  scale_shape_manual(values = c(22, 21, 24, 4)) +
  theme_void() + 
  theme(legend.position = c(0.9, 0.8),
        plot.title = element_text(hjust = 0.5)) +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br", 
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m", 
                              location = "br",
                              pad_y = unit(1, "cm")) +
  labs(title = "Gabriel graph")

g_delaunay <- ggplot() +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  geom_sf(data = ved_delaunay_lines, linetype = 1, size = 0.2, color = "gray40") +
  # geom_sf(data = ved_gabriel_lines, linetype = 3) +
  geom_sf(data = ved_sf, aes(shape = sex), fill = "white") +
  scale_shape_manual(values = c(22, 21, 24, 4)) +
  theme_void() + 
  theme(legend.position = c(0.9, 0.8),
        plot.title = element_text(hjust = 0.5)) +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br", 
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m", 
                              location = "br",
                              pad_y = unit(1, "cm")) +
  labs(title = "Delaunay triangulation")

g <- gridExtra::grid.arrange(g_delaunay, g_gabriel, nrow = 1)

ggsave(plot = g, filename = here("plots/plan_gabriel.pdf"), width = 10.5, height = 5)


# graph objects for delaunay and gabriel ----------------------------------
# gabriel
ved_g <- graph_from_edgelist(cbind(from = ved_gabriel$from, to = ved_gabriel$to), 
                             directed = FALSE)
V(ved_g)$name <- ved_sf$id_burial

# delaunay
ved_d <- graph_from_adj_list(ved_delaunay)
V(ved_d)$name <- ved_sf$id_burial

# variable vector of categorical variable --------------------------------
ved_sex <- ved_sf$sex
names(ved_sex) <- ved_sf$id_burial

# count mean neighbors ----------------------------------------------------
mean_neighbors_networks <- function(g, variable_vector) {
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

nb_sex_g <- mean_neighbors_networks(ved_g, ved_sex)
nb_sex_d <- mean_neighbors_networks(ved_d, ved_sex)

# randomization of sex for neighbors --------------------------------------
# sf - the attribute table - first column is ID
randomize_neighbors_network <- function(g, sf, variable, n_sim = 99) {
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
    res[[i]] <- mean_neighbors_networks(g, rand_vect)
  }
  res %>% bind_rows %>% 
    mutate(across(c(from, to), fct_relevel, lvls))
}

# simulation (200 iterations) ----------------------------------------------

ved_rand_sex_g <- randomize_neighbors_network(ved_g, ved_sf, "sex", n_sim = 200)
ved_rand_sex_d <- randomize_neighbors_network(ved_d, ved_sf, "sex", n_sim = 200)

# plots -------------------------------------------------------------------
g_d <- ved_rand_sex_g %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = nb_sex_g, aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Gabriel graph") +
  theme_universe

ggsave(here("plots/nb_gabriel_sex.pdf"), g_d, width = 12, height = 6)

g_g <- ved_rand_sex_d %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = nb_sex_d, aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Delaunay triangulation") +
  theme_universe

ggsave(here("plots/nb_delaunay_sex.pdf"), g_g,  width = 12, height = 6)


# neighborhood based on buffer --------------------------------------------
# buffer zone -------------------------------------------------------------
# 6 m (based on NN functions)
distance <- 6

mean_neighbors_buffer <- function(sf, variable, dist) {
  variable <- as.symbol(variable)
  lvls <- levels(factor(pull(sf, !!variable)))
  res <- matrix(nrow = length(lvls), ncol = length(lvls))
  colnames(res) <- lvls
  rownames(res) <- lvls
  # get buffer for each level of a variable
  for (i in seq_along(lvls)) {
    buffer <- sf %>% filter(!!variable == lvls[i]) %>% 
      st_buffer(dist = dist) %>% 
      st_geometry()
    all_nbs <- sum(st_intersects(buffer, st_geometry(sf), sparse = FALSE)) - length(buffer)
    for (j in seq_along(lvls)) {
      points <- sf %>% filter(!!variable == lvls[j]) %>% 
        st_geometry()
      intersect <- sum(st_intersects(buffer, points, sparse = FALSE))
      if (lvls[i] == lvls[j]) {
        res[i, j] <- (intersect - length(buffer)) / all_nbs
      } else {
        res[i, j] <- intersect / all_nbs
      }
    }
  }
  as.data.frame(res) %>% 
    rownames_to_column(var = "from") %>% 
    pivot_longer(-from, names_to = "to") %>% 
    mutate(across(c(from, to), fct_relevel, lvls)) %>% 
    rename(mean = value)
}

observed <- mean_neighbors_buffer(ved_sf, "sex", distance)
# mean_neighbors_buffer(ved_sf, "age_sim", distance)

# resampling --------------------------------------------------------------

randomize_neighbors_buffer <- function(sf, variable, dist, n_sim = 99) {
  n_bur <- nrow(sf)
  variable <- as.symbol(variable)
  lvls <- levels(factor(pull(sf, !!variable)))
  
  prob <- sf %>%
    st_drop_geometry() %>% 
    select(!!variable) %>%
    group_by(!!variable) %>%
    count() %>%
    mutate(prob = n / n_bur)
  
  res <- vector("list", n_sim)
  
  for (i in 1:n_sim) {
    rand_lvls <- randomizr::simple_ra(n_bur, prob_each = prob$prob, 
                                      conditions = pull(prob, !!variable))
    
    res[[i]] <- sf %>% bind_cols(rand = rand_lvls) %>% 
      mean_neighbors_buffer("rand", dist)
  }
  res %>% bind_rows()
}

expected <- randomize_neighbors_buffer(ved_sf, "sex", distance, n_sim = 200)

# plot
g_b <- expected %>% ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = observed, aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on buffer") +
  theme_universe

ggsave(here("plots", "nb_buffer_sex.pdf"), g_b, width = 12, height = 6)

# plot buffer -------------------------------------------------------------

# ved_buff <- st_buffer(ved_sf, dist = distance)
# 
# ggplot(data = ved_sf) +
#   geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
#   geom_sf(data = ved_buff, fill = "gray80", color = "gray80", alpha = 0.2) +
#   geom_sf(aes(shape = sex), fill = "white", size = 2) +
#   scale_shape_manual(values = c(22, 21, 24, 4)) +
#   ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
#                                     location = "br",
#                                     pad_y = unit(2, "cm")) +
#   ggspatial::annotation_scale(plot_unit = "m",
#                               location = "br",
#                               pad_y = unit(1, "cm")) +
#   labs(shape = "body sex") +
#   facet_wrap(vars(sex)) +
#   theme_void()



# export combined plot ----------------------------------------------------

g_full <- gridExtra::grid.arrange(g_b, g_d, g_g, nrow = 3)

ggsave(here("plots/nb_sex.pdf"), g_full, width = 12, height = 18)

# p-value -----------------------------------------------------------------
pval <- function(obs, exp) {
  lvls <- levels(obs$from)
  res <- obs %>% mutate(p = 0)
  for (i in seq_along(lvls)) {
    for (j in seq_along(lvls)) {
      obs_val <- obs %>% filter(from == lvls[i], to == lvls[j]) %>% 
        pull(mean)
      exp_val <- exp %>% filter(from == lvls[i], to == lvls[j]) %>% 
        pull(mean)
      res[(res[, "from"] == lvls[i]) & (res[, "to"] == lvls[j]), "p"] <- mean(exp_val >= obs_val)
    }
  }
  res <- res %>% mutate(signif = if_else(p <= 0.05, "*", ""),
                        signif = if_else(p <= 0.01, "**", signif),
                        signif = if_else(p <= 0.001, "***", signif),
                        across(where(is.numeric), round, 2),
                        p = str_c(p, signif, sep = "")) %>% 
    select(-signif)
  return(res)
}

p_d <- pval(nb_sex_d, ved_rand_sex_d) %>% select(-sum)
p_g <- pval(nb_sex_g, ved_rand_sex_g) %>% select(-sum)
p_b <- pval(observed, expected)

bind_rows(delaunay = p_d, gabriel = p_g, buffer = p_b, .id = "method") %>% 
  pivot_wider(values_from = c(mean, p), 
              names_from = method, 
              names_glue = "{method}_{.value}") %>% 
  select(from, to, 
         starts_with("buffer"), 
         starts_with("delaunay"), 
         starts_with("gabriel")) %>% 
  write_csv(here("data/temp/nb_sex_p.csv"))


