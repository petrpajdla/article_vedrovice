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
# ei <- read_csv(here("data/temp", "exceptionality.csv"), col_types = "cdff") %>% 
#   select(id_burial = burial, ei, ei_clust = fct)
ved_sf <- read_sf(here("data/temp", "layout.shp")) %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")),
         age_sim = fct_relevel(age_sim, c("juv.", "ad.", "mat.", "ind.")),
         origin = fct_relevel(origin, c("local", "non-local", "ind.")))
ved_layout <- ved_sf %>% 
  filter(pres != "dist.") %>% 
  st_coordinates()
ved_exc <- read_sf(here("data/temp", "window.shp"))

ved_sf_pres <- filter(ved_sf, pres != "dist.") %>% 
  mutate(age_sim = forcats::fct_drop(age_sim))

# create delaunay + gabriel neighborhoods ---------------------------------
ved_gabriel <- gabrielneigh(ved_layout)
ved_delaunay <- tri2nb(ved_layout)

# plot delaunay triangulation and gabriel graph ---------------------------
ved_gabriel_lines <- ved_gabriel %>% 
  graph2nb(row.names = ved_sf_pres$id_burial) %>% 
  nb2lines(coords = ved_layout) %>% 
  st_as_sf()

ved_delaunay_lines <- ved_delaunay %>% 
  nb2lines(coords = ved_layout) %>% 
  st_as_sf()

g_gabriel <- ggplot() +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  # geom_sf(data = ved_delaunay_lines, linetype = 1, size = 0.2) +
  geom_sf(data = ved_gabriel_lines, linetype = 1, size = 0.2, color = "gray40") +
  geom_sf(data = ved_sf, aes(shape = pres), fill = "white") +
  scale_shape_manual(values = c(15, 21)) +
  theme_void() + 
  theme(legend.position = c(0.9, 0.8),
        plot.title = element_text(hjust = 0.5)) +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br", 
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m", 
                              location = "br",
                              pad_y = unit(1, "cm")) +
  labs(title = "Gabriel graph", shape = "preservation")

g_delaunay <- ggplot() +
  geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
  geom_sf(data = ved_delaunay_lines, linetype = 1, size = 0.2, color = "gray40") +
  # geom_sf(data = ved_gabriel_lines, linetype = 3) +
  geom_sf(data = ved_sf, aes(shape = pres), fill = "white") +
  scale_shape_manual(values = c(15, 21)) +
  theme_void() + 
  theme(legend.position = c(0.9, 0.8),
        plot.title = element_text(hjust = 0.5)) +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br", 
                                    pad_y = unit(2, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m", 
                              location = "br",
                              pad_y = unit(1, "cm")) +
  labs(title = "Delaunay triangulation", shape = "preservation")

g <- gridExtra::grid.arrange(g_delaunay, g_gabriel, nrow = 1)

ggsave(plot = g, filename = here("plots/plan_graphs.pdf"), width = 10.5, height = 5)


# graph objects for delaunay and gabriel ----------------------------------
# gabriel
ved_g <- graph_from_edgelist(cbind(from = ved_gabriel$from, to = ved_gabriel$to), 
                             directed = FALSE)
V(ved_g)$name <- ved_sf_pres$id_burial

# delaunay
ved_d <- graph_from_adj_list(ved_delaunay)
V(ved_d)$name <- ved_sf_pres$id_burial

# variable vector of categorical variables --------------------------------
ved_sex <- ved_sf_pres$sex
names(ved_sex) <- ved_sf_pres$id_burial

ved_age <- ved_sf_pres$age_sim
names(ved_age) <- ved_sf_pres$id_burial

ved_origin <- ved_sf_pres$origin
names(ved_origin) <- ved_sf_pres$id_burial

# ved_ei <- ved_sf_pres$ei_clust
# names(ved_ei) <- ved_sf_pres$id_burial

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

nb_age_g <- mean_neighbors_networks(ved_g, ved_age)
nb_age_d <- mean_neighbors_networks(ved_d, ved_age)

nb_orig_g <- mean_neighbors_networks(ved_g, ved_origin)
nb_orig_d <- mean_neighbors_networks(ved_d, ved_origin)

# nb_ei_g <- mean_neighbors_networks(ved_g, ved_ei)
# nb_ei_d <- mean_neighbors_networks(ved_d, ved_ei)

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

# simulation (999 iterations) ----------------------------------------------

# ved_rand_sex_g <- randomize_neighbors_network(ved_g, ved_sf_pres, "sex", 
#                                               n_sim = 999)
# write_csv(ved_rand_sex_g, here("data/temp", "ved_rand_sex_g.csv"))
# 
# ved_rand_sex_d <- randomize_neighbors_network(ved_d, ved_sf_pres, "sex", 
#                                               n_sim = 999)
# write_csv(ved_rand_sex_d, here("data/temp", "ved_rand_sex_d.csv"))
# 
# ved_rand_age_g <- randomize_neighbors_network(ved_g, ved_sf_pres, "age_sim", 
#                                               n_sim = 999)
# write_csv(ved_rand_age_g, here("data/temp", "ved_rand_age_g.csv"))
# 
# ved_rand_age_d <- randomize_neighbors_network(ved_d, ved_sf_pres, "age_sim", 
#                                               n_sim = 999)
# write_csv(ved_rand_age_d, here("data/temp", "ved_rand_age_d.csv"))
# 
# ved_rand_orig_g <- randomize_neighbors_network(ved_g, ved_sf_pres, "origin", 
#                                                n_sim = 999)
# write_csv(ved_rand_orig_g, here("data/temp", "ved_rand_orig_g.csv"))
# 
# ved_rand_orig_d <- randomize_neighbors_network(ved_d, ved_sf_pres, "origin", 
#                                                n_sim = 999)
# write_csv(ved_rand_orig_d, here("data/temp", "ved_rand_orig_d.csv"))

# ved_rand_ei_g <- randomize_neighbors_network(ved_g, ved_sf, "ei_clust", n_sim = 999)
# write_csv(ved_rand_ei_g, here("data/temp", "ved_rand_ei_g.csv"))
# 
# ved_rand_ei_d <- randomize_neighbors_network(ved_d, ved_sf, "ei_clust", n_sim = 999)
# write_csv(ved_rand_ei_d, here("data/temp", "ved_rand_ei_d.csv"))

ved_rand_sex_g <- read_csv(here("data/temp", "ved_rand_sex_g.csv"))
ved_rand_sex_d <- read_csv(here("data/temp", "ved_rand_sex_d.csv"))
ved_rand_age_g <- read_csv(here("data/temp", "ved_rand_age_g.csv"))
ved_rand_age_d <- read_csv(here("data/temp", "ved_rand_age_d.csv"))
ved_rand_orig_g <- read_csv(here("data/temp", "ved_rand_orig_g.csv"))
ved_rand_orig_d <- read_csv(here("data/temp", "ved_rand_orig_d.csv"))
# ved_rand_ei_g <- read_csv(here("data/temp", "ved_rand_ei_g.csv"))
# ved_rand_ei_d <- read_csv(here("data/temp", "ved_rand_ei_d.csv"))


# plots -------------------------------------------------------------------
# sex
g_sex_g <- ved_rand_sex_g %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_sex_g, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Gabriel graph") +
  theme_universe

# ggsave(here("plots/nb_gabriel_sex.pdf"), g_sex_g, width = 9, height = 6)

g_sex_d <- ved_rand_sex_d %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_sex_d, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Delaunay triangulation") +
  theme_universe

# ggsave(here("plots/nb_delaunay_sex.pdf"), g_sex_d,  width = 9, height = 6)

pdf(here("plots/nb_sex.pdf"), width = 9, height = 12)
gridExtra::grid.arrange(g_sex_d, g_sex_g)
dev.off()

# age
g_age_g <- ved_rand_age_g %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_age_g, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Gabriel graph") +
  theme_universe

# ggsave(here("plots/nb_gabriel_age.pdf"), g_age_g, width = 9, height = 6)

g_age_d <- ved_rand_age_d %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_age_d, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Delaunay triangulation") +
  theme_universe

# ggsave(here("plots/nb_delaunay_age.pdf"), g_age_d,  width = 9, height = 6)

pdf(here("plots/nb_age.pdf"), width = 9, height = 12)
gridExtra::grid.arrange(g_age_d, g_age_g)
dev.off()

# localness
g_orig_g <- ved_rand_orig_g %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_orig_g, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Gabriel graph") +
  theme_universe

# ggsave(here("plots/nb_gabriel_orig.pdf"), g_orig_g, width = 6, height = 4)

g_orig_d <- ved_rand_orig_d %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_orig_d, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Delaunay triangulation") +
  theme_universe

# ggsave(here("plots/nb_delaunay_orig.pdf"), g_orig_d,  width = 6, height = 4)

pdf(here("plots/nb_orig.pdf"), width = 6, height = 8)
gridExtra::grid.arrange(g_orig_d, g_orig_g)
dev.off()

# # exceptionality index
# g_ei_g <- ved_rand_ei_g %>% 
#   filter(from != "ind.", to != "ind.") %>% 
#   ggplot(aes(mean)) +
#   geom_density(fill = "gray90", color = NA) +
#   geom_vline(data = filter(nb_ei_g, from != "ind.", to != "ind."), 
#              aes(xintercept = mean), size = 0.8) +
#   geom_rug() +
#   facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
#   scale_y_continuous(expand = c(0, 0, 0, 0)) +
#   scale_x_continuous(expand = c(0.01, 0)) +
#   labs(xlab = "mean neighbours", title = "Neighborhood based on Gabriel graph") +
#   theme_universe
# 
# ggsave(here("plots/nb_gabriel_ei.pdf"), g_ei_g, width = 6, height = 4)
# 
# g_ei_d <- ved_rand_ei_d %>% 
#   filter(from != "ind.", to != "ind.") %>% 
#   ggplot(aes(mean)) +
#   geom_density(fill = "gray90", color = NA) +
#   geom_vline(data = filter(nb_ei_d, from != "ind.", to != "ind."), 
#              aes(xintercept = mean), size = 0.8) +
#   geom_rug() +
#   facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
#   scale_y_continuous(expand = c(0, 0, 0, 0)) +
#   scale_x_continuous(expand = c(0.01, 0)) +
#   labs(xlab = "mean neighbours", title = "Neighborhood based on Delaunay triangulation") +
#   theme_universe
# 
# ggsave(here("plots/nb_delaunay_ei.pdf"), g_ei_d,  width = 6, height = 4)


# # neighborhood based on buffer --------------------------------------------
# # buffer zone -------------------------------------------------------------
# # 6 m (based on NN functions)
# distance <- 6
# 
# mean_neighbors_buffer <- function(sf, variable, dist) {
#   variable <- as.symbol(variable)
#   lvls <- levels(factor(pull(sf, !!variable)))
#   res <- matrix(nrow = length(lvls), ncol = length(lvls))
#   colnames(res) <- lvls
#   rownames(res) <- lvls
#   # get buffer for each level of a variable
#   for (i in seq_along(lvls)) {
#     buffer <- sf %>% filter(!!variable == lvls[i]) %>% 
#       st_buffer(dist = dist) %>% 
#       st_geometry()
#     all_nbs <- sum(st_intersects(buffer, st_geometry(sf), sparse = FALSE)) - length(buffer)
#     for (j in seq_along(lvls)) {
#       points <- sf %>% filter(!!variable == lvls[j]) %>% 
#         st_geometry()
#       intersect <- sum(st_intersects(buffer, points, sparse = FALSE))
#       if (lvls[i] == lvls[j]) {
#         res[i, j] <- (intersect - length(buffer)) / all_nbs
#       } else {
#         res[i, j] <- intersect / all_nbs
#       }
#     }
#   }
#   as.data.frame(res) %>% 
#     rownames_to_column(var = "from") %>% 
#     pivot_longer(-from, names_to = "to") %>% 
#     mutate(across(c(from, to), fct_relevel, lvls)) %>% 
#     rename(mean = value)
# }
# 
# observed <- mean_neighbors_buffer(ved_sf, "sex", distance)
# # mean_neighbors_buffer(ved_sf, "age_sim", distance)
# 
# # resampling --------------------------------------------------------------
# 
# randomize_neighbors_buffer <- function(sf, variable, dist, n_sim = 99) {
#   n_bur <- nrow(sf)
#   variable <- as.symbol(variable)
#   lvls <- levels(factor(pull(sf, !!variable)))
#   
#   prob <- sf %>%
#     st_drop_geometry() %>% 
#     select(!!variable) %>%
#     group_by(!!variable) %>%
#     count() %>%
#     mutate(prob = n / n_bur)
#   
#   res <- vector("list", n_sim)
#   
#   for (i in 1:n_sim) {
#     rand_lvls <- randomizr::simple_ra(n_bur, prob_each = prob$prob, 
#                                       conditions = pull(prob, !!variable))
#     
#     res[[i]] <- sf %>% bind_cols(rand = rand_lvls) %>% 
#       mean_neighbors_buffer("rand", dist)
#   }
#   res %>% bind_rows()
# }
# 
# expected <- randomize_neighbors_buffer(ved_sf, "sex", distance, n_sim = 200)
# 
# # plot
# g_b <- expected %>% ggplot(aes(mean)) +
#   geom_density(fill = "gray90", color = NA) +
#   geom_vline(data = observed, aes(xintercept = mean), size = 0.8) +
#   geom_rug() +
#   facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
#   scale_y_continuous(expand = c(0, 0, 0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   labs(xlab = "mean neighbours", title = "Neighborhood based on buffer") +
#   theme_universe
# 
# ggsave(here("plots", "nb_buffer_sex.pdf"), g_b, width = 12, height = 6)
# 
# # plot buffer -------------------------------------------------------------
# 
# # ved_buff <- st_buffer(ved_sf, dist = distance)
# # 
# # ggplot(data = ved_sf) +
# #   geom_sf(data = ved_exc, fill = NA, color = "gray90", size = 4) +
# #   geom_sf(data = ved_buff, fill = "gray80", color = "gray80", alpha = 0.2) +
# #   geom_sf(aes(shape = sex), fill = "white", size = 2) +
# #   scale_shape_manual(values = c(22, 21, 24, 4)) +
# #   ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
# #                                     location = "br",
# #                                     pad_y = unit(2, "cm")) +
# #   ggspatial::annotation_scale(plot_unit = "m",
# #                               location = "br",
# #                               pad_y = unit(1, "cm")) +
# #   labs(shape = "body sex") +
# #   facet_wrap(vars(sex)) +
# #   theme_void()


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
      res[(res[, "from"] == lvls[i]) & (res[, "to"] == lvls[j]), "p"] <- round(mean(exp_val >= obs_val), 2)
    }
  }
  res <- res %>% mutate(signif = if_else(p <= 0.05, "*", ""),
                        signif = if_else(p <= 0.01, "**", signif),
                        signif = if_else(p <= 0.001, "***", signif),
                        across(where(is.numeric), round, 2))
  # p = str_c(p, signif, sep = "")) %>% 
  # select(-signif)
  return(res)
}

# output of p-vals
p_vals <- list(sex = NA,
               age = NA,
               orig = NA)


# sex
p_sex_d <- pval(nb_sex_d, ved_rand_sex_d) %>% select(-sum)
p_sex_g <- pval(nb_sex_g, ved_rand_sex_g) %>% select(-sum)

p_vals$sex <- bind_rows(delaunay = p_sex_d, gabriel = p_sex_g, .id = "method") %>% 
  pivot_wider(values_from = c(mean, p, signif), 
              names_from = method, 
              names_glue = "{method}_{.value}") %>% 
  select(from, to, 
         starts_with("delaunay"), 
         starts_with("gabriel")) %>% 
  filter(from != "ind.", to != "ind.")

# age
p_age_d <- pval(nb_age_d, ved_rand_age_d) %>% select(-sum)
p_age_g <- pval(nb_age_g, ved_rand_age_g) %>% select(-sum)

p_vals$age <- bind_rows(delaunay = p_age_d, gabriel = p_age_g, .id = "method") %>% 
  pivot_wider(values_from = c(mean, p, signif), 
              names_from = method, 
              names_glue = "{method}_{.value}") %>% 
  select(from, to, 
         starts_with("delaunay"), 
         starts_with("gabriel"))

# origin
p_orig_d <- pval(nb_orig_d, ved_rand_orig_d) %>% select(-sum)
p_orig_g <- pval(nb_orig_g, ved_rand_orig_g) %>% select(-sum)

p_vals$orig <- bind_rows(delaunay = p_orig_d, gabriel = p_orig_g, .id = "method") %>% 
  pivot_wider(values_from = c(mean, p, signif), 
              names_from = method, 
              names_glue = "{method}_{.value}") %>% 
  select(from, to, 
         starts_with("delaunay"), 
         starts_with("gabriel")) %>% 
  filter(from != "ind.", to != "ind.")

write_rds(p_vals, here("data/temp/nb_pvals.csv"))

# # ei
# p_ei_d <- pval(nb_ei_d, ved_rand_ei_d) %>% select(-sum)
# p_ei_g <- pval(nb_ei_g, ved_rand_ei_g) %>% select(-sum)
# 
# bind_rows(delaunay = p_ei_d, gabriel = p_ei_g, .id = "method") %>% 
#   pivot_wider(values_from = c(mean, p), 
#               names_from = method, 
#               names_glue = "{method}_{.value}") %>% 
#   select(from, to, 
#          starts_with("delaunay"), 
#          starts_with("gabriel")) %>% 
#   write_csv(here("data/temp/nb_ei_p.csv"))
