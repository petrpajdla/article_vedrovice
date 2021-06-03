# Project "Vedrovice"
# Script nr. 3.1
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
  geom_sf(data = ved_exc, fill = "gray90", color = NA) +
  geom_sf(data = select(ved_sf, -pres), color = "gray") +
  geom_sf(data = ved_gabriel_lines, linetype = 1, size = 0.2) +
  geom_sf(data = filter(ved_sf, pres == "pres."), shape = 21, fill = "white") +
  theme_void() + 
  theme(legend.position = c(0.9, 0.8),
        plot.title = element_text(hjust = 0.5)) +
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
                                    location = "br",
                                    pad_y = unit(0.8, "cm"),
                                    pad_x = unit(-0.4, "cm"), 
                                    height = unit(1, "cm")) +
  ggspatial::annotation_scale(plot_unit = "m",
                              location = "br", 
                              height = unit(0.2, "cm")) +
  labs(title = "Gabriel graph", shape = "preservation", color = "preservation")

g_delaunay <- ggplot() +
  geom_sf(data = ved_exc, fill = "gray90", color = NA) +
  geom_sf(data = select(ved_sf, -pres), color = "gray") +
  geom_sf(data = ved_delaunay_lines, linetype = 1, size = 0.2) +
  geom_sf(data = filter(ved_sf, pres == "pres."), shape = 21, fill = "white",
          show.legend = FALSE) +
  theme_void() + 
  theme(legend.position = c(0.9, 0.8),
        plot.title = element_text(hjust = 0.5)) +
  # ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_minimal(),
  #                                   location = "br", 
  #                                   pad_y = unit(2, "cm")) +
  # ggspatial::annotation_scale(plot_unit = "m", 
  #                             location = "br",
  #                             pad_y = unit(1, "cm")) +
  labs(title = "Delaunay triangulation", shape = "preservation")

g <- gridExtra::grid.arrange(g_delaunay, g_gabriel, nrow = 1)

ggsave(plot = g, filename = here("plots/plan_graphs.pdf"), 
       width = 19, height = 8, units = "cm")


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

ved_sa <- ved_sf_pres$cat_sa
names(ved_sa) <- ved_sf_pres$id_burial

ved_os <- ved_sf_pres$cat_os
names(ved_os) <- ved_sf_pres$id_burial

ved_oa <- ved_sf_pres$cat_oa
names(ved_oa) <- ved_sf_pres$id_burial

ved_osa <- ved_sf_pres$cat_osa
names(ved_osa) <- ved_sf_pres$id_burial

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

nb_sa_g <- mean_neighbors_networks(ved_g, ved_sa)
nb_sa_d <- mean_neighbors_networks(ved_d, ved_sa)

nb_os_g <- mean_neighbors_networks(ved_g, ved_os)
nb_os_d <- mean_neighbors_networks(ved_d, ved_os)

nb_oa_g <- mean_neighbors_networks(ved_g, ved_oa)
nb_oa_d <- mean_neighbors_networks(ved_d, ved_oa)

nb_osa_g <- mean_neighbors_networks(ved_g, ved_osa)
nb_osa_d <- mean_neighbors_networks(ved_d, ved_osa)

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

ved_rand_sex_g <- randomize_neighbors_network(ved_g, ved_sf_pres, "sex",
                                              n_sim = 999)
write_csv(ved_rand_sex_g, here("data/temp", "ved_rand_sex_g.csv"))

ved_rand_sex_d <- randomize_neighbors_network(ved_d, ved_sf_pres, "sex",
                                              n_sim = 999)
write_csv(ved_rand_sex_d, here("data/temp", "ved_rand_sex_d.csv"))

ved_rand_age_g <- randomize_neighbors_network(ved_g, ved_sf_pres, "age_sim",
                                              n_sim = 999)
write_csv(ved_rand_age_g, here("data/temp", "ved_rand_age_g.csv"))

ved_rand_age_d <- randomize_neighbors_network(ved_d, ved_sf_pres, "age_sim",
                                              n_sim = 999)
write_csv(ved_rand_age_d, here("data/temp", "ved_rand_age_d.csv"))

ved_rand_orig_g <- randomize_neighbors_network(ved_g, ved_sf_pres, "origin",
                                               n_sim = 999)
write_csv(ved_rand_orig_g, here("data/temp", "ved_rand_orig_g.csv"))

ved_rand_orig_d <- randomize_neighbors_network(ved_d, ved_sf_pres, "origin",
                                               n_sim = 999)
write_csv(ved_rand_orig_d, here("data/temp", "ved_rand_orig_d.csv"))

ved_rand_sa_g <- randomize_neighbors_network(ved_g, ved_sf_pres, "cat_sa",
                                             n_sim = 999)
ved_rand_sa_d <- randomize_neighbors_network(ved_d, ved_sf_pres, "cat_sa",
                                             n_sim = 999)
write_csv(ved_rand_sa_g, here("data/temp", "ved_rand_sa_g.csv"))
write_csv(ved_rand_sa_d, here("data/temp", "ved_rand_sa_d.csv"))

ved_rand_os_g <- randomize_neighbors_network(ved_g, ved_sf_pres, "cat_os",
                                             n_sim = 999)
ved_rand_os_d <- randomize_neighbors_network(ved_d, ved_sf_pres, "cat_os",
                                             n_sim = 999)
write_csv(ved_rand_os_g, here("data/temp", "ved_rand_os_g.csv"))
write_csv(ved_rand_os_d, here("data/temp", "ved_rand_os_d.csv"))

ved_rand_oa_g <- randomize_neighbors_network(ved_g, ved_sf_pres, "cat_oa",
                                             n_sim = 999)
ved_rand_oa_d <- randomize_neighbors_network(ved_d, ved_sf_pres, "cat_oa",
                                             n_sim = 999)
write_csv(ved_rand_oa_g, here("data/temp", "ved_rand_oa_g.csv"))
write_csv(ved_rand_oa_d, here("data/temp", "ved_rand_oa_d.csv"))

ved_rand_osa_g <- randomize_neighbors_network(ved_g, ved_sf_pres, "cat_osa",
                                              n_sim = 999)
ved_rand_osa_d <- randomize_neighbors_network(ved_d, ved_sf_pres, "cat_osa",
                                              n_sim = 999)
write_csv(ved_rand_osa_g, here("data/temp", "ved_rand_osa_g.csv"))
write_csv(ved_rand_osa_d, here("data/temp", "ved_rand_osa_d.csv"))

ved_rand_sex_g <- read_csv(here("data/temp", "ved_rand_sex_g.csv"))
ved_rand_sex_d <- read_csv(here("data/temp", "ved_rand_sex_d.csv"))
ved_rand_age_g <- read_csv(here("data/temp", "ved_rand_age_g.csv"))
ved_rand_age_d <- read_csv(here("data/temp", "ved_rand_age_d.csv"))
ved_rand_orig_g <- read_csv(here("data/temp", "ved_rand_orig_g.csv"))
ved_rand_orig_d <- read_csv(here("data/temp", "ved_rand_orig_d.csv"))


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

# combined categories
# sa
g_sa_g <- ved_rand_sa_g %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_sa_g, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "free") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Gabriel graph") +
  theme_universe

g_sa_d <- ved_rand_sa_d %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_sa_d, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "free") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Delaunay triangulation") +
  theme_universe

pdf(here("plots/nb_sa.pdf"), width = 6, height = 8)
gridExtra::grid.arrange(g_sa_d, g_sa_g)
dev.off()

# os
g_os_g <- ved_rand_os_g %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_os_g, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Gabriel graph") +
  theme_universe

g_os_d <- ved_rand_os_d %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_os_d, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Delaunay triangulation") +
  theme_universe

pdf(here("plots/nb_os.pdf"), width = 6, height = 8)
gridExtra::grid.arrange(g_os_d, g_os_g)
dev.off()

# oa
g_oa_g <- ved_rand_oa_g %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_oa_g, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Gabriel graph") +
  theme_universe

g_oa_d <- ved_rand_oa_d %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_oa_d, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Delaunay triangulation") +
  theme_universe

pdf(here("plots/nb_oa.pdf"), width = 6, height = 8)
gridExtra::grid.arrange(g_oa_d, g_oa_g)
dev.off()

# osa
g_osa_g <- ved_rand_osa_g %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_osa_g, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Gabriel graph") +
  theme_universe

g_osa_d <- ved_rand_osa_d %>% 
  filter(from != "ind.", to != "ind.") %>% 
  ggplot(aes(mean)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = filter(nb_osa_d, from != "ind.", to != "ind."), 
             aes(xintercept = mean), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on Delaunay triangulation") +
  theme_universe

pdf(here("plots/nb_osa.pdf"), width = 6, height = 8)
gridExtra::grid.arrange(g_osa_d, g_osa_g)
dev.off()


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

p_vals_multicat <- list(sa = NA,
                        os = NA,
                        oa = NA,
                        osa = NA)

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

# sa
p_sa_d <- pval(nb_sa_d, ved_rand_sa_d) %>% select(-sum)
p_sa_g <- pval(nb_sa_g, ved_rand_sa_g) %>% select(-sum)

p_vals_multicat$sa <- bind_rows(delaunay = p_sa_d, gabriel = p_sa_g, .id = "method") %>% 
  pivot_wider(values_from = c(mean, p, signif), 
              names_from = method, 
              names_glue = "{method}_{.value}") %>% 
  select(from, to, 
         starts_with("delaunay"), 
         starts_with("gabriel")) %>% 
  filter(from != "ind.", to != "ind.")

# os
p_os_d <- pval(nb_os_d, ved_rand_os_d) %>% select(-sum)
p_os_g <- pval(nb_os_g, ved_rand_os_g) %>% select(-sum)

p_vals_multicat$os <- bind_rows(delaunay = p_os_d, gabriel = p_os_g, .id = "method") %>% 
  pivot_wider(values_from = c(mean, p, signif), 
              names_from = method, 
              names_glue = "{method}_{.value}") %>% 
  select(from, to, 
         starts_with("delaunay"), 
         starts_with("gabriel")) %>% 
  filter(from != "ind.", to != "ind.")

# oa
p_oa_d <- pval(nb_oa_d, ved_rand_oa_d) %>% select(-sum)
p_oa_g <- pval(nb_oa_g, ved_rand_oa_g) %>% select(-sum)

p_vals_multicat$oa <- bind_rows(delaunay = p_oa_d, gabriel = p_oa_g, .id = "method") %>% 
  pivot_wider(values_from = c(mean, p, signif), 
              names_from = method, 
              names_glue = "{method}_{.value}") %>% 
  select(from, to, 
         starts_with("delaunay"), 
         starts_with("gabriel")) %>% 
  filter(from != "ind.", to != "ind.")

# osa
p_osa_d <- pval(nb_osa_d, ved_rand_osa_d) %>% select(-sum)
p_osa_g <- pval(nb_osa_g, ved_rand_osa_g) %>% select(-sum)

p_vals_multicat$osa <- bind_rows(delaunay = p_osa_d, gabriel = p_osa_g, .id = "method") %>% 
  pivot_wider(values_from = c(mean, p, signif), 
              names_from = method, 
              names_glue = "{method}_{.value}") %>% 
  select(from, to, 
         starts_with("delaunay"), 
         starts_with("gabriel")) %>% 
  filter(from != "ind.", to != "ind.")

# output
write_rds(p_vals, here("data/temp/nb_pvals.csv"))
write_rds(p_vals_multicat, here("data/temp/nb_pvals_multicat.csv"))

# filtering significant values to assess...
# p_vals_multicat$osa %>% 
#   filter(str_detect(delaunay_signif, "\\*"))
# 
# p_vals_multicat$osa %>% 
#   filter(str_detect(gabriel_signif, "\\*"))

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
