# Project "Vedrovice"
# Script nr. 3.4
# SPATIAL NEIGHBORS: Buffer zone
# author: Petr Pajdla
# Prevalent sex at a certain distance (buffer)

library(here)
library(sf)
library(tidyverse)

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

ved_sf <- read_sf(here("data/temp", "layout.shp")) %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")))
ved_exc <- read_sf(here("data/temp", "window.shp"))

# buffer zone -------------------------------------------------------------
# 6 m (based on NN functions)
distance <- 6

mean_neighbors <- function(sf, variable, dist) {
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
    mutate(across(c(from, to), fct_relevel, lvls))
}

observed <- mean_neighbors(ved_sf, "sex", distance)
# mean_neighbors(ved_sf, "age_sim", distance)

# resampling --------------------------------------------------------------

randomize_neighbors <- function(sf, variable, dist, n_sim = 99) {
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
      mean_neighbors("rand", dist)
  }
  res %>% bind_rows()
}

expected <- randomize_neighbors(ved_sf, "sex", distance, n_sim = 200)

# plot
expected %>% ggplot(aes(value)) +
  geom_density(fill = "gray90", color = NA) +
  geom_vline(data = observed, aes(xintercept = value), size = 0.8) +
  geom_rug() +
  facet_grid(rows = vars(to), cols = vars(from), scales = "fixed") +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(xlab = "mean neighbours", title = "Neighborhood based on buffer") +
  theme_universe

ggsave(here("plots", "nb_buffer_sex.pdf"), width = 12, height = 6)

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


# p-value -----------------------------------------------------------------
pval <- function(obs, exp) {
  lvls <- levels(obs$from)
  res <- obs %>% mutate(p = 0)
  for (i in seq_along(lvls)) {
    for (j in seq_along(lvls)) {
      obs_val <- obs %>% filter(from == lvls[i], to == lvls[j]) %>% 
        pull(value)
      exp_val <- exp %>% filter(from == lvls[i], to == lvls[j]) %>% 
        pull(value)
      res[(res[, "from"] == lvls[i]) & (res[, "to"] == lvls[j]), "p"] <- mean(exp_val >= obs_val)
    }
  }
  res <- res %>% mutate(signif = if_else(p <= 0.05, TRUE, FALSE))
  return(res)
}

pval(observed, expected) %>% write_csv(here("data/temp/nb_p_buffer.csv"))
