# Project "Vedrovice"
# Script nr. 3.4
# SPATIAL ARRANGEMENT - Grouping of vars into clusters derived from percolation
# author: Petr Pajdla
# Identify deviations from random for different percolation clusters 

set.seed(42)

# number of simulations
nsim <- 999

# libs --------------------------------------------------------------------

library(tidyverse)
library(sf)


# data --------------------------------------------------------------------

explored_rad <- paste0(c(2.1, 2.8, 3.4, 4.1, 4.6, 5.6, 6.3, 8.8, 9.2))

perc_cl <- read_csv(here::here("data/temp", "perc_clusters.csv")) %>% 
  mutate(id_burial = as.character(id)) %>% 
  select(-id) %>% 
  filter(radius %in% explored_rad)

ved_sf <- sf::st_read(here::here("data/temp", "layout.geojson")) %>% 
  st_set_crs(NA) %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")),
         age_sim = fct_relevel(age_sim, c("juv.", "ad.", "mat.", "ind.")),
         origin = fct_relevel(origin, c("local", "non-local", "ind.")))

# prep ---------------------------------------------------------------------

input <- ved_sf %>% 
  select(id_burial, sex, age_sim, origin, 
         # starts_with("cat_")
  ) %>% 
  st_drop_geometry() %>% 
  left_join(perc_cl, by = c("id_burial")) %>% 
  filter(!is.na(radius), !is.na(cluster)) %>% 
  mutate(rc = paste0(radius, "-", cluster))

# count only for clusters with >= 4 members
clust_size <- input %>% 
  group_by(radius, cluster) %>% 
  count() %>% 
  # filter(n >= 4) %>% 
  mutate(rc = paste0(radius, "-", cluster))

input <- input %>% 
  filter(rc %in% clust_size$rc) %>% 
  select(-rc)

# functions ---------------------------------------------------------------

count_var_per_radius <- function(input, variable) {
  res <- input %>% 
    filter(!is.na(radius), !is.na(cluster)) %>%
    group_by(radius) %>% 
    nest() %>% 
    mutate(data = map(data, group_by, cluster, !!sym(variable)),
           data = map(data, count)) %>% 
    unnest(data) %>% 
    ungroup(radius) %>% 
    setNames(c("radius", "cluster", "variable", "n"))
  return(res)
}

randomize_var_per_radius <- function(input, variable, n_sim) {
  n_bur <- nrow(input)
  variable0 <- as.symbol(variable)
  
  prob <- input %>%
    select(!!variable0) %>%
    group_by(!!variable0) %>%
    count() %>%
    ungroup(!!variable0) %>% 
    mutate(prob = n / n_bur)
  
  lvls <- prob %>% 
    pull(!!variable0) %>% 
    levels()
  
  res <- vector("list", n_sim)
  # pb <- txtProgressBar(max = length(n_sim), style = 3)
  
  for (i in 1:n_sim) {
    rand_vect <- randomizr::simple_ra(n_bur, prob_each = prob$prob, 
                                      conditions = pull(prob, !!variable0)) %>% 
      as.character()
    
    input <- input %>% 
      select(id_burial, radius, cluster) %>% 
      mutate(temp = rand_vect)
    
    res[[i]] <- count_var_per_radius(input, "temp")
    # setTxtProgressBar(pb, i)
  }
  res %>% 
    bind_rows() # %>% 
  # group_by(radius, cluster, variable) %>%
  # summarise(n = mean(n, na.rm = TRUE), .groups = "drop") %>%
  # ungroup()
}

return_deviations_from_random <- function(expected, observed) {
  
  expected %>% 
    group_by(radius, cluster, variable) %>% 
    nest() %>%
    # mutate(data = map(data, \(x) bind_rows(x, tibble(n = rep(1, 999 - nrow(x)))))) %>% 
    left_join(observed, by = c("radius", "cluster", "variable")) %>% 
    mutate(n = if_else(is.na(n), 0L, n),
           larger = map2(data, n, \(x, y) filter(x, n >= y)),
           smaller = map2(data, n, \(x, y) filter(x, n <= y)),
           larger_pct = map2_dbl(data, larger, \(x, y) (nrow(y) / nrow(x))*100),
           smaller_pct = map2_dbl(data, smaller, \(x, y) (nrow(y) / nrow(x))*100)) %>% 
    filter(larger_pct < 5 | smaller_pct < 5) %>% 
    select(-data, -larger, -smaller)
  
  # expected %>% 
  #   group_by(radius, cluster, variable) %>% 
  #   summarise(m = mean(n, na.rm = TRUE),
  #             sd= sd(n, na.rm = TRUE),
  #             max = m + 2*sd,
  #             min = m - 2*sd, .groups = "drop") %>% 
  #   inner_join(observed, by = c("radius", "cluster", "variable")) %>% 
  #   mutate(dev = n - m,
  #          flt = (n >= max)|(n <= min)) %>% 
  #   filter(flt) %>% 
  #   select(-flt)
  
  # bind_rows(expected = expected, observed = observed, .id = "origin") %>% 
  #   group_by(radius) %>% 
  #   nest() %>%
  #   mutate(data = map(data, pivot_wider, names_from = origin, values_from = n),
  #          data = map(data, mutate, observed = if_else(is.na(observed), 0, observed)),
  #          data = map(data, mutate, deviation = observed - expected),
  #          data = map(data, group_by, variable),
  #          sums = map(data, summarise, 
  #                     m = mean(deviation),
  #                     sd = sd(deviation),
  #                     max = m + 2*sd,
  #                     min = m - 2*sd),
  #          full = map2(data, sums, left_join, by = "variable")) %>%
  #   select(full) %>% 
  #   unnest(cols = full) %>% 
  #   ungroup(radius) %>% 
  #   mutate(flt = (deviation >= max)|(deviation <= min)) %>% 
  #   filter(flt) %>% 
  #   select(radius, cluster, variable, expected, observed, deviation) %>% 
  #   arrange(desc(deviation))
}

deviates_from_random <- function(input, variable, n_sim) {
  x <- Sys.time()
  observed <- count_var_per_radius(input, variable)
  expected <- randomize_var_per_radius(input, variable, n_sim)
  res <- return_deviations_from_random(expected, observed)
  print(Sys.time() - x)
  return(res)
}


# analysis ----------------------------------------------------------------

# larger_pct gives a percentage of larger values than is expected
# smaller_pct gives a percentage of smaller numbers than is expected

devs_sex <- deviates_from_random(input, "sex", nsim)
devs_age <- deviates_from_random(input, "age_sim", nsim)
devs_loc <- deviates_from_random(input, "origin", nsim)

# # multicategorical variables
# devs_sa <- inner_join(deviates_from_random(input, "cat_sa", nsim),
#                       deviates_from_random(input, "cat_sa", nsim), 
#                       by = c("radius", "cluster", "variable")) %>% 
#   transmute(radius, cluster, variable,
#             observed = observed.x,
#             expected = (expected.x + expected.y)/2,
#             deviation = (deviation.x + deviation.y)/2) %>% 
#   inner_join(deviates_from_random(input, "cat_sa", nsim),
#              by = c("radius", "cluster", "variable")) %>% 
#   transmute(radius, cluster, variable,
#             observed = observed.x,
#             expected = (expected.x + expected.y)/2,
#             deviation = (deviation.x + deviation.y)/2)
# 
# devs_os <- inner_join(deviates_from_random(input, "cat_os", nsim),
#                       deviates_from_random(input, "cat_os", nsim), 
#                       by = c("radius", "cluster", "variable")) %>% 
#   transmute(radius, cluster, variable,
#             observed = observed.x,
#             expected = (expected.x + expected.y)/2,
#             deviation = (deviation.x + deviation.y)/2) %>% 
#   inner_join(deviates_from_random(input, "cat_os", nsim),
#              by = c("radius", "cluster", "variable")) %>% 
#   transmute(radius, cluster, variable,
#             observed = observed.x,
#             expected = (expected.x + expected.y)/2,
#             deviation = (deviation.x + deviation.y)/2)
# 
# devs_oa <- inner_join(deviates_from_random(input, "cat_oa", nsim),
#                       deviates_from_random(input, "cat_oa", nsim), 
#                       by = c("radius", "cluster", "variable")) %>% 
#   transmute(radius, cluster, variable,
#             observed = observed.x,
#             expected = (expected.x + expected.y)/2,
#             deviation = (deviation.x + deviation.y)/2) %>% 
#   inner_join(deviates_from_random(input, "cat_oa", nsim),
#              by = c("radius", "cluster", "variable")) %>% 
#   transmute(radius, cluster, variable,
#             observed = observed.x,
#             expected = (expected.x + expected.y)/2,
#             deviation = (deviation.x + deviation.y)/2)
# 
# devs_osa <- inner_join(deviates_from_random(input, "cat_osa", nsim),
#                        deviates_from_random(input, "cat_osa", nsim), 
#                        by = c("radius", "cluster", "variable")) %>% 
#   transmute(radius, cluster, variable,
#             observed = observed.x,
#             expected = (expected.x + expected.y)/2,
#             deviation = (deviation.x + deviation.y)/2) %>% 
#   inner_join(deviates_from_random(input, "cat_osa", nsim),
#              by = c("radius", "cluster", "variable")) %>% 
#   transmute(radius, cluster, variable,
#             observed = observed.x,
#             expected = (expected.x + expected.y)/2,
#             deviation = (deviation.x + deviation.y)/2)


# output ------------------------------------------------------------------

out <- list(sex = devs_sex,
            age = devs_age,
            orig = devs_loc)

# out_multicat <- list(sa = devs_sa,
#                      os = devs_os,
#                      oa = devs_oa,
#                      osa = devs_osa)

write_rds(out, here::here("data/temp", "perc_deviations.RDS"))
# write_rds(out_multicat, here::here("data/temp", "perc_deviations_multicat.RDS"))

tab01 <- bind_rows(
  Sex = out$sex, 
  Age = x$age, 
  Origin = x$orig, .id = "Category") %>% 
  select(
    Radius = radius, 
    Cluster = cluster, 
    Variable = Category, 
    Value = variable, 
    Observed = n, 
    Percent = larger_pct) %>% 
  arrange(Radius, Cluster, desc(Observed)) %>% 
  mutate(Percent = round(Percent, 1))

tab01 %>% 
  filter(Percent != 100) %>% 
  write_excel_csv(here::here("tables", "Table_01.csv"))
  
tab01 %>% 
  filter(Percent == 100) %>% 
  write_excel_csv(here::here("tables", "Table_02.csv"))

# exploration -------------------------------------

# x <- read_rds(here::here("data/temp", "perc_deviations.RDS"))
# 
# # explore clusters
# fuubar <- perc_cl %>% 
#   filter(radius == 5.6) %>% 
#   full_join(ved_sf) %>% 
#   st_as_sf
# 
# 
# ggplot() +
#   # geom_sf(data = ved_exc) +
#   geom_sf(data = fuubar, aes(color = factor(cluster))) +
#   # geom_sf(aes(color = factor(cluster))) +
#   ggrepel::geom_text_repel(
#     data = fuubar, 
#     aes(geometry = geometry, label = cluster), 
#     stat = "sf_coordinates") +
#   theme_void()



