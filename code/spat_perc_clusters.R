# Project "Vedrovice"
# Script nr. 
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

explored_rad <- paste0(c(2.8, 3.9, 4.6, 5.5, 6.5, 7.8))

perc_cl <- read_csv(here::here("data/temp", "perc_clusters.csv")) %>% 
  mutate(id_burial = as.character(id)) %>% 
  select(-id) %>% 
  filter(radius %in% explored_rad)

ved_sf <- sf::st_read(here::here("data/temp", "layout.shp")) %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")),
         age_sim = fct_relevel(age_sim, c("juv.", "ad.", "mat.", "ind.")),
         origin = fct_relevel(origin, c("local", "non-local", "ind."))) %>% 
  filter(pres != "dist.")

# prep ---------------------------------------------------------------------

input <- ved_sf %>% select(id_burial, sex, age_sim, origin) %>% 
  st_drop_geometry() %>% 
  left_join(perc_cl, by = c("id_burial"))


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
  
  for (i in 1:n_sim) {
    rand_vect <- randomizr::simple_ra(n_bur, prob_each = prob$prob, 
                                      conditions = pull(prob, !!variable0)) %>% 
      as.character()
    
    input <- input %>% 
      select(id_burial, radius, cluster) %>% 
      mutate(temp = rand_vect)
    
    res[[i]] <- count_var_per_radius(input, "temp")
  }
  res %>% 
    bind_rows() %>% 
    group_by(radius, cluster, variable) %>% 
    summarise(n = mean(n, na.rm = TRUE), .groups = "drop") %>% 
    ungroup()
}

return_deviations_from_random <- function(expected, observed) {
  bind_rows(expected = expected, observed = observed, .id = "origin") %>% 
    group_by(radius) %>% 
    nest() %>%
    mutate(data = map(data, pivot_wider, names_from = origin, values_from = n),
           data = map(data, mutate, observed = if_else(is.na(observed), 0, observed)),
           data = map(data, mutate, deviation = observed - expected),
           data = map(data, group_by, variable),
           sums = map(data, summarise, 
                      m = mean(deviation),
                      sd = sd(deviation),
                      max = m + 2*sd,
                      min = m - 2*sd),
           full = map2(data, sums, left_join, by = "variable")) %>%
    select(full) %>% 
    unnest(cols = full) %>% 
    ungroup(radius) %>% 
    mutate(flt = (deviation >= max)|(deviation <= min)) %>% 
    filter(flt) %>% 
    select(radius, cluster, variable, expected, observed, deviation) %>% 
    arrange(desc(deviation))
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

devs_sex <- inner_join(deviates_from_random(input, "sex", nsim),
                       deviates_from_random(input, "sex", nsim), 
                       by = c("radius", "cluster", "variable")) %>% 
  transmute(radius, cluster, variable,
            observed = observed.x,
            expected = (expected.x + expected.y)/2,
            deviation = (deviation.x + deviation.y)/2) %>% 
  inner_join(deviates_from_random(input, "sex", nsim), 
             by = c("radius", "cluster", "variable")) %>% 
  transmute(radius, cluster, variable,
            observed = observed.x,
            expected = (expected.x + expected.y)/2,
            deviation = (deviation.x + deviation.y)/2)

devs_age <- inner_join(deviates_from_random(input, "age_sim", nsim),
                       deviates_from_random(input, "age_sim", nsim), 
                       by = c("radius", "cluster", "variable")) %>% 
  transmute(radius, cluster, variable,
            observed = observed.x,
            expected = (expected.x + expected.y)/2,
            deviation = (deviation.x + deviation.y)/2) %>% 
  inner_join(deviates_from_random(input, "age_sim", nsim),
             by = c("radius", "cluster", "variable")) %>% 
  transmute(radius, cluster, variable,
            observed = observed.x,
            expected = (expected.x + expected.y)/2,
            deviation = (deviation.x + deviation.y)/2)

devs_loc <- inner_join(deviates_from_random(input, "origin", nsim),
                       deviates_from_random(input, "origin", nsim), 
                       by = c("radius", "cluster", "variable")) %>% 
  transmute(radius, cluster, variable,
            observed = observed.x,
            expected = (expected.x + expected.y)/2,
            deviation = (deviation.x + deviation.y)/2) %>% 
  inner_join(deviates_from_random(input, "origin", nsim),
             by = c("radius", "cluster", "variable")) %>% 
  transmute(radius, cluster, variable,
            observed = observed.x,
            expected = (expected.x + expected.y)/2,
            deviation = (deviation.x + deviation.y)/2)


# output ------------------------------------------------------------------

out <- list(sex = devs_sex,
            age = devs_age,
            orig = devs_loc)

write_rds(out, here::here("data/temp", "perc_deviations.RDS"))

