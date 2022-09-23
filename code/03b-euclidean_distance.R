# Project "Vedrovice"
# Script nr. 3B
# ED between pairs of burials
# 
# author: Petr Pajdla

set.seed(42)

n_sim <- 9999

# libs --------------------------------------------------------------------

library(tidyverse)
library(infer)

theme_set(theme_minimal())

# data --------------------------------------------------------------------

ved <- readRDS(here::here("data/temp", "vedrovice_dataset.RDS"))

# v_statistic <- read_csv(here("data/temp", "cooc_v_stats.csv"))
# v_normalized <- read_csv(here("data/temp", "cooc_v_normal.csv"))

non_rand_vars_abbrv <- read_lines(here::here("data/temp", "cooc_non_rand_vars.txt"))
non_rand_vars <- ved$var_names$full %>% 
  filter(abbrv %in% non_rand_vars_abbrv) %>% 
  pull(vnames)

# var_clusters <- read_csv(here("data/temp", "cooc_var_clusters.csv"))


# modified input ----------------------------------------------------------

# goal: count mean euclidean distance for groups of burials 
ved_input <- ved$metadata %>% 
  filter(analysis) %>% 
  select(id_burial, sex, age = age_sim, orig = origin) %>% 
  left_join(as_tibble(ved$bin_vars$count_mat[, non_rand_vars], 
                      rownames = "id_burial"), by = c("id_burial"))


# functions ---------------------------------------------------------------

#' Euclidean distance
#' 
#' Get Euclidean distance between pairs of burials.
#'
#' @param df A data frame with a grouping variable column and 
#' columns defining the space.
#' @param grouping_var A grouping variable, factor.
#' @param distance_vars Numeric variables defining the euclidean space.
#'
#' @return Melted matrix of Euclidean distances between 
#' @export
#'
#' @examples
edistance <- function(df, grouping_var, distance_vars) {
  var <- sym(grouping_var)
  
  mat <- df %>%
    select(id_burial, !!var, all_of(distance_vars)) %>%
    group_by(!!var) %>%
    nest() %>%
    mutate(data = map(data, column_to_rownames, "id_burial"),
           data = map(data, as.matrix))
  
  lvls <- as.character(pull(mat, !!var))
  comb <- gtools::combinations(n = length(lvls),
                               r = 2,
                               v = lvls,
                               repeats.allowed = TRUE) %>%
    as_tibble() %>%
    mutate(name = str_c(V1, V2, sep = " - "))
  
  dist <- vector("list", nrow(comb)) %>% 
    setNames(comb$name)
  
  for (i in 1:nrow(comb)) {
    x <- comb$V1[i]
    y <- comb$V2[i]
    if (x == y) {
      dist[[i]] <- as.matrix(dist(filter(mat, !!var == x)$data[[1]]))
    } else {
      dist[[i]] <- as.matrix(pdist::pdist(
        X = filter(mat, !!var == x)$data[[1]],
        Y = filter(mat, !!var == y)$data[[1]]))
    }
  }
  
  map(dist, mean) %>%
    as_tibble() %>%
    pivot_longer(everything())
  # filter(value != 0)
}

#' Randomizate Euclidean distances
#'
#' @param df A data frame.
#' @param grouping_var Grouping variable.
#' @param distance_vars Distance variables.
#' @param n_sim Number of simulations.
#'
#' @return A data frame with mean group Euclidean distances 
#' for each of the iterations.
#' @export
#'
#' @examples
randomize_edistance <- function(df, grouping_var, distance_vars, n_sim) {
  n_bur <- nrow(df)
  var <- as.symbol(grouping_var)
  
  prob <- df %>%
    select(!!var) %>%
    group_by(!!var) %>%
    count() %>%
    ungroup(!!var) %>% 
    mutate(prob = n / n_bur)
  
  lvls <- prob %>% 
    pull(!!var) %>% 
    levels()
  
  res <- vector("list", n_sim)
  
  for (i in 1:n_sim) {
    rand_vect <- randomizr::simple_ra(n_bur, prob_each = prob$prob, 
                                      conditions = pull(prob, !!var)) %>% 
      as.character()
    
    input <- select(df, id_burial, all_of(distance_vars)) %>% 
      bind_cols(temp = rand_vect)
    
    res[[i]] <- edistance(input, "temp", distance_vars)
    # group_by(comb) %>% 
    # summarise(mean = mean(value))
  }
  res %>% 
    # map(group_by, comb) %>% 
    # map(summarize, mean = mean(value)) %>% 
    bind_rows(.id = "iter")
}

#' Test proability of observed ED
#'
#' @param experimental Monte-Carlo simulated values of ED.
#' @param observed The observed values of ED.
#'
#' @return
#' @export
#'
#' @examples
test_randomness <- function(experimental, observed) {
  exp <- experimental %>% 
    select(name, value) %>% 
    group_by(name) %>% 
    nest()
  
  exp %>% 
    left_join(observed, by = "name") %>% 
    mutate(larger = map2(data, value, \(x, y) filter(x, value >= y)),
           smaller = map2(data, value, \(x, y) filter(x, value <= y)),
           larger_pct = map2_dbl(data, larger, \(x, y) (nrow(y) / nrow(x))*100),
           smaller_pct = map2_dbl(data, smaller, \(x, y) (nrow(y) / nrow(x))*100))
  # filter(larger_pct < 5 | smaller_pct < 5)
}

#' Master function to deliver all results
#'
#' @param df 
#' @param grouping_var 
#' @param distance_vars 
#' @param n_sim 
#'
#' @return
#' @export
#'
#' @examples
get_edistance <- function(df, grouping_var, distance_vars, n_sim) {
  begin <- Sys.time()
  res <- list(obs = NA,
              exp = NA,
              test = NA)
  res$obs <- edistance(df, grouping_var, distance_vars)
  res$exp <- randomize_edistance(df, grouping_var, distance_vars, n_sim)
  res$test <- test_randomness(res$exp, res$obs)
  
  # res$mean <- res$obs %>% group_by(comb) %>% 
  #   summarise(mean = mean(value))
  print(Sys.time() - begin)
  return(res)
}


# analysis ----------------------------------------------------------------

# edistance(ved_input, "sex", non_rand_vars)
# randomize_edistance(ved_input, "sex", non_rand_vars, 10)
# test_randomness(randomize_edistance(ved_input, "sex", non_rand_vars, 10), 
#                 edistance(ved_input, "sex", non_rand_vars))

ed_sex <- get_edistance(ved_input, "sex", non_rand_vars, n_sim)
ed_age <- get_edistance(ved_input, "age", non_rand_vars, n_sim)
ed_orig <- get_edistance(ved_input, "orig", non_rand_vars, n_sim)


# output ------------------------------------------------------------------

ed <- list(sex = ed_sex,
           age = ed_age,
           orig = ed_orig)

# write_rds(ed, here::here("data/temp", "cooc_ed.RDS"))
# ed <- read_rds(here("data/temp", "cooc_ed.RDS"))


# output table ------------------------------------

# table 00

# bind_rows(
#   sex = select(ed$sex$test, name, value, ends_with("pct")),
#   age = select(ed$age$test, name, value, ends_with("pct")),
#   orig = select(ed$orig$test, name, value, ends_with("pct")), 
#   .id = "var"
#   ) %>% 
#   write_csv(here::here("tables/Table_00.csv"))
  

# plots -------------------------------------------------------------------

# ggplot() +
#   geom_density(data = rand_sex, aes(x = m)) +
#   geom_vline(data = ed_sex, aes(xintercept = mean(value))) +
#   facet_wrap(vars(comb))

# plot_densities <- function(ed_obj) {
#   ed_obj$exp %>%
#     ggplot(aes(value)) +
#     geom_density() +
#     # geom_rug(color = "gray", alpha = 0.2) +
#     geom_vline(data = ed_obj$obs, aes(xintercept = value)) +
#     facet_wrap(vars(name), scales = "free")
# }
# 
plot_densities <- function(ed_obj) {
  ed_obj$exp %>%
    ggplot(aes(value, name)) +
    geom_boxplot(color = "gray",
                 width = 0.1) +
    geom_violin(fill = NA) +
    geom_point(data = ed_obj$obs,
               shape = 21,
               size = 2.4, fill = "black")
  # facet_wrap(vars(name), scales = "free")
}

plot_densities(ed_sex)
plot_densities(ed_age)
plot_densities(ed_orig)

# z %>% 
#   ggplot(aes(mean)) +
#   geom_density() +
#   facet_wrap(~comb) +
#   geom_vline(data = y, aes(xintercept = mean))
