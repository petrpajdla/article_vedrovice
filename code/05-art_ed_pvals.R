# Project "Vedrovice"
# Script nr. 
# 
# author: Petr Pajdla
# 

set.seed(42)

n_sim <- 9999

# libs --------------------------------------------------------------------

library(tidyverse)
library(here)
library(infer)

theme_set(theme_minimal())

# data --------------------------------------------------------------------

ved <- readRDS(here("data/temp", "vedrovice_dataset.RDS"))

# v_statistic <- read_csv(here("data/temp", "cooc_v_stats.csv"))
# v_normalized <- read_csv(here("data/temp", "cooc_v_normal.csv"))

non_rand_vars_abbrv <- read_lines(here("data/temp", "cooc_non_rand_vars.txt"))
non_rand_vars <- ved$var_names$full %>% 
  filter(abbrv %in% non_rand_vars_abbrv) %>% 
  pull(vnames)

# var_clusters <- read_csv(here("data/temp", "cooc_var_clusters.csv"))


# modified input ----------------------------------------------------------

# goal: count mean euclidean distance for groups of burials 
ved_input <- ved$metadata %>% 
  filter(pres != "dist.") %>% 
  select(id_burial, sex, age = age_sim, orig = origin, starts_with("cat_")) %>% 
  left_join(as_tibble(ved$bin_vars$bin_mat[, non_rand_vars], 
                      rownames = "id_burial"), by = c("id_burial"))


# functions ---------------------------------------------------------------

# euclidean distance between burials in groups
edistance <- function(df, grouping_var, distance_vars) {
  var <- sym(grouping_var)
  
  mat <- df %>% 
    select(id_burial, !!var, all_of(distance_vars)) %>% 
    group_by(!!var) %>% 
    nest() %>% 
    mutate(data = map(data, column_to_rownames, "id_burial"),
           data = map(data, as.matrix)) %>% 
    filter(!!var != "ind.")
  
  lvls <- as.character(pull(mat, !!var))
  comb <- gtools::combinations(n = length(lvls), 
                               r = 2, 
                               v = lvls, 
                               repeats.allowed = FALSE) %>% 
    as_tibble() %>% 
    mutate(name = str_c(V1, V2, sep = " - "))
  
  dist <- vector("list", nrow(comb))
  names(dist) <- comb$name
  
  for (i in 1:nrow(comb)) {
    dist[[i]] <- pdist::pdist(
      X = filter(mat, !!var == comb$V1[i])$data[[1]],
      Y = filter(mat, !!var == comb$V2[i])$data[[1]])@dist
  }
  
  map(dist, as_tibble) %>% 
    bind_rows(.id = "comb") %>% 
    group_by(comb) %>% 
    summarize(value = mean(value), .groups = "drop")
}

# randomization experiment
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
  }
  res %>% 
    # map(group_by, comb) %>% 
    # map(summarize, mean = mean(value)) %>% 
    bind_rows()
}

# # t test on random distributions
# test_randomness <- function(experimental, observed) {
#   bind_rows(exp = experimental, obs = observed, .id = "parent") %>% 
#     group_by(comb) %>% 
#     nest() %>% 
#     mutate(ttest = map(data, t_test, value ~ parent, order = c("exp", "obs")),
#            wtest = map2()) %>% 
#     select(-data) %>% 
#     unnest(cols = ttest) %>% 
#     ungroup() %>% 
#     mutate(signif = if_else(p_value <= 0.05, "*", ""),
#            signif = if_else(p_value <= 0.01, "**", signif),
#            signif = if_else(p_value <= 0.001, "***", signif))
# }
# 
# # combination function to deliver all results
# get_edistance <- function(df, grouping_var, distance_vars, n_sim) {
#   begin <- Sys.time()
#   res <- list(obs = NA,
#               exp = NA,
#               ttest = NA,
#               mean = NA)
#   
#   res$obs <-edistance(df, grouping_var, distance_vars)
#   
#   res$exp <- randomize_edistance(df, grouping_var, distance_vars, n_sim)
#   
#   res$ttest <- test_randomness(res$exp, res$obs)
#   
#   res$mean <- res$obs %>% group_by(comb) %>% 
#     summarise(mean = mean(value))
#   print(Sys.time() - begin)
#   return(res)
# }

get_p <- function(df, grouping_var, distance_vars, n_sim) {
  obs <- edistance(df, grouping_var, distance_vars)
  exp <- randomize_edistance(df, grouping_var, distance_vars, n_sim)
  
  res <- exp %>% group_by(comb) %>% 
    nest() %>% 
    full_join(obs, by = "comb") %>% 
    mutate(obs = map(value, as_tibble),
           larger = map2(data, obs, function(x, y) x$value > y$value),
           pval = map(larger, mean)) %>% 
    select(comb, value, pval) %>% 
    ungroup() %>% 
    unnest(pval)
  return(res)
}

# analysis ----------------------------------------------------------------

ed_sex <- get_p(ved_input, "sex", non_rand_vars, n_sim)
ed_age <- get_p(ved_input, "age", non_rand_vars, n_sim)
ed_orig <- get_p(ved_input, "orig", non_rand_vars, n_sim)
ed_sa <- get_p(ved_input, "cat_sa", non_rand_vars, n_sim)
ed_os <- get_p(ved_input, "cat_os", non_rand_vars, n_sim)
ed_oa <- get_p(ved_input, "cat_oa", non_rand_vars, n_sim)
ed_osa <- get_p(ved_input, "cat_osa", non_rand_vars, n_sim)

ed_sex %>% mutate(signif = if_else(pval <= 0.05, "*", ""),
                  signif = if_else(pval <= 0.01, "**", signif),
                  signif = if_else(pval <= 0.001, "***", signif))
ed_age %>% mutate(signif = if_else(pval <= 0.05, "*", ""),
                  signif = if_else(pval <= 0.01, "**", signif),
                  signif = if_else(pval <= 0.001, "***", signif))
ed_orig %>% mutate(signif = if_else(pval <= 0.05, "*", ""),
                   signif = if_else(pval <= 0.01, "**", signif),
                   signif = if_else(pval <= 0.001, "***", signif))
ed_sa %>% mutate(signif = if_else(pval <= 0.05, "*", ""),
                  signif = if_else(pval <= 0.01, "**", signif),
                  signif = if_else(pval <= 0.001, "***", signif))
ed_os %>% mutate(signif = if_else(pval <= 0.05, "*", ""),
                 signif = if_else(pval <= 0.01, "**", signif),
                 signif = if_else(pval <= 0.001, "***", signif))
ed_oa %>% mutate(signif = if_else(pval <= 0.05, "*", ""),
                 signif = if_else(pval <= 0.01, "**", signif),
                 signif = if_else(pval <= 0.001, "***", signif))
ed_osa %>% mutate(signif = if_else(pval <= 0.05, "*", ""),
                  signif = if_else(pval <= 0.01, "**", signif),
                  signif = if_else(pval <= 0.001, "***", signif))

# output ------------------------------------------------------------------

output_pvals <- list(sex = ed_sex,
                     age = ed_age,
                     orig = ed_orig,
                     sa = ed_sa,
                     os = ed_os,
                     oa = ed_oa,
                     osa = ed_osa)

write_rds(output_pvals, here("data/temp", "cooc_ed.RDS"))

ed_pvals <- read_rds(here("data/temp", "cooc_ed.RDS"))

ggplot(exp, aes(value)) +
  geom_density() +
  geom_vline(data = obs, aes(xintercept = value)) +
  facet_wrap(~comb)


# plots -------------------------------------------------------------------

plot_densities <- function(ed_obj) {
  ed_obj$obs %>% 
    ggplot(aes(value)) +
    geom_density(data = ed_obj$exp, color = "gray") +
    geom_rug(data = ed_obj$exp, color = "gray", alpha = 0.2) +
    geom_density() +
    geom_rug(alpha = 0.2) +
    geom_vline(data = ed_obj$mean, aes(xintercept = mean)) +
    facet_wrap(vars(comb))
}

plot_densities(ed_orig)

z %>% 
  ggplot(aes(mean)) +
  geom_density() +
  facet_wrap(~comb) +
  geom_vline(data = y, aes(xintercept = mean))



# data centric p ----------------------------------------------------------

eeeh <- function(observed, expected) {
  obs <- observed %>%
    group_by(comb) %>%
    nest() %>%
    rename(obs = data)
  expected %>%
    group_by(comb) %>%
    nest() %>%
    rename(exp = data) %>%
    full_join(obs, by = c("comb")) %>%
    mutate(across(c(exp, obs), map, as.matrix),
           # wtest = map2(obs, exp, wilcox.test),
           # wtest = map(wtest, broom::tidy),
           kstest = map2(obs, exp, ks.test),
           kstest = map(kstest, broom::tidy)
    ) %>%
    select(-exp, -obs) %>%
    ungroup() %>%
    unnest(wtest) %>% 
    mutate(signif = if_else(p.value <= 0.05, "*", ""),
           signif = if_else(p.value <= 0.01, "**", signif),
           signif = if_else(p.value <= 0.001, "***", signif))
}
