# Project "Vedrovice"
# Script nr. 2.1
# ANALYSIS OF ARTEFACT CO-OCCURRENCES RANDOMNESS
# author: Petr Pajdla
# Randomization test on the input binary matrix is performed in order to 
# find structure in artefact coocurences
# generalized Monte Carlo procedure...

set.seed(42)

# packages =====================================================================
library(here)
library(dplyr)
library(purrr)
library(igraph)
library(ggplot2)

# theme for ggplot graphics ====================================================
theme_universe <- theme(panel.border = element_rect(colour = "black", 
                                                    fill = NA, 
                                                    size = 0.8),
                        panel.background = element_blank(),
                        line = element_blank(),
                        strip.background = element_blank(), 
                        # strip.text = element_text(face = "italic"),
                        axis.text.y = element_blank(), 
                        axis.title.y = element_blank(), 
                        panel.spacing = unit(1, "lines"))

# functions ====================================================================
# count coocurences in an occurence (binary) matrix
# @var matrix = binary matrix of co-occurrences
count_cooccurrence <- function(matrix) {
  colpairs <- expand.grid(colnames(matrix), colnames(matrix))
  result <- data.frame(colpairs[, 1:2],
                       rep(0, nrow(colpairs)))
  names(result) <- c("var1", "var2", "nr.cooc")
  for (i in 1:nrow(colpairs)) {
    result[i, "nr.cooc"] <- sum(unname(matrix[, unname(colpairs[i, 1])] == 1) & 
                                  unname(matrix[, unname(colpairs[i, 2])] == 1))
  }
  return(result)
}

# count v values for a data frame (based on articles by Manly)
# where: R...nr of variables, df...data frame of coocurences, 
# expected = df of expected  coocurences (med.coocs)
# @var df = dataframe of sums of co-occurences
# @var expected = data frame of expected co-occurences
# @var R = number of variables
get_v <- function(df, expected, R) {
  v_vals <- rep(NA, R)
  temp_df <- left_join(df, expected, by = c("var1", "var2"))
  v_vals <- temp_df %>% group_by(var1) %>% 
    summarise(v = sum((nr.cooc - mean.cooc)^2 / R))
  return(v_vals)
}

# get p value (percentage of values of v statistic as large or larger than
# observed value of v statistic in a randomly changed occurrence matrices)
# @var v_exp = expected values of statistic v (data frame)
# @var v_obs = observed v values (data frame)
# @var size = number of v values derived by randomization
get_p <- function(v_exp, v_obs) {
  mean(v_exp >= v_obs)
}

# normalize to 0-1 range
normalize01 <- function(x) {
  rg <- range(x, na.rm = TRUE)
  (x - rg[1]) / (rg[2] - rg[1])
}

# read data ====================================================================
ved <- readRDS(here("data/temp", "vedrovice_dataset.RDS"))

# counting co-occurrences in observed matrix ===================================
input_matrix <- ved$bin_vars$bin_mat[, ved$bin_vars$over5]
ncol_input <- ncol(input_matrix)

cooc_obs <- count_cooccurrence(matrix = input_matrix)
cooc_obs <- cooc_obs %>% filter(var1 != var2)
# cooc_obs %>% tidyr::spread(var2, nr.cooc)

# randomization of co-occurrence matrix - list of many matrices ----------------
n_permutations <-  9999
# ======
rand_mat <- vegan::permatfull(input_matrix,
                              fixedmar = "both",
                              mtype = "prab",
                              times = n_permutations)
# ======

# going parallel to speed up randomization...
# detecting cores for parallel
no_cores <- parallel::detectCores() - 1

# cooccurences on random matrices ==============================================
# ======
# note to self: on Debian without everything was smooth even without specifying 
# setup_strategy, but apparently, this is Rstudio vs R version problem
# see: https://github.com/rstudio/rstudio/issues/6692
cl <- parallel::makeCluster(no_cores, setup_strategy = "sequential")
cooc_rand <- parallel::parLapply(cl, rand_mat$perm, count_cooccurrence)
parallel::stopCluster(cl)

# save counts of co-occurences on random matrices
readr::write_rds(cooc_rand, here("data/temp", "cooc_random_mat.RDS"))
# ======

# load from temporary data
cooc_rand <- readr::read_rds(here("data/temp", "cooc_random_mat.RDS"))

# summarization of random occurences into expected cooccurences ----------------
cooc_exp <- cooc_rand %>% bind_rows() %>% 
  group_by(var1, var2) %>% 
  summarise(mean.cooc = mean(nr.cooc),
            med.cooc = median(nr.cooc))
# cooc_exp <- cooc_exp  %>% filter(var1 != var2)

# v values =====================================================================
# observed v values ------------------------------------------------------------
v_observed <- tidyr::spread(get_v(cooc_obs, cooc_exp, 
                                  R = ncol_input),
                            var1, v)

# counting expected v values ---------------------------------------------------
# creating cluster with forked global environment
cl <- parallel::makeCluster(no_cores, type = "FORK")
v_experimental <- parallel::parLapply(cl, cooc_rand,
                                      get_v, cooc_exp, 
                                      ncol_input)
v_experimental <- bind_rows(parallel::parLapply(cl, v_experimental, 
                                                tidyr::spread, "var1", "v"))
parallel::stopCluster(cl)

# P-val: percentage of experimental v values that are larger than observed v ---
v_p_values <- map2(v_experimental, v_observed, get_p) %>% 
  unlist()

v_statistic <- tibble(variable = names(v_observed),
                      v = t(unname(v_observed))[, 1],
                      p = unname(v_p_values)) %>% 
  mutate(
    # across(is.numeric, round, 2),
    signif = if_else(p <= 0.05, "*", ""),
    signif = if_else(p <= 0.01, "**", signif),
    signif = if_else(p <= 0.001, "***", signif)) %>% 
  mutate(abbrv = unname(ved$var_names$short[variable]),
         long = unname(ved$var_names$long[variable]))

readr::write_csv(v_statistic, here("data/temp", "cooc_v_stats.csv"))
# v_statistic <- readr::read_csv(here("data/temp", "cooc_v_stats.csv"))

# visualizing of v values by small multiples -----------------------------------
v_exp_g <- v_experimental %>% tidyr::gather(variable, value) %>% 
  mutate(abbrv = forcats::as_factor(unname(ved$var_names$short[variable])),
         long = forcats::as_factor(unname(ved$var_names$long[variable])))
v_obs_g <- v_observed %>% tidyr::gather(variable, value) %>% 
  mutate(abbrv = forcats::as_factor(unname(ved$var_names$short[variable])),
         long = forcats::as_factor(unname(ved$var_names$long[variable])))

v_annot <- v_statistic %>% select(long, p, signif) %>%
  mutate(p = round(p, 2),
         long = forcats::as_factor(long),
         txt = paste0(p, signif))

v_plot <- ggplot(v_exp_g, mapping = aes(x = value)) +
  geom_density(fill = "gray80", color = NA, alpha = 0.6) +
  geom_rug(alpha = 0.2) +
  geom_vline(data = v_obs_g, mapping = aes(xintercept = value), size = 0.8) +
  facet_wrap(vars(long), scales = "free", ncol = 3) +
  xlab("v statistic") +
  geom_text(data = v_annot, aes(x = Inf, y = Inf, label = txt), size = 3,
            hjust = +1.2, vjust = +1.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0, 0, 0.2)) +
  theme_universe

v_plot

ggsave(here("plots", "cooc_vstats.pdf"), 
       plot = v_plot, 
       width = 19, height = 20, units = "cm")


# normalizing to range 0-1 ------------------------------------------------
v_normalized <- bind_rows(exp = v_exp_g, obs = v_obs_g, .id = "orig") %>% 
  group_by(variable) %>% 
  tidyr::nest() %>% 
  mutate(data = map(data, mutate, value_norm = normalize01(value))) %>% 
  tidyr::unnest(cols = data)

v_normalized_plot <- ggplot(filter(v_normalized, orig == "exp"), 
                            mapping = aes(x = value_norm)) +
  geom_density(fill = "gray80", color = NA, alpha = 0.6) +
  geom_rug(alpha = 0.2) +
  geom_vline(data = filter(v_normalized, orig == "obs"), 
             mapping = aes(xintercept = value_norm), size = 0.8) +
  facet_wrap(vars(long), scales = "free", ncol = 3) +
  xlab("v statistic (normalized to 0 - 1 range)") +
  geom_text(data = v_annot, aes(x = Inf, y = Inf, label = txt), size = 3,
            hjust = +1.2, vjust = +1.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0, 0, 0.2)) +
  theme_universe

# ggsave(here("plots", "cooc_vstats_norm.pdf"), plot = v_normalized_plot, 
#        width = 6, height = 6)

readr::write_csv(v_normalized, here("data/temp", "cooc_v_normal.csv"))


# S statistic ==================================================================
s_observed <- sum(v_observed / ncol_input)

s_experimental <- rowSums(apply(v_experimental, 2, 
                                "/", ncol_input))

s_p_value <- get_p(s_experimental, s_observed)
s_p_value

s_plot <- ggplot(as_tibble(s_experimental), aes(value)) +
  geom_density(fill = "gray80", color = NA, alpha = 0.6) +
  geom_rug(alpha = 0.2) +
  geom_vline(xintercept = s_observed, size = 0.8) +
  geom_text(data = tibble(p = paste0(round(s_p_value, 3), "**")), 
            aes(x = Inf, y = Inf, label = p), 
            size = 3,
            hjust = +1.1, vjust = +1.4) +
  xlab("S statistic") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0, 0, 0.2)) +
  theme_universe

s_plot

ggsave(here("plots", "cooc_sstat.pdf"), plot = s_plot, 
       width = 6, height = 4, units = "cm")

# The S statistic is significantly larger than experimental S and
# occures in less then 5% of cases (p is smaller than 0.05), i.e. there is
# strong evidence for non-random artefact distribution

# quazi corrplot to show distribution / co-occurrence of artefacts =============
# observed co-occurence - expected co-occurrence (mean) = dif
# if (dif > 0) [blue] {
#       overrepresentation in observed data
#   } else i <- theme(panel.border = element_rect(colour = "black",
# fill = NA, 
# size = 0.8),
# panel.background = element_blank(),
# line = element_blank(),
# strip.background = element_blank(), 
# # strip.text = element_text(face = "italic"),
# axis.text.y = element_blank(), 
# axis.title.y = element_blank(), 
# f (dif < 0) [red] {
# #       underrepresentaion in observed data
#   }

cooccurrence <- left_join(cooc_obs, cooc_exp, by = c("var1", "var2")) %>% 
  as_tibble() %>% 
  mutate(dif = nr.cooc - mean.cooc,
         var1 = unname(ved$var_names$short[as.character(var1)]),
         var2 = unname(ved$var_names$short[as.character(var2)])) %>% 
  select(var1, var2, dif) %>% 
  reshape2::acast(var1 ~ var2, value.var = "dif")
cooccurrence[is.na(cooccurrence)] <- 0

# bw <- colorRampPalette(colors = c("gray90", "gray20"))

# pdf(here("plots", "cooc_corrplot.pdf"), width = 7, height = 7)
# corrplot::corrplot(round(cooccurrence, 2), 
#                    is.corr = FALSE, method = "pie", 
#                    type = "upper", diag = FALSE, 
#                    order = "FPC", 
#                    tl.col = "gray20", tl.cex = 0.8, 
#                    cl.cex = 0.6, 
#                    # col = bw(12)
# )
# dev.off()

# pdf(here("plots", "cooccurrence_alphabet.pdf"), width = 7, height = 7)
# corrplot::corrplot(cooccurrence, 
#                    is.corr = FALSE, method = "pie", 
#                    type = "upper", diag = FALSE, 
#                    order = "alphabet", 
#                    tl.col = "gray20", tl.cex = 0.8, 
#                    cl.cex = 0.6)
# dev.off()

# co-occurrence networks for:
#   1) non random variables only
#   2) all variables *
# non random variables =========================================================
# non_rand_vars <- v_statistic %>% filter(signif == TRUE) %>% pull(abbrv)
# non_rand_cooc <- list(pos = cooccurrence[non_rand_vars, non_rand_vars],
#                       neg = cooccurrence[non_rand_vars, non_rand_vars])
# non_rand_cooc$pos[non_rand_cooc$pos < 0] <- 0
# non_rand_cooc$neg[non_rand_cooc$neg > 0] <- 0
# non_rand_cooc$neg <- abs(non_rand_cooc$neg)
# 
# g_cooc_positive <- simplify(graph_from_adjacency_matrix(non_rand_cooc$pos, 
#                                                         mode = "undirected", 
#                                                         weighted = TRUE))
# g_cooc_negative <- simplify(graph_from_adjacency_matrix(non_rand_cooc$neg, 
#                                                         mode = "undirected",
#                                                         weighted = TRUE))

# create co-occurrence / co-absence matrix
cooc_posneg <- list(pos = cooccurrence,
                    neg = cooccurrence)
cooc_posneg$pos[cooccurrence < 0] <- 0
cooc_posneg$neg[cooccurrence > 0] <- 0
cooc_posneg$neg <- cooc_posneg$neg * -1

# filter for co-occurrences / co-absences over 1
cooc_posneg$pos[cooc_posneg$pos <= 1] <- 0
cooc_posneg$neg[cooc_posneg$neg <= 1] <- 0

# create networks, simplify edges, delete vertices with 0 edges
g_cooc_positive <- simplify(graph_from_adjacency_matrix(cooc_posneg$pos,
                                                        mode = "undirected",
                                                        weighted = TRUE))
g_cooc_positive <- delete.vertices(g_cooc_positive, 
                                   degree(g_cooc_positive) == 0)

g_cooc_negative <- simplify(graph_from_adjacency_matrix(cooc_posneg$neg,
                                                        mode = "undirected",
                                                        weighted = TRUE))
g_cooc_negative <- delete.vertices(g_cooc_negative, 
                                   degree(g_cooc_negative) == 0)

non_rand_vars <- v_statistic %>% filter(stringr::str_detect(signif, "\\*")) %>% pull(abbrv)

pdf(here("plots", "cooc_nets.pdf"), width = 12)
par(mfrow = c(1, 2), mar = c(rep(2, 4)))
plot(g_cooc_positive, 
     vertex.shape = "circle", 
     vertex.color = "white",
     # vertex.size = eigen_centrality(g_cooc_positive)$vector*20,
     vertex.label.family = "sans", 
     vertex.label.cex = .6,
     vertex.label.font = if_else(vertex_attr(g_cooc_positive)$name %in% non_rand_vars, 2, 3),
     vertex.label.color = "black",
     edge.width = edge_attr(g_cooc_positive)$weight^1.2,
     edge.color = "black",
     mark.groups = cluster_fast_greedy(g_cooc_positive),
     mark.col = "gray90",
     mark.border = "white")
title(main = "present", cex.main = 1, family = "sans", font.main = 1)
plot(g_cooc_negative, 
     vertex.shape = "circle", 
     vertex.color = "white",
     # vertex.size = eigen_centrality(g_cooc_negative)$vector*20,
     vertex.label.family = "sans", 
     vertex.label.cex = .6,
     vertex.label.font = if_else(vertex_attr(g_cooc_negative)$name %in% non_rand_vars, 2, 3),
     vertex.label.color = "black",
     edge.width = edge_attr(g_cooc_negative)$weight^1.2,
     edge.color = "black",
     mark.groups = cluster_fast_greedy(g_cooc_negative),
     mark.col = "gray90",
     mark.border = "white")
title(main = "absent", cex.main = 1, family = "sans", font.main = 1)
dev.off()

# clusters of variables
cl_absent <- cluster_fast_greedy(g_cooc_negative) %>% 
  groups()
cl_present <- cluster_fast_greedy(g_cooc_positive) %>% 
  groups()

extract_network_groups <- function(lst, type) {
  len <- length(lst)
  temp <- vector("list", len)
  for (i in 1:len) {
    temp[[i]] <- tibble(cluster = paste0(type, i), abbrv = lst[[i]])
  }
  bind_rows(temp)
}

var_clusters <- bind_rows(extract_network_groups(cl_absent, "absent"), 
                          extract_network_groups(cl_present, "present"))

# output =======================================================================
writeLines(non_rand_vars, here("data/temp", "cooc_non_rand_vars.txt"))
readr::write_csv(var_clusters, here("data/temp", "cooc_var_clusters.csv"))
