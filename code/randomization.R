# Project "Vedrovice"
# Script nr. 2
# ANALYSIS OF ARTEFACT CO-OCCURRENCES 
# author: Petr Pajdla
# Randomization test on the input binary matrix is performed in order to 
# find structure in artefact coocurences
# generalized Monte Carlo procedure...

set.seed(42)

# packages =====================================================================
library(here)
library(dplyr)
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
get_p <- function(v_exp, v_obs, size) {
  res <- (sum(v_exp >= v_obs)/size)
  return(res)
}

# read data ====================================================================
ved <- readRDS(here("data/temp", "vedrovice_dataset.RDS"))

# counting co-occurrences in observed matrix ===================================
cooc_obs <- count_cooccurrence(matrix = ved$bin_vars$bin_mat)
cooc_obs <- cooc_obs %>% filter(var1 != var2)
# cooc_obs %>% tidyr::spread(var2, nr.cooc)

# randomization of co-occurrence matrix - list of many matrices ----------------
n_permutations <-  999
# # ======
# rand_mat <- vegan::permatfull(ved$bin_vars$bin_mat,
#                               fixedmar = "both",
#                               mtype = "prab",
#                               times = n_permutations)
# # ======

# going parallel to speed up randomization...
# detecting cores for parallel
no_cores <- parallel::detectCores() - 1

# cooccurences on random matrices ==============================================
# # ======
# cl <- parallel::makeCluster(no_cores)
# cooc_rand <- parallel::parLapply(cl, rand_mat$perm, count_cooccurrence)
# parallel::stopCluster(cl)
# 
# # save counts of co-occurences on random matrices
# readr::write_rds(cooc_rand, here("data/temp", "cooc_random_mat.RDS"))
# # ======

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
                                  R = length(colnames(ved$bin_vars$bin_mat))),
                            var1, v)

# counting expected v values ---------------------------------------------------
# creating cluster with forked global environment
cl <- parallel::makeCluster(no_cores, type = "FORK")
v_experimental <- parallel::parLapply(cl, cooc_rand,
                                      get_v, cooc_exp, 
                                      length(colnames(ved$bin_vars$bin_mat)))
v_experimental <- bind_rows(parallel::parLapply(cl, v_experimental, 
                                                tidyr::spread, "var1", "v"))
parallel::stopCluster(cl)

# P-val: percentage of experimental v values that are larger than observed v ---
v_p_values <- purrr::map2(v_experimental, v_observed, get_p, n_permutations) %>% unlist()

v_statistic <- tibble(variable = names(v_observed),
                      v = t(unname(v_observed))[, 1],
                      p = unname(v_p_values),
                      signif = unname(v_p_values < (5 / length(v_p_values)))) %>% 
  mutate(abbrv = unname(ved$var_names$short[variable]),
         long = unname(ved$var_names$long[variable]))

readr::write_csv(v_statistic, here("data/temp", "v_statistics.csv"))

# visualizing of v values by small multiples -----------------------------------
v_exp_g <- v_experimental %>% tidyr::gather(variable, value) %>% 
  mutate(abbrv = forcats::as_factor(unname(ved$var_names$short[variable])),
         long = forcats::as_factor(unname(ved$var_names$long[variable])))
v_obs_g <- v_observed %>% tidyr::gather(variable, value) %>% 
  mutate(abbrv = forcats::as_factor(unname(ved$var_names$short[variable])),
         long = forcats::as_factor(unname(ved$var_names$long[variable])))

v_annot <- v_statistic %>% select(long, p, signif) %>% 
  mutate(signif = if_else(signif, "*", "")) %>% 
  mutate(p = paste0(round(p, 2)*100, "% ", signif))

ggplot(v_exp_g, mapping = aes(x = value)) +
  geom_density(fill = "gray80", color = NA, alpha = 0.6) +
  geom_rug() +
  geom_vline(data = v_obs_g, mapping = aes(xintercept = value), size = 0.8) +
  facet_wrap(vars(long), scales = "free", nrow = 4) +
  xlab("v statistic") +
  # geom_text(data = v_annot, aes(x = Inf, y = Inf, label = p), size = 3.4,
  #           hjust = +1.1, vjust = +1.2) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_universe

ggsave(here("plots", "v_values.pdf"), width = 12, height = 8)

# S statistic ==================================================================
s_observed <- sum(v_observed / length(names(ved$bin_vars)))

s_experimental <- rowSums(apply(v_experimental, 2, 
                                "/", length(names(ved$bin_vars))))

ggplot(as_tibble(s_experimental), aes(value)) +
  geom_density(fill = "gray80", color = NA, alpha = 0.6) +
  geom_rug() +
  geom_vline(xintercept = s_observed, size = 0.8) +
  xlab("S statistic") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_universe

ggsave(here("plots", "s_statistic.pdf"), width = 4, height = 2)

s_p_value <- get_p(s_experimental, s_observed, n_permutations + 1)
s_p_value < 0.05

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

pdf(here("plots", "cooccurrence.pdf"), width = 7, height = 7)
corrplot::corrplot(cooccurrence, 
                   is.corr = FALSE, method = "pie", 
                   type = "upper", diag = FALSE, 
                   order = "FPC", 
                   tl.col = "gray20", tl.cex = 0.8, 
                   cl.cex = 0.6)
dev.off()

# pdf(here("plots", "cooccurrence_alphabet.pdf"), width = 7, height = 7)
# corrplot::corrplot(cooccurrence, 
#                    is.corr = FALSE, method = "pie", 
#                    type = "upper", diag = FALSE, 
#                    order = "alphabet", 
#                    tl.col = "gray20", tl.cex = 0.8, 
#                    cl.cex = 0.6)
# dev.off()

# non random variables =========================================================
non_rand_vars <- v_statistic %>% filter(signif == TRUE) %>% pull(abbrv)
non_rand_cooc <- list(pos = cooccurrence[non_rand_vars, non_rand_vars],
                      neg = cooccurrence[non_rand_vars, non_rand_vars])
non_rand_cooc$pos[non_rand_cooc$pos < 0] <- 0
non_rand_cooc$neg[non_rand_cooc$neg > 0] <- 0
non_rand_cooc$neg <- abs(non_rand_cooc$neg)

g_cooc_positive <- simplify(graph_from_adjacency_matrix(non_rand_cooc$pos, 
                                                        mode = "undirected", 
                                                        weighted = TRUE))
g_cooc_negative <- simplify(graph_from_adjacency_matrix(non_rand_cooc$neg, 
                                                        mode = "undirected",
                                                        weighted = TRUE))
pdf(here("plots", "cooc_networks.pdf"), width = 12)
par(mfrow = c(1, 2), mar = c(rep(2, 4)))
plot(g_cooc_positive, 
     vertex.shape = "circle", 
     vertex.color = "white",
     # vertex.size = eigen_centrality(g_cooc_positive)$vector*20,
     vertex.label.family = "sans", 
     vertex.label.cex = .6,
     vertex.label.color = "black",
     edge.width = edge_attr(g_cooc_positive)$weight^1.2,
     mark.groups = cluster_fast_greedy(g_cooc_positive),
     mark.col = "gray90",
     mark.border = NA)
title(main = "copresence", cex.main = 1, family = "sans", font.main = 1)
plot(g_cooc_negative, 
     vertex.shape = "circle", 
     vertex.color = "white",
     # vertex.size = eigen_centrality(g_cooc_negative)$vector*20,
     vertex.label.family = "sans", 
     vertex.label.cex = .6,
     vertex.label.color = "black",
     edge.width = edge_attr(g_cooc_negative)$weight^1.2,
     mark.groups = cluster_fast_greedy(g_cooc_negative),
     mark.col = "gray90", 
     mark.border = NA)
title(main = "coabsence", cex.main = 1, family = "sans", font.main = 1)
dev.off()

# output =======================================================================
writeLines(non_rand_vars, here("data/temp", "non_random_vars.txt"))
