# Project "Vedrovice"
# Script nr. 2.2
# EXCEPTIONAITY INDEX FOR BURIALS
# author: Petr Pajdla
# Exceptionality index of individual burials in order to create groups...

# packages =====================================================================
library(tidyverse)
library(ggdendro)

# functions ====================================================================
# root mean square (quadratic mean) for vector x
rms <- function(x) {
  x <- sqrt(mean(x^2))
  return(x)
}

# # scaling
# normalize01 <- function(x) {
#   rg <- range(x, na.rm = TRUE)
#   (x - rg[1]) / (rg[2] - rg[1])
# }

# data =========================================================================
ved <- read_rds(here::here("data/temp/vedrovice_dataset.RDS"))
# non_rand_vars <- readLines(here("data/temp", "non_random_vars.txt"))
v_vals <- read_csv(here::here("data/temp/v_statistics.csv")) %>% 
  select(variable, p) %>% 
  mutate(weight = 1 - p)

var_weights <- v_vals %>% pull(weight)
names(var_weights) <- v_vals %>% pull(variable)


# weighting variables -----------------------------------------------------
# variables are weighted by the p value (measure of departure from norm) of v stats
# m_cont <- matrix(nrow = nrow(ved$bin_vars$count_mat), 
#                  ncol = ncol(ved$bin_vars$count_mat))

m_bin <- matrix(nrow = nrow(ved$bin_vars$bin_mat), 
                ncol = ncol(ved$bin_vars$bin_mat))

# for (i in 1:ncol(m_cont)) {
#   m_cont[, i] <- ved$bin_vars$count_mat[, i] * var_weights[i]
# }

for (i in 1:ncol(m_bin)) {
  m_bin[, i] <- ved$bin_vars$bin_mat[, i] * var_weights[i]
}

# colnames(m_cont) <- colnames(ved$bin_vars$count_mat)
colnames(m_bin) <- colnames(ved$bin_vars$bin_mat)


# exceptionality index =========================================================
# art_pca_c <- prcomp(m_cont, scale. = TRUE, center = TRUE)
art_pca_b <- prcomp(m_bin)

# summary(art_pca_c)
# ggbiplot::ggbiplot(art_pca_c)

summary(art_pca_b)
# ggbiplot::ggbiplot(art_pca_b)

pc_threshold <- which(summary(art_pca_b)$importance[3, ] >= 0.9)[1]

# counting ei for individual burials
ei_b <- apply(art_pca_b$x[, 1:pc_threshold], 1, rms) / 
  max(apply(art_pca_b$x[, 1:pc_threshold], 1, rms))
hist(as.matrix(ei_b))


# explore ei --------------------------------------------------------------

ei_tbl <- tibble(burial = ved$id_burials, ei = ei_b)
ei_tbl %>% arrange(desc(ei))

# plots
ggplot(ei_tbl, aes(x = forcats::fct_reorder(factor(burial), ei), y = ei)) +
  geom_point() +
  coord_flip() +
  xlab("burial nr.")

# relationship with amount of grave goods
ei_gg <- bind_cols(ei_tbl,
                   n_gg = rowSums(ved$bin_vars$count_mat),
                   n_gg_types = rowSums(ved$bin_vars$bin_mat))

ei_gg %>% ggplot(aes(ei, n_gg)) +
  geom_point() +
  stat_smooth()

ei_gg %>% ggplot(aes(ei, n_gg_types)) +
  geom_point() +
  stat_smooth()

cor(ei_gg$ei, ei_gg$n_gg)
cor(ei_gg$ei, ei_gg$n_gg_types)
# there is a significant correlation between exceptionality and:
# amount of grave goods - 0.78
# amount of grave good types - 0.95

# groups in ei ------------------------------------------------------------
ei_mat <- as.matrix(ei_tbl$ei)
rownames(ei_mat) <- ei_tbl$burial

ei_hcl <- hclust(dist(ei_mat), method = "ward.D2")

h <- 0.18
dendro <- ggdendrogram(ei_hcl, rotate = FALSE) +
  scale_y_sqrt(expand = expansion(c(0, 0))) +
  geom_hline(yintercept = h, color = "gray", size = 2, alpha = 0.4)

dendro

ggsave(here::here("plots", "ei_dendrogram.pdf"), width = 9, height = 3)

# cut dendrogram
ei_clusters <- cutree(ei_hcl, h = h)

ggplot(ei_tbl, aes(x = forcats::fct_reorder(factor(burial), ei), y = ei)) +
  geom_point(aes(color = factor(ei_clusters))) +
  coord_flip() +
  xlab("burial nr.") +
  theme_bw()

# output data - ei groups
output <- ei_tbl %>% bind_cols(ei_cluster = ei_clusters) %>% 
  dplyr::mutate(fct = fct_reorder(factor(ei_cluster), ei, .desc = TRUE),
                fct = fct_relabel(fct, ~ paste0(letters[1:7])))

ggplot(ei_tbl, aes(x = forcats::fct_reorder(factor(burial), ei), y = ei)) +
  geom_point(aes(color = output$fct)) +
  coord_flip() +
  xlab("burial nr.") +
  theme_bw()

ei_gg %>% ggplot(aes(ei, n_gg)) +
  geom_point(aes(color = output$fct))

ei_gg %>% ggplot(aes(ei, n_gg_types)) +
  geom_point(aes(color = output$fct))

write_csv(output, here::here("data/temp", "exceptionality.csv"))
