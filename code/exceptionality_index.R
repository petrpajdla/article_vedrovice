# Project "Vedrovice"
# Script nr. 3
# COUNT EXCEPTIONAITY INDEX OF BURIALS
# author: Petr Pajdla
# Exceptionality index of individual burials in order to create groups...

# packages =====================================================================
library(dplyr)

# functions ====================================================================
# root mean square (quadratic mean) for vector x
rms <- function(x) {
  x <- sqrt(mean(x^2))
  return(x)
}

# data =========================================================================
ved <- readRDS("./data/vedrovice_dataset.RDS")
non_rand_vars <- readLines("./data/output/non_random_vars.txt")
# exceptionality index =========================================================
# non random variables
art_pca <- prcomp(ved$bin_vars$bin_mat[, non_rand_vars], 
                  scale. = FALSE, center = FALSE)
summary(art_pca)
ggbiplot::ggbiplot(art_pca)

# # for_pca <- ved %>% select(c("pit_len", "pit_wid", "pit_dep", vars_binary)) %>% 
# #   mutate_at(c("pit_len", "pit_wid", "pit_dep"), scale) %>% 
# #   mutate_all(replace_na)
# 
# # pca object
# pca <- princomp(for_pca, scores = T)
# # summary(pca)

# counting ei for individual burials
ei <- apply(art_pca$x[, 1:9], 1, rms) / max(apply(art_pca$x[, 1:9], 1, rms))
hist(as.matrix(ei))

ei_tbl <- tibble(burial = names(ei), ei)
ei_tbl %>% arrange(desc(ei))

ggplot(ei_tbl, aes(x = forcats::fct_reorder(burial, ei), y = ei)) +
  geom_point() +
  coord_flip() +
  xlab("burial nr.")

ei_tbl %>% arrange(desc(ei)) %>% View()

output <- ei_tbl %>% mutate(group = if_else(ei > 0.6, "group1",
                                  if_else((ei < 0.6 & ei > 0.5), "group2", 
                                          if_else((ei < 0.5 & ei > 0.4), "group3",
                                                  if_else((ei < 0.4 & ei > 0.2), "group4", "group5")))))

readr::write_csv(output, "./data/output/exceptionality.csv")
