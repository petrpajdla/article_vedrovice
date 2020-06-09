# Project "Vedrovice"
# Script nr. 6
# CREATE GROUPS BASED ON COMMON ARTEFACTS
# author: Petr Pajdla
# Create groups based on common artefacts

# pckg
library(tidyverse)

# data
ved <- read_rds("./data/vedrovice_dataset.RDS")

non_rand_vars <- read_lines("./data/output/non_random_vars.txt")

ei <- read_csv("./data/output/exceptionality.csv")

# experiment

ved_bin_dist_nr <- dist(ved$bin_vars$bin_mat[, non_rand_vars], method = "binary")
ved_bin_dist <- dist(ved$bin_vars$bin_mat, method = "binary")

ved_hclust_nr <- hclust(ved_bin_dist_nr, method = "ward.D2")
ved_hclust <- hclust(ved_bin_dist, method = "ward.D2")

par(mfcol = c(1, 2))
plot(ved_hclust)
plot(ved_hclust_nr)

groups_common_items <- cutree(ved_hclust_nr, k = 6)

ved_gathered <- bind_cols(id = ved$metadata$id_burial, 
                          as_tibble(ved$bin_vars$bin_mat), 
                          gr_item = groups_common_items,
                          gr_ei = ei$group) %>% 
  gather(-gr_item, -gr_ei, -id, key = "item", value = "presence")

# ad hoc function
prop_item <- function(df) {
  temp_total <- df %>% group_by(item) %>% 
    summarise(tot = sum(presence))
  temp_groups <- df %>% group_by(gr_item, item) %>% 
    summarise(sum = sum(presence)) %>% 
    group_map(~ mutate(.x, prop = sum/temp_total$tot))
  return(temp_groups)
}

item_prop <- prop_item(ved_gathered) %>% 
  bind_rows(.id = "gr_item") %>% 
  mutate(gr_item = str_c("group", gr_item))

# randomly distributed items
ggplot(item_prop, aes(y = prop, x = item)) +
  geom_col(fill = "white", color = "black") +
  facet_wrap(~ gr_item, nrow = 1) +
  coord_flip() +
  labs(x = "artefact category", y = "proportion")

item_prop %>% group_by(item) %>% summarise(sumx = sum(prop))

# non-randomly distributed items
item_prop %>% filter(item %in% non_random_vars) %>% 
  ggplot(aes(y = prop, x = item)) +
  geom_col(fill = "white", color = "black") +
  facet_wrap(~ gr_item, nrow = 1) +
  coord_flip() +
  labs(x = "artefact category", y = "proportion")








