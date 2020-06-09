# Project "Vedrovice"
# Script nr. 5
# EXPLORE EXCEPTIONALITY
# author: Petr Pajdla
# Explore exceptionality of burials, what does it cause etc...

# libraries
library(tidyverse)

# data 
base <- read_rds("./data/vedrovice_dataset.RDS")

ei <- read_csv("./data/output/exceptionality.csv")

non_random_vars <- read_lines("./data/output/non_random_vars.txt")

# explore...
# spatial distribution
spatial_clustering <- hclust(dist(base$layout, method = "euclidean"))
spatial_groups <- cutree(spatial_clustering, k = 8)


ggplot(as_tibble(base$layout), aes(layout_x, layout_y, color = factor(spatial_groups))) +
  geom_point() + 
  ggforce::geom_mark_hull(expand = unit(2.2, "mm"))

ggplot(as_tibble(base$layout), aes(layout_x, layout_y, color = factor(spatial_groups))) +
  geom_point() + 
  ggforce::geom_mark_ellipse(expand = unit(2, "mm"))

ved_spatial <- bind_cols(as_tibble(base$layout), ei = ei$group, spat_gr = spatial_groups)

ggplot(ved_spatial, aes(layout_x, layout_y, color = ei)) +
  geom_point() + 
  facet_wrap(~factor(spat_gr))

ggplot() +
  geom_point(ved_spatial[, c("layout_x", "layout_y")], mapping = aes(layout_x, layout_y), alpha = 0.1) +
  geom_point(ved_spatial, mapping = aes(layout_x, layout_y, color = ei)) +
  # stat_density2d() +
  facet_wrap(~ei)

# spat distr. using kmeans (experiment)
spat_kmeans <- kmeans(base$layout, centers = 8, nstart = 200)

scree1 <- function(df) {
  x <- vector(mode = "double", length = 12)
   for (i in 1:12) {
     x[i] <- kmeans(df, centers = i, nstart = 20)$tot.withinss
   }
  return(x)
}

plot(scree1(base$layout), type = "b")

ggplot(ved_spatial, aes(layout_x, layout_y, color = factor(spat_kmeans$cluster))) +
  geom_point() +
  ggforce::geom_mark_ellipse(expand = unit(2, "mm"))


# groups items...
ved_items <- bind_cols(id = ei$burial, 
                       as_tibble(base$bin_vars$bin_mat), 
                       ei = ei$group)

# ad hoc function
prop_item <- function(df) {
  df_gathered <- df %>% gather(-id, -ei, key = "item", value = "presence")
  temp_total <- df_gathered %>% group_by(item) %>% 
    summarise(tot = sum(presence))
  temp_groups <- df_gathered %>% group_by(ei, item) %>% 
    summarise(sum = sum(presence)) %>% 
    group_map(~ mutate(.x, prop = sum/temp_total$tot))
  return(temp_groups)
}

item_prop <- prop_item(ved_items) %>% 
  bind_rows(.id = "ei") %>% 
  mutate(ei = str_c("group", ei))

# randomly distributed items
ggplot(item_prop, aes(y = prop, x = factor(item))) +
  geom_col(fill = "white", color = "black") +
  facet_wrap(~ ei, nrow = 1) +
  coord_flip() +
  labs(x = "artefact category", y = "proportion")

item_prop %>% group_by(item) %>% summarise(sumx = sum(prop))

# non-randomly distributed items
item_prop %>% filter(item %in% non_random_vars) %>% 
  ggplot(aes(y = prop, x = factor(item))) +
  geom_col(fill = "white", color = "black") +
  facet_wrap(~ ei, nrow = 1) +
  coord_flip() +
  labs(x = "artefact category", y = "proportion")


          