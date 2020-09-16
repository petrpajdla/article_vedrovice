# Project "Vedrovice"
# Script nr. 
# EXPLORE EXCEPTIONALITY
# author: Petr Pajdla
# Explore exceptionality of burials, what does it cause etc...

# libraries
library(tidyverse)
library(here)

# data 
base <- read_rds(here("data/temp", "vedrovice_dataset.RDS"))

ei <- read_csv(here("data/temp", "exceptionality.csv"), col_types = "cdff")

ei_fct_order <- ei %>% group_by(fct) %>% 
  summarise(mean_ei = mean(ei))

v_vals <- read_csv(here("data/temp", "v_statistics.csv"))


# grave goods in various ei  ----------------------------------------------
ved_gg <- left_join(ei, 
                    as_tibble(base$bin_vars$count_mat[, base$bin_vars$over5], 
                              rownames = "burial"),
                    by = c("burial")) %>% 
  select(-ei, -ei_cluster)

# proportion of grave goods in a given group
gg_prop <- function(df, id, group) {
  temp <- df %>% pivot_longer(-c(!!sym(id), !!sym(group)), 
                              names_to = "gg", 
                              values_to = "presence")
  temp_sums <- temp %>% group_by(gg) %>% 
    summarise(total = sum(presence)) %>% 
    ungroup()
  temp %>% group_by(!!sym(group), gg) %>% 
    summarise(sum = sum(presence)) %>% 
    group_map(~ mutate(.x, prop = sum/temp_sums$total), .keep = TRUE) %>% 
    bind_rows()
}

# reordering factors
ved_gg_prop <- gg_prop(ved_gg, "burial", "fct") %>% 
  left_join(ei_fct_order) %>% 
  mutate(fct = fct_reorder(fct, desc(mean_ei)),
         gg_order = as.integer(fct_relevel(gg, v_vals$variable)),
         gg_long = unname(base$var_names$long[gg]),
         gg_long = fct_reorder(gg_long, desc(gg_order))) %>% 
  select(-mean_ei)


# plotting
ved_gg_prop %>% ggplot(aes(y = prop, x = gg_long)) +
  geom_col(fill = "white", color = "black") +
  facet_wrap(vars(fct), nrow = 1) +
  coord_flip() +
  labs(x = "artefact category", y = "proportion") + 
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     minor_breaks = c(0.25, 0.75), 
                     labels = c(0, 0.5, 1), 
                     limits = c(0, 1))

ggsave(here("plots", "exc_gg_prop.pdf"), width = 8, height = 5)

ved_gg_prop %>% group_by(fct, gg) %>% 
  summarise(prop = sum(prop)) %>% 
  write_csv(here("data/temp", "exc_gg_prop.csv"))

