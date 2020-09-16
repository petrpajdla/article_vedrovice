# Project "Vedrovice"
# Script nr. 1.2
# LOCAL AND NON-LOCAL Sr RANGES
# author: Petr Pajdla
# Simplify and clean age groups 
# Assign bodies into local vs non-local based on Sr range

library(tidyverse)
library(here)


# data --------------------------------------------------------------------

ved <- read_rds(here("data/temp/vedrovice_dataset.RDS"))


# edit age ----------------------------------------------------------------

age <- ved$orig_dataset %>% 
  select(id_burial, age_interval) %>% 
  mutate(age = if_else(age_interval == "Adult", "20-40", age_interval),
         age = if_else(age_interval == "±1", "0-2", age),
         age = if_else(age_interval == "20+", "20-30", age),
         age = if_else(age_interval == "50+", "50-60", age),
         age = if_else(age_interval == "60+", "60-70", age),
         age = if_else(age_interval == "9", "8-10", age),
         age = if_else(age_interval == "30", "25-35", age),
         age = if_else(age_interval == "20", "15-25", age),
         age_min = as.numeric(str_extract(age, "^[:digit:]*(?=-|–)")),
         age_max = as.numeric(str_extract(age, "(?<=-|–)[:digit:]*$")),
         id_burial = as.character(id_burial)) %>% 
  select(-age_interval, -age)


# group by age ------------------------------------------------------------

ved_sr <- ved$metadata %>% 
  select(id_burial, sex, age_cat) %>% 
  left_join(tibble(rownames_to_column(ved$cont_vars)), 
            by = c("id_burial" = "rowname")) %>% 
  left_join(age) %>% 
  mutate(age_mean = (age_min + age_max) / 2,
         group = if_else(age_mean < 18, "juvenile", as.character(sex))
         ) %>% 
  filter(!is.na(group), !group %in% c("n. a.", "ind."))

local <- ved_sr %>% 
  filter(group == "juvenile",
         !is.na(sr)) %>%
  summarise(n(), 
            mean = mean(sr, na.rm = TRUE), 
            sd = sd(sr, na.rm = TRUE)) %>% 
  mutate(min = mean - 2 * sd,
         max = mean + 2 * sd)


# local range of sr -------------------------------------------------------

ved_local <- ved_sr %>% 
  mutate(origin = if_else(sr > local$min & sr < local$max, "local", "non-local"))
  
ved_local %>% 
  ggplot(aes(x = sr)) +
  annotate("rect", xmin = local$min, xmax = local$max, 
           ymin = -Inf, ymax = +Inf,
           alpha = .2) +
  # geom_point(aes(y = age_min, color = local)) +
  # geom_point(aes(y = age_max, color = local)) +
  geom_linerange(aes(ymin = age_min, ymax = age_max)) +
  facet_wrap(~group, nrow = 3) +
  ggthemes::theme_clean()

ggsave(here("plots", "local_range.pdf"))

# export ------------------------------------------------------------------

# local range defined as 2 SD from the mean of juveniles
ved_local %>% 
  select(id_burial, sr, age_min, age_max, age_mean, group, origin) %>%
  write_rds(here("data/temp/local.rds"))


