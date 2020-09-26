# Project "Vedrovice"
# Script nr. 3.2
# SPATIAL ARRANGEMENT
# author: Petr Pajdla
# Spatial organization within the cemetery

library(here)
library(spatstat)
library(sf)
library(tidyverse)

# number of simulations ---------------------------------------------------
nsim = 999
nrank = 5

# k...nrank
# M...nsim
# probability 2 * k / M + 1
# thus for 999 iterations 2 * 25 / 999 + 1 = 0.05 (5% probability)

# theme for ggplot graphics
theme_universe <- theme(panel.border = element_rect(colour = "black", 
                                                    fill = NA, 
                                                    size = 0.8),
                        panel.background = element_blank(),
                        line = element_blank(),
                        strip.background = element_blank(), 
                        strip.text = element_text(face = "italic"),
                        axis.text.y = element_blank(), 
                        axis.title.y = element_blank(), 
                        panel.spacing = unit(1, "lines"))

# read data and layouts --------------------------------------------------------
ved <- read_rds(here("data/temp", "vedrovice_dataset.RDS"))

ved_sf <- read_sf(here("data/temp", "layout.shp")) %>% 
  filter(pres != "dist.")
ved_exc <- read_sf(here("data/temp", "window.shp")) %>% 
  as.owin(ved_exc)

ved_sf <- bind_cols(ved_sf, 
                    x = st_coordinates(ved_sf)[, 1], 
                    y = st_coordinates(ved_sf)[, 2])

# marks ------------------------------------------------------------------------
ved_marks <- ved_sf %>% 
  select(sex, age = age_sim, origin, starts_with("cat")) %>% 
  st_drop_geometry() %>% 
  mutate(sex = factor(sex, levels = c("n. a.", "F", "M", "ind.")),
         age = factor(age, levels = c("juv.", "ad.", "mat.", "ind.")),
         origin = factor(origin, levels = c("local", "non-local", "ind.")),
         across(starts_with("cat_"), factor))

# create ppp object ------------------------------------------------------------
# ved_pp <- ppp(x = ved_layout$layout_x, y = ved_layout$layout_y,
#               window = window_poly,
#               marks = ved_marks)
ved_pp <- ppp(x = ved_sf$x, y = ved_sf$y, 
              window = ved_exc,
              marks = ved_marks)

# plot(ved_pp, use.marks = FALSE)

# point pattern analysis
# functions ====================================================================

# function L
# scaled centered Ripley's K, based on O'Sullivan & Unwin 2010, p. 147
# @x = ppp object
estimate_L <- function(x) {
  stopifnot(verifyclass(x, "ppp"))
  l <- envelope(x, Kest, nrank = nrank, nsim = nsim) %>%
    as_tibble() %>%
    mutate(r = r,
           obs = sqrt(obs/pi) - r,
           theo = sqrt(theo/pi) - r,
           lo = sqrt(lo/pi) - r,
           hi = sqrt(hi/pi) - r)
  return(l)
}

# note to self: in most of the sources, the distance r is not subtracted, 
# see Nakoinz & Knitter 2016, p. 138
# estimate_L <- function(x) {
#   stopifnot(verifyclass(x, "ppp"))
#   l <- envelope(x, Kest, nrank = 2, nsim = 99) %>%
#     fix_data_frame() %>%
#     mutate(r = r,
#            obs = sqrt(obs/pi),
#            theo = sqrt(theo/pi),
#            lo = sqrt(lo/pi),
#            hi = sqrt(hi/pi))
#   return(l)
# }

# plot function estimate
# @x = tidy tibble of estimated values
# @fun = function name
plot_estimate <- function(x, fun) {
  stopifnot(verifyclass(x, "tbl_df"))
  stopifnot(all(names(x) %in% c("r", "obs", "theo", "lo", "hi")))
  ylabel = paste0(fun, "(r)")
  p1 <- x %>% 
    ggplot() +
    geom_ribbon(aes(r, ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
    geom_line(aes(r, theo), linetype = 2) +
    geom_line(aes(r, obs), size = 0.8) + 
    annotate("text", -Inf, +Inf, label = ylabel, 
             hjust = -0.2, vjust = 1.4, 
             size = 6, fontface = "italic") +
    # labs(y = ylabel) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_universe +
    theme(axis.title.x = element_blank())
  # ggsave(here("plots", paste0("pointprocess_", fun, ".pdf")))
  p1
}

plot_estimate_facet <- function(x) {
  stopifnot(verifyclass(x, "tbl_df"))
  stopifnot(all(names(x) %in% c("fun", "r", "obs", "theo", "lo", "hi")))
  p1 <- x %>% 
    ggplot() +
    geom_ribbon(aes(r, ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
    geom_line(aes(r, theo), linetype = 2) +
    geom_line(aes(r, obs), size = 0.8) + 
    facet_wrap(vars(fun), nrow = 1, scales = "free") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "r (m)") +
    theme_universe
  p1
}

# overall structure ============================================================
# density based ----------------------------------------------------------------
# quadrat test
# quadrat.test(ved_pp)
# quadrat.test(ved_pp, alternative = "clustered")
# quadrat.test(ved_pp, alternative = "regular")
# quadrat test suggests clustering

# distance based ---------------------------------------------------------------
# estimating function G
ved_g <- envelope(ved_pp, Gest, nrank = nrank, nsim = nsim) %>%
  as_tibble()

plot_estimate(ved_g, "G")

# # F function
ved_f <- envelope(ved_pp, Fest, nrank = nrank, nsim = nsim) %>%
  as_tibble()

plot_estimate(ved_f, "F")

# estimating Ripley's K function
# ved_k <- Kest(ved_pp)
ved_k <- envelope(ved_pp, Kest, nrank = nrank, nsim = nsim) %>%
  as_tibble()

plot_estimate(ved_k, "K")

# estimating L function
ved_l <- estimate_L(ved_pp)

plot_estimate(ved_l, "L")

# estimate T fun/stat
ved_t <- envelope(ved_pp, Tstat, nrank = nrank, nsim = nsim) %>%
  as_tibble()

plot_estimate(ved_t, "T")

# export combined figure -------------------------------------------------------
# using grid and individual plots
grid_fns <- gridExtra::grid.arrange(plot_estimate(ved_g, "G"),
                                    plot_estimate(ved_f, "F"),
                                    plot_estimate(ved_k, "K"),
                                    plot_estimate(ved_t, "T"),
                                    nrow = 1)

# using facets
grid_fns <- bind_rows("G(r)" = ved_g,
                      "F(r)" = ved_f,
                      "K(r)" = ved_k,
                      "T(r)" = ved_t, .id = "fun") %>%
  mutate(fun = as_factor(fun)) %>%
  plot_estimate_facet()

grid_fns

ggsave(here("plots", "pp_funs.pdf"), width = 12, height = 3)

# for different marks ==========================================================
multiplecrossfunction <- function(pp, var, crossfun, nsim, nrank) {
  input <- subset(pp, select = var, drop = TRUE)
  input <- input[input$marks != "ind."]
  input$marks <- droplevels(input$marks)
  # cross function
  crossf <- alltypes(input, crossfun, nsim = nsim, nrank = nrank, envelope = TRUE)
  # proper names
  plot_names <- map(crossf[["fns"]], attr, "fname") %>% 
    map(2) %>% 
    map(as_tibble) %>% 
    bind_rows() %>% 
    mutate(value = str_remove(value, "^list\\("),
           value = str_remove(value, "\\)$"),
           value = str_replace(value, ",", " - "),
           value = str_remove_all(value, "\\.(?=\\.)")) %>% 
    pull(value)
  # output
  res <- bind_rows(map(crossf[["fns"]], as_tibble), .id = "id") %>% 
    mutate(name = plot_names[as.integer(id)]) %>% 
    separate(name, into = c("from", "to"), sep = " - ")
  
  return(res)
}

# sex
mppk_sex <- multiplecrossfunction(ved_pp, "sex", "Kcross", nsim, nrank)
# mppg_sex <- multiplecrossfunction(ved_pp, "sex", "Gcross", nsim, nrank)

# age
mppk_age <- multiplecrossfunction(ved_pp, "age", "Kcross", nsim, nrank)
# mppg_age <- multiplecrossfunction(ved_pp, "age", "Gcross", nsim, nrank)

# origin
mppk_orig <- multiplecrossfunction(ved_pp, "origin", "Kcross", nsim, nrank)
# mppg_orig <- multiplecrossfunction(ved_pp, "origin", "Gcross", nsim, nrank)

# K funs
g1 <- mppk_sex %>% 
  ggplot(aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
  geom_line(aes(y = theo), linetype = 2) +
  geom_line(aes(y = obs), size = 0.8) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_universe +
  facet_grid(from ~ to) +
  labs(x = "r (m)", title = "Sex")

g2 <- mppk_age %>% 
  ggplot(aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
  geom_line(aes(y = theo), linetype = 2) +
  geom_line(aes(y = obs), size = 0.8) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_universe +
  facet_grid(from ~ to) +
  labs(x = "r (m)", title = "Age category")
  
g3 <- mppk_orig %>% 
  ggplot(aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
  geom_line(aes(y = theo), linetype = 2) +
  geom_line(aes(y = obs), size = 0.8) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_universe +
  facet_grid(from ~ to) +
  labs(x = "r (m)", title = "Origin")

pdf(here::here("plots", "pp_crossk.pdf"), width = 9, height = 14)
gridExtra::grid.arrange(top = grid::textGrob("Multitype K(r) function"), 
                        g1, g2, g3, ncol = 1)
dev.off()

# # G funs
# mppg_sex %>% 
#   ggplot(aes(x = r)) +
#   geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
#   geom_line(aes(y = theo), linetype = 2) +
#   geom_line(aes(y = obs), size = 0.8) + 
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme_universe +
#   facet_grid(from ~ to) +
#   labs(x = "r (m)", title = "Multitype G(r)")
# 
# mppg_age %>% 
#   ggplot(aes(x = r)) +
#   geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
#   geom_line(aes(y = theo), linetype = 2) +
#   geom_line(aes(y = obs), size = 0.8) + 
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme_universe +
#   facet_grid(from ~ to) +
#   labs(x = "r (m)", title = "Multitype G(r)")
#   
# mppg_orig %>% 
#   ggplot(aes(x = r)) +
#   geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
#   geom_line(aes(y = theo), linetype = 2) +
#   geom_line(aes(y = obs), size = 0.8) + 
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme_universe +
#   facet_grid(from ~ to) +
#   labs(x = "r (m)", title = "Multitype G(r)")


# combined vars -----------------------------------------------------------

mppk_sa <- multiplecrossfunction(ved_pp, "cat_sa", "Kcross", nsim, nrank)
mppk_os <- multiplecrossfunction(ved_pp, "cat_os", "Kcross", nsim, nrank)
mppk_oa <- multiplecrossfunction(ved_pp, "cat_oa", "Kcross", nsim, nrank)
mppk_osa <- multiplecrossfunction(ved_pp, "cat_osa", "Kcross", nsim, nrank)

mppk_osa %>%
  ggplot(aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
  geom_line(aes(y = theo), linetype = 2) +
  geom_line(aes(y = obs), size = 0.8) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_universe +
  facet_grid(from ~ to) +
  labs(x = "r (m)", title = "Multitype K(r)")
