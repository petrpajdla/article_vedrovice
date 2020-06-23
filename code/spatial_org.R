# Project "Vedrovice"
# Script nr. ?
# SPATIAL ARRANGEMENT
# author: Petr Pajdla
# Spatial organization within the cemetery

library(here)
library(tidyverse)
library(broom)
library(spatstat)
library(sf)

# theme for ggplot graphics
theme_universe <- theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
                        panel.background = element_blank(),
                        line = element_blank(),
                        strip.background = element_blank(), 
                        strip.text = element_text(face = "italic"),
                        axis.text.y = element_blank(), 
                        axis.title.y = element_blank(), 
                        panel.spacing = unit(1, "lines"))

# read data
ved <- read_rds(here("data/temp", "vedrovice_dataset.RDS"))

# data prep ====================================================================
# solving problems with negative coordinaes of polygon window
ved_layout <- as_tibble(ved$layout + 10)

# create window ----------------------------------------------------------------
# not a rectangle, but a polygon closely bounding the burials using convex hull
ved_sf <- st_as_sf(bind_cols(ved_layout, ved$metadata), 
                   coords = c("layout_x", "layout_y")) %>% 
  mutate(sex = fct_relevel(sex, c("n. a.", "F", "M", "ind.")))
# ved_sf_buffer <- st_buffer(ved_sf, dist = 1.6)
# ved_sf_buffer <- st_simplify(ved_sf_buffer, dTolerance = 1)
ved_sf_hull <- st_convex_hull(st_union(st_geometry(ved_sf)))

# excavation polygon
ved_exc <- matrix(c(c(0.5, 0.5, 0, 0, -1, -1, 1.6, 1.6, 7.6, 7.6, 9.1, 9.1, 
                      12.3, 12.3, 13.5, 13.5, 12.3, 12.3, 8.6, 8.6),
                    c(-1, 2.5, 2.5, 10.3, 10.3, 14.3, 14.3, 13.6, 13.6, 
                      12.8, 12.8, 13.4, 13.4, 12.4, 12.4, 4.5, 4.5, 1.5, 1.5, -1)), 
                  ncol = 2, byrow = FALSE)
ved_exc <- st_cast(st_as_sf(st_geometry(st_multipoint(ved_exc + 10))), "POLYGON")

# window_rect <- owin(xrange = c(10, 22), yrange = c(10, 22.7))
window_poly <- as.owin(st_cast(ved_sf_hull, "POLYGON"))
window_exc <- as.owin(ved_exc)

ggplot(data = ved_sf) +
  geom_sf(data = ved_exc, fill = "gray90", color = NA) +
  # geom_sf(data = ved_sf_hull, fill = NA, color = "gray40", linetype = 3) +
  geom_sf(aes(shape = sex), fill = "white") + 
  scale_shape_manual(values = c(22, 21, 24, 4)) +
  # ggsflabel::geom_sf_text_repel(aes(label = id_burial)) +
  theme_void()

ggsave(here("plots", "plan_vedrovice.pdf"), scale = 2)

# marks ------------------------------------------------------------------------
ved_marks <- ved$metadata %>% 
  select(sex) %>% 
  mutate(sex = if_else(is.na(sex), "ind.", as.character(sex)),
         sex = as_factor(sex))

# create ppp object ------------------------------------------------------------
# ved_pp <- ppp(x = ved_layout$layout_x, y = ved_layout$layout_y,
#               window = window_poly,
#               marks = ved_marks)
ved_pp <- ppp(x = ved_layout$layout_x, y = ved_layout$layout_y, 
              window = window_exc,
              marks = ved_marks)

# plot(ved_pp, use.marks = FALSE)

# point pattern analysis
# functions ====================================================================

# function L
# scaled centered Ripley's K, based on O'Sullivan & Unwin 2010, p. 147
# @x = ppp object
estimate_L <- function(x) {
  stopifnot(verifyclass(x, "ppp"))
  l <- envelope(x, Kest, nrank = 2, nsim = 99) %>%
    fix_data_frame() %>%
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
    theme_universe
  p1
}

# overall structure ============================================================
# density based ----------------------------------------------------------------
# quadrat test
quadrat.test(ved_pp)
quadrat.test(ved_pp, alternative = "clustered")
quadrat.test(ved_pp, alternative = "regular")
# quadrat test suggests clustering

# estimating function G
ved_g <- envelope(ved_pp, Gest, nrank = 2, nsim = 99) %>% 
  fix_data_frame()

plot_estimate(ved_g, "G")

# # F function
ved_f <- envelope(ved_pp, Fest, nrank = 2, nsim = 99) %>%
  fix_data_frame()

plot_estimate(ved_f, "F")

# distance based ---------------------------------------------------------------
# estimating Ripley's K function
# ved_k <- Kest(ved_pp)
ved_k <- envelope(ved_pp, Kest, nrank = 2, nsim = 99) %>% 
  fix_data_frame()

plot_estimate(ved_k, "K")

# estimating L function
ved_l <- estimate_L(ved_pp)

plot_estimate(ved_l, "L")

# estimate T fun/stat
ved_t <- envelope(ved_pp, Tstat, nrank = 2, nsim = 99) %>% 
  fix_data_frame()

plot_estimate(ved_t, "T")

# export combined figure

# using grid and individual plots
# grid_fns <- gridExtra::grid.arrange(plot_estimate(ved_g, "G"),
#                                     plot_estimate(ved_f, "F"),
#                                     plot_estimate(ved_k, "K"),
#                                     plot_estimate(ved_t, "T"),
#                                     nrow = 1)

# using facets
grid_fns <- bind_rows("G(r)" = ved_g, 
                      "F(r)" = ved_f, 
                      "K(r)" = ved_k, 
                      "T(r)" = ved_t, .id = "fun") %>% 
  mutate(fun = as_factor(fun)) %>% 
  plot_estimate_facet()

ggsave(plot = grid_fns, here("plots", "pointprocess_fun.pdf"), device = "pdf",
       width = 12, height = 3)

# for different marks ==========================================================
ved_pp_split <- split(ved_pp)
plot(ved_pp_split)
# density based ----------------------------------------------------------------
estimate_L(ved_pp_split[["M"]]) %>% plot_estimate("fuu")

estimate_L(ved_pp_split[["F"]]) %>%  plot_estimate("fuu")

estimate_L(ved_pp_split[["n. a."]]) %>%  plot_estimate("fuu")

# density
plot(ved_pp_split)
plot(density(ved_pp_split), equal.ribbon = TRUE)

# pp with marks ----------------------------------------------------------------
plot(envelope(split(ved_pp)[[2]], Kest))
plot(envelope(split(ved_pp)[[3]], Kest))
plot(envelope(split(ved_pp)[[4]], Kest))

plot(density(split(ved_pp)[2:4]))

k <- Kcross(ved_pp, "M", "F")
plot(k)
plot(envelope(ved_pp, Kcross, funargs = list(i = "M", j = "F")))
