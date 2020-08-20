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

ved_sf <- read_sf(here("data/temp", "layout.shp"))
ved_exc <- read_sf(here("data/temp", "window.shp")) %>% 
  as.owin(ved_exc)

ved_sf <- bind_cols(ved_sf, 
                    x = st_coordinates(ved_sf)[, 1], 
                    y = st_coordinates(ved_sf)[, 2])

# marks ------------------------------------------------------------------------
ved_marks <- ved$metadata %>% 
  select(sex) %>% 
  mutate(sex = if_else(is.na(sex), "ind.", as.character(sex)),
         sex = as_factor(sex))

# create ppp object ------------------------------------------------------------
# ved_pp <- ppp(x = ved_layout$layout_x, y = ved_layout$layout_y,
#               window = window_poly,
#               marks = ved_marks)
ved_pp <- ppp(x = ved_sf$x, y = ved_sf$y, 
              window = ved_exc,
              marks = ved_marks)

plot(ved_pp, use.marks = FALSE)

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
quadrat.test(ved_pp)
quadrat.test(ved_pp, alternative = "clustered")
quadrat.test(ved_pp, alternative = "regular")
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

ggsave(here("plots", "pointprocess_fun.pdf"), width = 12, height = 3)

# for different marks ==========================================================
# ved_pp_split <- split(ved_pp)
# plot(ved_pp_split)
# # density based ----------------------------------------------------------------
# # estimate_L(ved_pp_split[["M"]])
# # 
# # estimate_L(ved_pp_split[["F"]])
# # 
# # estimate_L(ved_pp_split[["n. a."]])
# 
# ved_m <- ved_pp_split[["M"]]
# ved_f <- ved_pp_split[["F"]]
# ved_na <- ved_pp_split[["n. a."]]
# 
# estimate_all <- function(pp, nrank, nsim) {
#   g <- envelope(pp, Gest, nrank = nrank, nsim = nsim) %>%
#     as_tibble()
#   f <- envelope(pp, Fest, nrank = nrank, nsim = nsim) %>%
#     as_tibble()
#   k <- envelope(pp, Kest, nrank = nrank, nsim = nsim) %>%
#     as_tibble()
#   t <- envelope(pp, Tstat, nrank = nrank, nsim = nsim) %>%
#     as_tibble()
#   
#   grid_fns <- bind_rows("G(r)" = g,
#                         "F(r)" = f,
#                         "K(r)" = k,
#                         "T(r)" = t, .id = "fun") %>%
#     mutate(fun = as_factor(fun))
#   
#   return(grid_fns)
# }
# 
# est_m <- estimate_all(ved_m, nrank = 1, nsim = 99)
# est_f <- estimate_all(ved_f, nrank = 1, nsim = 99)
# est_na <- estimate_all(ved_na, nrank = 1, nsim = 99)
# 
# p_m <- est_m %>% plot_estimate_facet() + labs(title = "M")
# p_f <- est_f %>% plot_estimate_facet() + labs(title = "F")
# p_na <- est_na %>% plot_estimate_facet() + labs(title = "n. a.")
# 
# p_all <- gridExtra::grid.arrange(p_m, p_f, p_na)
# 
# ggsave(here("plots/pointprocess_fun_sex.pdf"), plot = p_all, 
#        width = 12, height = 9)

# Cross functions  --------------------------------------------------------
# plot(envelope(ved_pp, Kcross, funargs = list(i = "M", j = "F")))
# plot(envelope(ved_pp, Kcross, funargs = list(i = "F", j = "M")))
# plot(envelope(ved_pp, Kcross, funargs = list(i = "M", j = "n. a.")))
# plot(envelope(ved_pp, Kcross, funargs = list(i = "n. a.", j = "M")))
# plot(envelope(ved_pp, Kcross, funargs = list(i = "F", j = "n. a.")))
# plot(envelope(ved_pp, Kcross, funargs = list(i = "n. a.", j = "F")))

all_k <- alltypes(ved_pp, Kcross, nsim = nsim, nrank = nrank, envelope = TRUE)
all_g <- alltypes(ved_pp, Gcross, nsim = nsim, nrank = nrank, envelope = TRUE)

# plot(all_k)
# plot(all_g)

all_k_tibble <- bind_rows(lapply(all_k[["fns"]][c(6:8, 10:12, 14:16)], 
                                 as_tibble), .id = "id")
all_g_tibble <- bind_rows(lapply(all_g[["fns"]][c(6:8, 10:12, 14:16)], 
                                 as_tibble), .id = "id")

panelnames <- c("F-F", "F-M", "F-n.a.",
                "M-F", "M-M", "M-n.a.",
                "n.a.-F", "n.a.-M", "n.a.-n.a.")

multiK <- all_k_tibble %>% 
  mutate(id = panelnames[as.integer(id)]) %>% 
  separate(id, into = c("from", "to"), sep = "-") %>% 
  ggplot(aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
  geom_line(aes(y = theo), linetype = 2) +
  geom_line(aes(y = obs), size = 0.8) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "r (m)", title = "Multitype K(r)") +
  theme_universe +
  facet_grid(from ~ to)

multiG <- all_g_tibble %>% 
  mutate(id = panelnames[as.integer(id)]) %>% 
  separate(id, into = c("from", "to"), sep = "-") %>% 
  ggplot(aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.6) +
  geom_line(aes(y = theo), linetype = 2) +
  geom_line(aes(y = obs), size = 0.8) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "r (m)", title = "Multitype G(r)") +
  theme_universe +
  facet_grid(from ~ to)
  
ggsave(filename = here("plots/pointprocess_multiK_sex.pdf"), plot = multiK, 
       width = 9, height = 9)
ggsave(filename = here("plots/pointprocess_multiG_sex.pdf"), plot = multiG, 
       width = 9, height = 9)


# # fuuu --------------------------------------------------------------------
# 
# # density
# plot(density(ved_pp_split), equal.ribbon = TRUE)
# plot(density(ved_pp_split), equal.ribbon = TRUE)
# 
# str(density(ved_pp_split))
# 
# plot(density(split(ved_pp)[2:4]))
