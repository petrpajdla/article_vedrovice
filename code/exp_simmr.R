library(tidyverse)
library(simmr)

# vignette("quick_start", package = "simmr")
# vignette("simmr", package = "simmr")

# load data - vedrovice
x <- read_rds("./data/vedrovice_dataset.RDS")

# prepare data (matrix, grouping factor)
mix_in <- cbind(x$metadata$sex, x$cont_vars$cont_vars[, c("d13c", "d15n")]) %>% 
  drop_na()

mix <- as.matrix(mix_in[, 2:3])
grp <- mix_in$`x$metadata$sex`

# create sources data (Bickle, Whittle 2013, p. 121)
source_data <- tribble(
  ~spec, ~Meand13c, ~SDd13c, ~Meand15n, ~SDd15n, ~n,
  "cattle",-20.2,  0.29, 6.2,    1.04, 3,
  "sheep", -19.8,  0.18, 5.9,    0.38, 3,
  "pig",   -20.4,  0.53, 8.2,    0.86, 2
)

s_names <- source_data$spec
s_means <- source_data[, c("Meand13c", "Meand15n")]
s_sd <- source_data[, c("SDd13c", "SDd15n")]

# create simmr object
simmr_obj <- simmr_load(mix, s_names, s_means, s_sd, group = grp)

# overall plot
plot(simmr_obj, group = 1:3)

# fit model
simmr_out <- simmr_mcmc(simmr_obj)

# validate - values should be close to 1
summary(simmr_out, type='diagnostics', group = 1:3)

# blue y values should fall within light blue y_rep intervals
posterior_predictive(simmr_out, group = 1)
posterior_predictive(simmr_out, group = 2)
posterior_predictive(simmr_out, group = 3)

# prior prob. distr. viz
prior_viz(simmr_out, group = 1)
prior_viz(simmr_out, group = 2)
prior_viz(simmr_out, group = 3)

summary(simmr_out, type='statistics', group = 1:3)
summary(simmr_out, type='quantiles', group = 1:3)

# plots
plot(simmr_out, type = 'density', group = 1)
plot(simmr_out, type = 'density', group = 2)
plot(simmr_out, type = 'density', group = 3)

# matrices
# large  negative correlations - model is not able to differentiate the two groups
plot(simmr_out, type = 'matrix', group = 1)
plot(simmr_out, type = 'matrix', group = 2)
plot(simmr_out, type = 'matrix', group = 3)

# boxplots
plot(simmr_out, type = 'boxplot', group = 1, title='simmr output group 1')
plot(simmr_out, type = 'boxplot', group = 2, title='simmr output group 2')
plot(simmr_out, type = 'boxplot', group = 3, title='simmr output group 3')

# comparison of sources
compare_sources(simmr_out, source_names=c("cattle", "sheep", "pig"), group = 1)
compare_sources(simmr_out, source_names=c("cattle", "sheep", "pig"), group = 2)
compare_sources(simmr_out, source_names=c("cattle", "sheep", "pig"), group = 3)

# comparison of groups
compare_groups(simmr_out, source = 'cattle', groups=2:3)
compare_groups(simmr_out, source = 'pig', groups=2:3)
compare_groups(simmr_out, source = 'sheep', groups=2:3)

