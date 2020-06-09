library(tidyverse)
library(MixSIAR)

x <- read_rds("./data/vedrovice_dataset.RDS")

write_csv(bind_cols(x$metadata[, c("id_burial", "sex")],
                    x$cont_vars$cont_vars[, c("d13c", "d15n")]), 
            "./data/output/mixing.csv")
          

source_data <- tribble(
  ~spec, ~Meand13c, ~SDd13c, ~Meand15n, ~SDd15n, ~n,
  "cattle",-20.2,  0.29, 6.2,    1.04, 3,
  "sheep", -19.8,  0.18, 5.9,    0.38, 3,
  "pig",   -20.4,  0.53, 8.2,    0.86, 2,
  "fishy", -18, 0.3, 13, 2.5, 2
)

write_csv(source_data, "./data/output/source_data.csv")

discr <- tribble(
  ~spec, ~Meand13c, ~SDd13c, ~Meand15n, ~SDd15n,
  "cattle", 0, 0, 0, 0,
  "sheep", 0, 0, 0, 0,
  "pig", 0, 0, 0, 0,
  "fishy", 0, 0, 0, 0
)

write_csv(discr, "./data/output/discr_data.csv")

mix <- load_mix_data(filename = "./data/output/mixing.csv", 
              iso_names = c("d13c", "d15n"), factors = c("sex"),
              fac_random = FALSE, fac_nested = FALSE, cont_effects = NULL)

srcs <- load_source_data("./data/output/source_data.csv", source_factors = NULL,
                        conc_dep = FALSE, data_type = "means", mix = mix)

discr <- load_discr_data("./data/output/discr_data.csv", mix = mix)

plot_data(mix = mix, source = srcs, discr = discr)
calc_area(srcs, mix, discr)
plot_prior(alpha.prior = 1, srcs)
write_JAGS_model("./data/output/mixsiar_model.txt", resid_err = TRUE, process_err = FALSE, mix = mix, source = srcs)

jags.1 <- run_model(run="test", mix, srcs, discr, "./data/output/mixsiar_model.txt", 
                    alpha.prior = 1, TRUE, FALSE)

