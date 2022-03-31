# Project "Vedrovice"
# Script nr. 1.1
# DATA INPUT AND MANIPULATION
# author: Petr Pajdla
# Script manipulates input data and saves a list with manipulated dataset (processed data)

# packages =====================================================================
library(dplyr)
library(here)
library(stringr)

# functions ====================================================================
# continuous vector to dichotomous (binarization)
binarize <- function(x, threshold = 0) {
  x[x > threshold] <- 1
  x[x <= threshold] <- 0
  return(x)
}

# data input and manipulation ==================================================
input <- readr::read_csv(here("data", "data_vedrovice_v02.csv"), skip = 3) %>% 
  filter(space)

# list to store the data
ved <- list(bin_vars = list(count_mat = NA,
                            bin_mat = NA,
                            complete_mat = NA,
                            over5 = NA), # variables (grave goods) present in more than 5% of graves
            cont_vars = NA,
            cat_vars = NA,
            metadata = NA,
            layout = NA,
            var_names = list(full = NA,
                             short = NA,
                             long = NA),
            id_burials = NA,
            orig_dataset = input)

# filling the data -------------------------------------------------------------

# # burial id
ved$id_burials <- input %>%
  filter(analysis == TRUE) %>%
  pull(id_burial)

# binary variables - artefact coocurences
ved$bin_vars$count_mat <- input %>% 
  filter(analysis == TRUE) %>% 
  select(pigment, pol_adze, pol_axe, grinding, pebble, 
         lit_local, lit_nonlocal, 
         pot_special, pot_bowl, pot_globular, pot_bottle, 
         # pot_head, pot_mid, pot_leg, 
         bone_tool, # bone_awl,
         pendant, # pendant_L, pendant_U, pendant_I, 
         buckle_O, 
         bracelet_spond, beads_spond, beads_marble, 
         deer_teeth, antler, shell) %>% 
  as.matrix()

rownames(ved$bin_vars$count_mat) <- ved$id_burials

ved$bin_vars$bin_mat <- apply(ved$bin_vars$count_mat, 2, binarize)

# complete indicator matrix
ved$bin_vars$complete_mat <- ved$bin_vars$bin_mat %>% 
  as_tibble() %>% 
  mutate(across(everything()), abs(. - 1)) %>% 
  select_all(list(~ paste0(., ".neg"))) %>% 
  bind_cols(as_tibble(ved$bin_vars$bin_mat)) %>% 
  select(order(colnames(.))) %>% as.matrix()

rownames(ved$bin_vars$complete_mat) <- ved$id_burials

# grave goods present in more than 5% of graves
n_burs <- nrow(filter(input, analysis))
ved$bin_vars$over5 <- (colSums(ved$bin_vars$bin_mat) / n_burs) > 0.05
table(ved$bin_vars$over5)

# continuous variables
ved$cont_vars <- input %>% 
  filter(analysis == TRUE) %>% 
  select(
    # pit_len, pit_wid, pit_dep,
    # d13c, d15n, 
    sr, sr_ppm, 
    # body_height
  ) %>% 
  as.data.frame()

rownames(ved$cont_vars) <- ved$id_burials

# factor variables
# ved$cat_vars <- input %>% 
#   filter(undisturbed == TRUE) %>% 
#   select(pit_orient_cat, head_orient_cat, body_side) %>%
#   mutate_all(factor) %>% 
#   as.data.frame()
# 
# rownames(ved$cat_vars) <- ved$id_burials

# metadata
# lookup_sex <- c("n. a." = "gold", "F" = "tomato", "M" = "steelblue")

ved$metadata <- input %>% 
  select(id_burial, analysis, sex, age_cat, dat) %>% 
  mutate(id_burial = as.character(id_burial),
         pres = if_else(analysis, "pres.", "dist."),
         sex = factor(sex),
         # sex_col = unname(lookup_sex[input$sex]),
         age_cat = ordered(age_cat, 
                           levels = c("I1", "I", "I2", "J",
                                      "A1", "A", "A2", "AM", 
                                      "M1", "M", "M2", "S")),
         age_sim = ordered(if_else(age_cat %in% c("I1", "I", "I2", "J"), "juv.",
                                   if_else(age_cat %in% c("A1", "A", "A2"), "ad.",
                                           if_else(age_cat %in% c("AM", "M1", "M", 
                                                                  "M2", "S"), "mat.", 
                                                   "ind."))), 
                           levels = c("juv.", "ad.", "mat.", "ind.")),
         dat = ordered(dat, levels = c("un", "Ib", "IIa")))

# layout
ved$layout <- input %>% 
  select(layout_x, layout_y) %>% 
  as.matrix()

rownames(ved$layout) <- input$id_burial

# variable names ----------------------------------------------------------

ved$var_names$full <- tibble(vnames = names(input),
                             abbrv = c("id", "row", "id.unified", "id.full", "analysis", "space", "note",
                                       "x", "y", "sex", "sex.d", "sex.p",
                                       "age.cat", "age.int", "pit.l", "pit.w", "pit.d",
                                       "pit.o", "pit.o.cat", "head.o", "head.o.cat",
                                       "body.side",
                                       "pigm", "pol.sum", "adz", "axe", "grnd", "peb",
                                       "l.sum", "l.loc", "l.non", "l.kl", "l.skj", "l.oth", "l.rad", "l.sgs",
                                       "p.sum", "p.bow", "p.glo", "p.bot", "p.oth", "p.frag",
                                       "p.top", "p.mid", "p.leg","p.fill",
                                       "bon.t", "bon.awl", "pen", "pen.L", "pen.U", "pen.I", 
                                       "buck", "brac", "neck", "b.sp", "b.mar", "b.div",
                                       "b.cyl", "b.circ", "b.olive", "d.tooth", "ant", "shell",
                                       "dat", "d13C", "d15N", "Sr", "Sr.ppm", "body.h"),
                             long = c("Burial ID", "row", "id.unified", "id.full", "analysis", "space", "note",
                                      "x", "y", "sex", "sex.d", "sex.p",
                                      "Age", "age.int", "Pit length", "Pit width", "Pit depth",
                                      "Pit orientation", "Pit orientation", "Head orientation", "Head orientation",
                                      "Side",
                                      "Pigment", "PST (sum)", "Adze", "Axe", "Grinding tool", "Pebble",
                                      "Lithics (sum)", "Local lithics", "Non-local lithics", "Lithics (KF)", 
                                      "Lithics (SKJ)", "Lithics (other)", "Lithics (Rad.)", "Lithics (SGS)",
                                      "Pottery (sum)", "Bowl", "Globular pot", "Bottle", "Other pottery", "Pottery fragment",
                                      "Pot around head", "Pot in the middle ", "Pot around legs","Pot in fill",
                                      "Bone tool", "Bone awl", "Pendant", "L-shaped pendant", "U-shaped pendant", "I-shaped pendant", 
                                      "O-shaped Buckle", "Bracelet", "Necklace", "Spondylus bead", "Marble bead", "Dividing bead",
                                      "Cylinder bead", "Circular bead", "Olive-shaped bead", "Deer tooth", "Antler", "Shell",
                                      "Relative chrono.", "d13C", "d15N", "Sr", "Sr.ppm", "Body height"))

ved$var_names$short <- ved$var_names$full$abbrv
names(ved$var_names$short) <- ved$var_names$full$vnames

ved$var_names$long <- ved$var_names$full$long
names(ved$var_names$long) <- ved$var_names$full$vnames


# origin ------------------------------------------------------------------

local_range <- left_join(ved$metadata, 
                         as_tibble(ved$cont_vars, rownames = "id_burial"),
                         by = "id_burial") %>% 
  filter(age_cat %in% c("I", "I1", "I2"), !is.na(sr)) %>% 
  summarise(n = n(), 
            mean = mean(sr, na.rm = TRUE), 
            sd = sd(sr, na.rm = TRUE)) %>% 
  mutate(min = mean - 2 * sd,
         max = mean + 2 * sd)

ved$metadata <- left_join(ved$metadata, 
                          as_tibble(ved$cont_vars, rownames = "id_burial"),
                          by = "id_burial") %>% 
  mutate(origin = if_else(sr > local_range$min & sr < local_range$max, 
                          "local", "non-local"),
         origin = if_else(is.na(origin), "ind.", origin),
         origin = factor(origin, levels = c("local", "non-local", "ind."))) %>% 
  select(-sr, -sr_ppm)


# # mixed categories --------------------------------------------------------
# detect_ind <- function(x) {
#   if_else(str_detect(x, "ind."), "ind.", x)
# }
# 
# # only categories that are common for at least 5% of burials are kept
# threshold <- 0.05 * nrow(ved$bin_vars$count_mat)
# 
# vec_origin <- c("local" = "loc.", "non-local" ="non.", "ind." = "ind.")
# 
# comb_cats <- ved$metadata %>% 
#   filter(pres != "dist.") %>% 
#   mutate(
#     origin = vec_origin[origin],
#     sex = str_remove(sex, "\\s"),
#     cat_sa = str_c(sex, age_sim, sep = "/"),
#     cat_os = str_c(origin, sex, sep = "/"),
#     cat_oa = str_c(origin, age_sim, sep = "/"),
#     cat_osa = str_c(origin, sex, age_sim, sep = "/"),
#     across(starts_with("cat_"), detect_ind)
#   ) %>% 
#   select(id_burial, starts_with("cat"))
# 
# get_perv_lvls <- function(x) {
#   lvls <- table(x)
#   names(lvls[lvls >= threshold])
# }
# 
# comb_cats_over5 <- comb_cats %>% mutate(
#   cat_sa = if_else(cat_sa %in% get_perv_lvls(cat_sa), cat_sa, "ind."),
#   cat_os = if_else(cat_os %in% get_perv_lvls(cat_os), cat_os, "ind."),
#   cat_oa = if_else(cat_oa %in% get_perv_lvls(cat_oa), cat_oa, "ind."),
#   cat_osa = if_else(cat_osa %in% get_perv_lvls(cat_osa), cat_osa, "ind.")
# )
# 
# ved$metadata <- left_join(ved$metadata, comb_cats_over5, by = c("id_burial")) %>%
#   mutate(across(starts_with("cat_"), function(x) if_else(is.na(x), "ind.", x)))

# saving processed dataset =====================================================
if (!dir.exists(here("data", "temp"))) {
  dir.create(here("data", "temp"))
}

saveRDS(ved, here("data/temp", "vedrovice_dataset.RDS"))

# ved <- readRDS(here("data/temp", "vedrovice_dataset.RDS"))
