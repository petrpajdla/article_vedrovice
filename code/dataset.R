# Project "Vedrovice"
# Script nr. 1
# DATA INPUT AND MANIPULATION
# author: Petr Pajdla
# Script manipulates input data and saves a list with manipulated dataset

# packages =====================================================================
library(dplyr)
library(here)

# functions ====================================================================
# continuous vector to dichotomous (binarization)
binarize <- function(x, threshold = 0) {
  x[x > threshold] <- 1
  x[x <= threshold] <- 0
  return(x)
}

# data input and manipulation ==================================================
input <- readr::read_csv(here("data", "data_vedrovice.csv"), skip = 3) %>% 
  filter(undisturbed == TRUE)

# list to store the data
ved <- list(bin_vars = list(count_mat = NA,
                            bin_mat = NA),
            cont_vars = list(cont_vars = NA),
            cat_vars = list(cat_vars = NA),
            metadata = NA,
            layout = NA,
            var_names = list(full = NA,
                             short = NA),
            id_burials = input$id_burial,
            orig_dataset = input)

# filling the data -------------------------------------------------------------
# binary variables - artefact coocurences
ved$bin_vars$count_mat <- input %>% 
  select(pigment, pol_adze, pol_axe, grinding, pebble, 
         lit_other, lit_kl, lit_skj, 
         pot_special, pot_bowl, pot_globular, pot_bottle, 
         pot_head, pot_mid, pot_leg, bone_tool, 
         pendant_L, pendant_U, pendant_I, buckle_O, 
         bracelet_spond, beads_spond, beads_marble) %>% 
  as.matrix()

rownames(ved$bin_vars$count_mat) <- ved$id_burials

ved$bin_vars$bin_mat <- apply(ved$bin_vars$count_mat, 2, binarize)

# continuous variables
ved$cont_vars$cont_vars <-  input %>% select(pit_len, pit_wid, pit_dep,
                                             d13c, d15n, sr, sr_ppm, body_height) %>% 
  as.data.frame()

rownames(ved$cont_vars$cont_vars) <- ved$id_burials

# factor variables
ved$cat_vars$cat_vars <- input %>% 
  select(pit_orient_cat, head_orient_cat, body_side) %>%
  mutate_all(factor) %>% 
  as.data.frame()

rownames(ved$cat_vars$cat_vars) <- ved$id_burials

# metadata
lookup_sex <- c("C" = "gold", "F" = "tomato", "M" = "steelblue")

ved$metadata <- input %>% select(id_burial, sex, age_cat, dat) %>% 
  mutate(sex = factor(sex),
         sex_col = unname(lookup_sex[input$sex]),
         age_cat = ordered(age_cat, 
                           levels = c("I1", "I", "I2", "J",
                                      "A1", "A", "A2", "AM", 
                                      "M1", "M", "M2", "S")),
         age_sim = ordered(if_else(age_cat %in% c("I1", "I", "I2", "J"), "you",
                                   if_else(age_cat %in% c("A1", "A", "A2"), "mid",
                                           if_else(age_cat %in% c("AM", "M1", "M", 
                                                                  "M2", "S"), "old", 
                                                   "NA"))), 
                           levels = c("you", "mid", "old")),
         dat = ordered(dat, levels = c("un", "Ib", "IIa")))

# layout
ved$layout <- input %>% select(layout_x, layout_y) %>% 
  as.matrix()

rownames(ved$layout) <- input$id_burial

ved$var_names$full <- tibble(vnames = names(input),
                        abbrv = c("id", "row", "id.unified", "id.full", "undist",
                                  "x", "y", "sex", "sex.d", "sex.p",
                                  "age.cat", "age.int", "pit.l", "pit.w", "pit.d",
                                  "pit.o", "pit.o.cat", "head.o", "head.o.cat",
                                  "body.side",
                                  "pigm", "pol.sum", "adz", "axe", "grd", "peb",
                                  "l.sum", "l.kl", "l.skj", "l.oth", "l.rad", "l.sgs",
                                  "p.sum", "p.bow", "p.glo", "p.bot", "p.spec", "p.frag",
                                  "p.top", "p.mid", "p.leg","p.fill",
                                  "bon.t", "bon.awl", "pen", "pen.L", "pen.U", "pen.I", 
                                  "buck", "brac", "neck", "b.sp", "b.mar", "b.div",
                                  "b.cyl", "b.circ", "b.olive", "d.tooth", "ant", "shell",
                                  "dat", "d13C", "d15N", "Sr", "Sr.ppm", "body.h"))

ved$var_names$short <- ved$var_names$full$abbrv
names(ved$var_names$short) <- ved$var_names$full$vnames

# saving processed dataset =====================================================
saveRDS(ved, here("data/temp", "vedrovice_dataset.RDS"))

# ved <- readRDS(here("data/temp", "vedrovice_dataset.RDS"))
