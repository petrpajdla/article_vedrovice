# Project "Vedrovice"
# Script nr. 1.1
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
input <- readr::read_csv(here("data", "data_vedrovice_v01.csv"), skip = 3) %>% 
  filter(lbk == TRUE)

# list to store the data
ved <- list(bin_vars = list(count_mat = NA,
                            bin_mat = NA),
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
  filter(undisturbed == TRUE) %>%
  pull(id_burial)

# binary variables - artefact coocurences
ved$bin_vars$count_mat <- input %>% 
  filter(undisturbed == TRUE) %>% 
  select(pigment, pol_adze, pol_axe, grinding, pebble, 
         lit_local, lit_nonlocal, 
         pot_special, pot_bowl, pot_globular, pot_bottle, 
         # pot_head, pot_mid, pot_leg, 
         bone_tool, 
         pendant_L, pendant_U, pendant_I, buckle_O, 
         bracelet_spond, beads_spond, beads_marble, antler) %>% 
  as.matrix()

rownames(ved$bin_vars$count_mat) <- ved$id_burials

ved$bin_vars$bin_mat <- apply(ved$bin_vars$count_mat, 2, binarize)

# continuous variables
ved$cont_vars <- input %>% 
  filter(undisturbed == TRUE) %>% 
  select(pit_len, pit_wid, pit_dep,
         d13c, d15n, sr, sr_ppm, 
         body_height) %>% 
  as.data.frame()

rownames(ved$cont_vars) <- ved$id_burials

# factor variables
ved$cat_vars <- input %>% 
  filter(undisturbed == TRUE) %>% 
  select(pit_orient_cat, head_orient_cat, body_side) %>%
  mutate_all(factor) %>% 
  as.data.frame()

rownames(ved$cat_vars) <- ved$id_burials

# metadata
lookup_sex <- c("n. a." = "gold", "F" = "tomato", "M" = "steelblue")

ved$metadata <- input %>% select(id_burial, sex, age_cat, dat) %>% 
  mutate(id_burial = as.character(id_burial),
         sex = factor(sex),
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

# variable names ----------------------------------------------------------

ved$var_names$full <- tibble(vnames = names(input),
                             abbrv = c("id", "row", "id.unified", "id.full", "undist", "lbk",
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
                             long = c("Burial ID", "row", "id.unified", "id.full", "undist", "lbk",
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

# saving processed dataset =====================================================
saveRDS(ved, here("data/temp", "vedrovice_dataset.RDS"))

# ved <- readRDS(here("data/temp", "vedrovice_dataset.RDS"))
