# Project "Vedrovice"
# Script nr. 4
# DISTANCES BETWEEN BURIALS
# author: Petr Pajdla
# what is here?

set.seed(42)

# packages =====================================================================
library(dplyr)
library(cluster)
library(igraph)

# functions ====================================================================
# continuous vector to dichotomous (binarization)
binarize <- function(x, threshold = 0) {
  x[x > threshold] <- 1
  x[x <= threshold] <- 0
  return(x)
}

# replace missing (NA) data with mean for the vector
replace_na <- function(x) {
  avg <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- avg
  return(x)
}

# perform chi-squre test on multiple columns of a data frame df
# colpairs contains pairs of column name combinations that chi-squared test is 
# performed on, created by grid.expand() on a vector of column names
multiple.chisq.test <- function(df, colpairs) {
  p_val <- data.frame(colpairs[, 1:2], 
                      rep(NA, nrow(colpairs)),
                      rep(NA, nrow(colpairs)))
  colnames(p_val) <- c("var1", "var2", "p.val", "independent")
  for (i in 1:nrow(colpairs)) {
    cont_tab <- table(df[, colpairs[i, 1:2]])
    p_val[i, "p.val"] <- round(chisq.test(cont_tab)$p.value, 3)
    if (p_val[i, "p.val"] <= 0.05) {
      p_val[i, "independent"] <- FALSE
    } else {
      p_val[i, "independent"] <- TRUE
    }
  }
  return(p_val)
}

# count v values for a data frame (based on articles by Manly)
# where: R...nr of variables, df...data frame of coocurences, 
# expected = df of expected  coocurences (med.coocs)
get_v <- function(df, expected = cooc_exp, R = length(vars_binary)) {
  v_vals <- rep(NA, R)
  w_df <- left_join(df, expected, by = c("var1", "var2"))
  v_vals <- w_df %>% group_by(var1) %>% 
    summarise(v = sum((nr.cooc - mean.cooc)^2 / R))
  return(v_vals)
}

# expand grid function to contain only unique combinations of vector values
# expand.grid.unique <- function(x, y, include.equals=FALSE) {
#   x <- unique(x)
#   y <- unique(y)
#   g <- function(i) {
#     z <- setdiff(y, x[seq_len(i-include.equals)])
#     if(length(z)) cbind(x[i], z, deparse.level=0)
#   }
#   do.call(rbind, lapply(seq_along(x), g))
# }

# counting coocurences in an occurence matrix
count_coocurence <- function(matrix, colpairs) {
  cooc <- data.frame(colpairs[, 1:2],
                     rep(0, nrow(colpairs)))
  names(cooc) <- c("var1", "var2", "nr.cooc")
  for (i in 1:nrow(colpairs)) {
    cooc[i, "nr.cooc"] <- sum(unname(matrix[, unname(colpairs[i, 1])] == 1) & 
                                unname(matrix[, unname(colpairs[i, 2])] == 1))
  }
  return(cooc)
}

# data input and manipulation ==================================================
input <- readr::read_csv("./data/data_vedrovice_v2.csv", skip = 3) %>% 
  filter(undisturbed == TRUE)
# names(input)

# selecting relevant attributes
ved <- input %>% select(pit_len, pit_wid, pit_dep, 
                        pit_orient_cat, head_orient_cat, body_side,
                        pigment, pol_adze, pol_axe, grinding, pebble, 
                        lit_other, lit_kl, lit_skj, 
                        pot_special, pot_bowl, pot_globular, pot_bottle, 
                        pot_head, pot_mid, pot_leg, bone_tool, 
                        L_pendant, U_pendant, I_pendant, O_buckle, 
                        spond_bracelet, spond_beads, marble_beads, 
                        d13c, d15n, sr, body_height) %>%
  mutate_at(c("pit_orient_cat", "head_orient_cat", "body_side"), factor) %>% 
  as.data.frame()

rownames(ved) <- input$id_burial

# variable types for subsetting ------------------------------------------------
vars_binary <- c("pigment", "pol_adze", "pol_axe", "grinding", "pebble", 
                 "lit_other", "lit_kl", "lit_skj", 
                 "pot_special", "pot_bowl", "pot_globular", "pot_bottle", 
                 "pot_head", "pot_mid", "pot_leg", "bone_tool", 
                 "L_pendant", "U_pendant", "I_pendant", "O_buckle", 
                 "spond_bracelet", "spond_beads", "marble_beads")
vars_cont <- c("pit_len", "pit_wid", "pit_dep", "body_height", 
               "d13c", "d15n", "sr")

# data frame with binarized variables ------------------------------------------
# for distance calculation, for other reasons, non binarized variant is used
binarized <- ved %>% mutate_at(vars_binary, binarize)
rownames(binarized) <- rownames(ved)

# identifying dependent variables  =============================================

# This part partly answers Q1

# correlation of continuous variables ------------------------------------------
psych::pairs.panels(ved[, vars_cont], lm = TRUE)
corrplot::corrplot(cor(ved[, vars_cont], use = "p"), 
                   type = "lower", diag = FALSE)
# only some mild correlations, no probs...

# chi-squared test of categorical and binary variables -------------------------
# H0: x1 independent on x2 with p = 0.05, if p < 0.05, H0 rejected
chisq.test(table(ved$pit_orient_cat, ved$head_orient_cat)) # rejected = dep
chisq.test(table(ved$pit_orient_cat, ved$body_side)) # not rejected = indep
chisq.test(table(ved$body_side, ved$head_orient_cat)) # not rejected = indep
# pit_orient and head_orient are dependent, thus only one will be used further

# colpairs <- expand.grid.unique(vars_binary, vars_binary)
colpairs <- as.matrix(expand.grid(c(vars_binary, "pit_orient_cat", "body_side"), 
                              c(vars_binary, "pit_orient_cat", "body_side")))

# variables that are dependent with p = 0.05
var_dependent <- multiple.chisq.test(binarized, colpairs)

var_dependent <- var_dependent %>% filter(var1 != var2, 
                                     independent == FALSE)

# network of dependent variables -----------------------------------------------
g_var_dependent <- simplify(graph.data.frame(d = var_dependent, 
                                                 directed = FALSE))

var_labels <- data.frame(name = unique(colpairs[, 1]), 
                         abbrv = c("pig",
                                   "adz",
                                   "ax",
                                   "grd",
                                   "peb",
                                   "l.oth",
                                   "l.kl",
                                   "l.skj",
                                   "p.oth",
                                   "p.bow",
                                   "p.glo",
                                   "p.bot",
                                   "p.top",
                                   "p.mid",
                                   "p.leg",
                                   "bon.t",
                                   "pen.L",
                                   "pen.U",
                                   "pen.I",
                                   "buck",
                                   "brac",
                                   "b.sp",
                                   "b.mar",
                                   "pit.o",
                                   "side"),
                         col = c("gray80", 
                                 rep("deepskyblue", 7), 
                                 rep("greenyellow", 7),
                                 "gray80",
                                 rep("gold", 7),
                                 rep("gray80", 2)))

rownames(var_labels) <- var_labels$name

V(g_var_dependent)$abbrv <- as.character(var_labels[
  V(g_var_dependent)$name, "abbrv"])
V(g_var_dependent)$col <- as.character(var_labels[
  V(g_var_dependent)$name, "col"])
V(g_var_dependent)$size <- unname(eigen_centrality(g_var_dependent)$vector*20)
V(g_var_dependent)$shp <- if_else(V(g_var_dependent)$size > 5, "circle", "none")

lay_vars <- layout_with_fr(g_var_dependent)

# plot(g_var_dependent,
#      layout = lay_vars,
#      vertex.shape = "rectangle",
#      vertex.size = 12, vertex.size2 = 6,
#      vertex.color = V(g_var_dependent)$col,
#      vertex.frame.color = V(g_var_dependent)$col,
#      vertex.label = V(g_var_dependent)$abbrv,
#      vertex.label.family = "mono",
#      vertex.label.color = "black",
#      vertex.label.cex = 0.8,
#      edge.width = 0.6,
#      mark.groups = cluster_edge_betweenness(g_var_dependent),
#      mark.col = NA)

# pdf("./plots/dependent_variables.pdf", width = 14, height = 7)
# par(mfrow = c(1, 2), mar = c(rep(1, 4)))
# plot(g_var_dependent,
#      layout = lay_vars,
#      vertex.shape = "none",
#      vertex.label = V(g_var_dependent)$abbrv,
#      vertex.label.family = "mono",
#      vertex.label.cex = 1,
#      vertex.label.color = "black",
#      vertex.color = NA,
#      edge.width = 0.8,
#      mark.groups = cluster_edge_betweenness(g_var_dependent),
#      mark.col = NA,
#      mark.border = "gray80")
# plot(g_var_dependent,
#      layout = lay_vars,
#      vertex.shape = "circle",
#      vertex.label = NA,
#      vertex.label.family = "mono",
#      vertex.label.cex = 0.8,
#      vertex.label.color = "black",
#      vertex.size = eigen_centrality(g_var_dependent)$vector*20,
#      vertex.color = NA,
#      edge.width = 0.8)
# graphics.off()

pdf("./plots/dependent_variables.pdf", width = 20, height = 20)
par(mar = rep(0, 4))
plot(g_var_dependent,
     layout = lay_vars,
     vertex.shape = V(g_var_dependent)$shp,
     vertex.label = V(g_var_dependent)$abbrv,
     vertex.label.family = "mono",
     vertex.label.cex = 2,2,
     vertex.label.color = "black",
     vertex.size = V(g_var_dependent)$size,
     vertex.color = "white",
     edge.width = 0.8,
     edge.color = "gray60",
     mark.groups = cluster_edge_betweenness(g_var_dependent),
     mark.col = "gray90",
     mark.border = NA,
     mark.expand = 14)
dev.off()

rm(list = c("colpairs", "g_var_dependent", "lay_vars", "var_dependent"))

# removing dependent variables -------------------------------------------------
ved <- ved %>% select(-head_orient_cat)
binarized <- binarized %>% select(-head_orient_cat)



# metadata tables ==============================================================
lookup_sex <- c("C" = "gold", "F" = "tomato", "M" = "steelblue")

meta <- input %>% select(id_burial, sex, age_cat, dat) %>% 
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

meta[is.na(meta$sex_col), "sex_col"] <- "gray50"

meta$ei <- ei

# spatial layout of x~y
lay_space <- input %>% select(layout_x, layout_y) %>% 
  as.matrix()

rownames(lay_space) <- input$id_burial

# plot(lay_space, col = meta$sex_col, pch = 20, cex = meta$ei*4, asp = c(1,1),
#      ann = F, axes = F)
# # spatial distance matrix
# d <- dist(lay_space, method = "euclidean")
# 
# plot(meta$ei ~ meta$sex)
# plot(meta$ei ~ factor(meta$dat))
# plot(meta$ei ~ meta$age_sim)
# plot(meta$ei ~ meta$age_cat)

rm(list = c("for_pca", "ei", "lookup_sex"))

# randomization of grave goods =================================================
colpairs2 <- as.matrix(expand.grid(c(vars_binary), 
                                   c(vars_binary)))

# original observed matrix
orig_mat <- as.matrix(binarized %>% select(vars_binary))

# randomization of occurence matrix - list of many matrices
n_permutations = 10000
rand_mat <- vegan::permatfull(orig_mat, fixedmar = "both",
                              mtype = "prab", times = n_permutations)$perm

# count coocurences ------------------------------------------------------------
# coocurences in an observed matrix
cooc_obs <- count_coocurence(matrix = orig_mat, colpairs = colpairs2)

# detecting cores for parallel
no_cores <- parallel::detectCores() - 1

# # ---------------------------------
# cooccurences on random matrices
# going parallel, works on linux only I presume
cl <- parallel::makeCluster(no_cores)
cooc_rand <- parallel::parLapply(cl, rand_mat, count_coocurence, colpairs2)
parallel::stopCluster(cl)
# 
# save(cooc_rand, file = "./data/temp/coocurence_random_matrix.Rdata")
# # ---------------------------------
# cooc_rand <- load("./data/temp/coocurence_random_matrix.Rdata", verbose = TRUE)

# summarization of random occurences into expected cooccurences
cooc_exp <- cooc_rand %>% bind_rows() %>% 
  group_by(var1, var2) %>% 
  summarise(mean.cooc = mean(nr.cooc),
            med.cooc = median(nr.cooc))

v_observed <- tidyr::spread(get_v(cooc_obs), var1, v)

# creating cluster with forked global environment
cl <- parallel::makeCluster(no_cores - 2, type = "FORK")
v_experimental <- parallel::parLapply(cl, cooc_rand, get_v)
v_experimental <- bind_rows(parallel::parLapply(cl, v_experimental, 
                                                tidyr::spread, "var1", "v"))
parallel::stopCluster(cl)

hist(v_experimental$pol_adze)
abline(v = v_observed$pol_adze)

hist(v_experimental$bone_tool)
abline(v = v_observed$bone_tool)

# p-val: nr of occurences when vi is equall or exceeds experimental value
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
v_observed[1] >= v_experimental$bone_tool

10000/(sum(v_experimental$bone_tool <= v_observed[[1]])*100)

get_P <- function(v_exp, v_obs = v_observed, size = 
                    n_permutations) {
  size / (sum(v_exp <= v_obs) * 100)
}

median(v_experimental$pol_adze)
v_observed$pol_adze


# s statistic for whole dataset????
s_stat <- c("s_obs" = NA, "s_exp" = NA)
s_stat["s_obs"] <- sum(v_observed / length(vars_binary))
median(v_experimental / length(vars_binary))

# table of variables for the article ===========================================

bind_cols(variable = vars_binary, 
      abbreviation = as.character(var_labels[, "abbrv"])[1:23],
      description = c(
        "presence of ochre pigment",
        "polished stone adze",
        "polished stone axe",
        "ground stone tool",
        "hammerstone or pebble",
        "lithics, various raw materials",
        "lithics, Krumlov forest chers",
        "lithics, Krakow-Czestochowa flint",
        "other/special shapes of pottery",
        "pottery - bowl",
        "pottery - globular shape",
        "pottery - bottle",
        "pottery positioned around head",
        "pottery in the middle of the grave",
        "pottery close to legs",
        "bone tools",
        "L-shaped pendant",
        "U-shaped pendant",
        "I-shaped pendant",
        "O-shaped buckle",
        "bracelet made of spondylus shell",
        "beads made of spondylus shell",
        "beads made of marble"),
      "v stat. exp." = rep(NA, length(vars_binary)),
      "v stat. obs." = )

# counting distance matrix =====================================================

var_weights <- c(rep(1/3, 3), rep(1/2, 2), rep(1, 28))

# result of daisy is a dissimilarity (distance) matrix
# i.e. 0 = similar vs 1 = dissimilar
dist <- daisy(binarized, metric = "gower", 
              type = list(asymm = vars_binary), weights = var_weights)

summary(dist)
# hist(dist, xlim = c(0, 1), breaks = 10)

# ???????????????????????, corrplot
corrplot::corrplot(as.matrix(dist), 
                   type = "upper", order = "FPC", is.corr = FALSE)

# melting distance matrix into a df of individual edges ------------------------
dist_mx <- as.matrix(dist)
dist_mx[upper.tri(dist_mx)] <- NA

dist_df <- reshape2::melt(dist_mx)
names(dist_df) <- c("from", "to", "distance")

dist_df <- dist_df %>% filter(from != to, !is.na(distance))

# hist(dist_df$distance, breaks = 10)
quantile(dist_df$distance)[2]

# cutoff value to create edges -------------------------------------------------
dist_df %>% filter(distance < quantile(dist_df$distance)[2])

plot(graph.data.frame(d = dist_df %>% filter(distance < quantile(dist_df$distance)[2]), 
                 vertices = meta, 
                 directed = FALSE), vertex.label = NA, vertex.size = (meta$ei)*10, vertex.color = meta$sex_col)




