# recoding data, recode and reimport to spreadsheet...

input <- read_csv("./input/data_vedrovice.csv")

# recoding head_orient =========================================================
# output <- input %>% transmute(id_burial, 
#                     head_orient_new = ifelse(input$head_orient > 90 - (45/2) & input$head_orient < 90 + (45/2), "E", 
#                                           ifelse(input$head_orient > 135 - (45/2) & input$head_orient < 135 + (45/2), "SE", 
#                                                  ifelse(input$head_orient > 180 - (45/2) & input$head_orient < 180 + (45/2), "S", 
#                                                         ifelse(input$head_orient > 315 - (45/2) - 2 & input$head_orient < 315 + (45/2), "NW", "NOPE")))))
# 
# write_csv(output, "./output/head_orient.csv")

# recoding pit_orient ==========================================================
# # CHECK DATA !!! 
# input %>% group_by(pit_orient_cat) %>% summarise(min = min(pit_orient, na.rm = T), max = max(pit_orient, na.rm = T), n())
# hist(input$pit_orient)
# ggplot(input, aes(pit_orient, factor(pit_orient_cat))) + geom_jitter()
# ggplot(input, aes(pit_orient)) + geom_histogram(binwidth = 2)

# binarizing values ============================================================

# binary <- biclust::binarize(input[, 14:36], threshold = 0)
# 
# write_csv(binary, "./output/binary.csv")

