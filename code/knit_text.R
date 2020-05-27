# Project "Vedrovice"
# Script nr. 00
# KNIT OUTPUT FILE
# author: Petr Pajdla
# Create the text of the article...

# knitr::knit(input = here::here("text", "article-vedrovice.Rmd"),  
#             output = here::here("text", "article-vedrovice.md"))

rmarkdown::render(input = here::here("text", "article_vedrovice.Rmd"))
