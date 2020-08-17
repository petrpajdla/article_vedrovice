# Project "Vedrovice"
# Script nr. 00
# KNIT OUTPUT FILE
# author: Petr Pajdla
# Create the text of the article...

# knitr::knit(input = here::here("text", "article-vedrovice.Rmd"),  
#             output = here::here("text", "article-vedrovice.md"))

# knit main text
rmarkdown::render(input = here::here("text", "article_vedrovice.Rmd"), 
                  output_format = bookdown::pdf_document2())
rmarkdown::render(input = here::here("text", "article_vedrovice.Rmd"), 
                  output_format = bookdown::word_document2())

# knit si
rmarkdown::render(input = here::here("text", "si_vedrovice.Rmd"), 
                  output_format = bookdown::pdf_document2())

