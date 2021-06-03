# Project "Vedrovice"
# Script nr. 00
# KNIT OUTPUT FILE
# author: Petr Pajdla
# Create the text of the article in DOCX for further processing

# knit main text
rmarkdown::render(input = here::here("text", "article_vedrovice.Rmd"), 
                  output_format = bookdown::word_document2(
                    toc = FALSE,
                    fig_caption = TRUE,
                    number_sections = FALSE
                  ))