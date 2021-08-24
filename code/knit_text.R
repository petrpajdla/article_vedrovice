# Project "Vedrovice"
# Script nr. 00
# KNIT OUTPUT FILE
# author: Petr Pajdla
# Create the text of the article in pdf for comments etc.

# knit main text
rmarkdown::render(input = here::here("text", "article_vedrovice.Rmd"), 
                  output_format = bookdown::pdf_document2(
                    toc = FALSE,
                    latex = "pdflatex",
                    fig_caption = TRUE,
                    number_sections = TRUE,
                    global_numbering = FALSE,
                    keep_tex = FALSE
                  ))

