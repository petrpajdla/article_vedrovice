# Project "Vedrovice"
# Script nr. 00
# KNIT OUTPUT FILE 
# author: Petr Pajdla
# Create the SI

# knit si
rmarkdown::render(input = here::here("text", "si_vedrovice.Rmd"),
                  output_format = bookdown::pdf_document2(
                    toc = FALSE,
                    latex = "pdflatex",
                    fig_caption = TRUE,
                    number_sections = FALSE,
                    keep_tex = FALSE
                  ))

