all: text/article_vedrovice.pdf text/article_vedrovice.docx

text/article_vedrovice.pdf: text/article_vedrovice.Rmd text/bib_vedrovice.bib code/knit_text.R
	cd code; R CMD BATCH knit_text.R

text/article_vedrovice.docx: text/article_vedrovice.Rmd text/bib_vedrovice.bib code/knit_text_word.R
	cd code; R CMD BATCH knit_text_word.R

.PHONY: clean
clean:
	cd code; rm Rplots.pdf .RData
