all: text/article_vedrovice.pdf

text/article_vedrovice.pdf: text/article_vedrovice.Rmd text/bib_vedrovice.bib code/knit_text.R
	cd code; R CMD BATCH knit_text.R

.PHONY: clean
clean:
	cd code; rm Rplots.pdf .RData
