all: hw3.pdf hw3_solns.html hw3_solns.pdf

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

hw3_solns.html: hw3_solns.Rmd
	R ${R_OPTS} -e 'library(knitr);knit2html("hw3_solns.Rmd")'

hw3_solns.pdf: hw3_solns.Rmd
	pandoc --smart hw3_solns.md -o hw3_solns.pdf

hw3.pdf: hw3.tex
	pdflatex hw3.tex
