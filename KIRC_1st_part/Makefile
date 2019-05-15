all : index.html

index.html : index.Rmd
	Rscript -e "rmarkdown::render('$<')"
