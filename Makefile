install: data
	R CMD INSTALL .

data: data-raw
	Rscript --vanilla data-raw/DATASET.R

