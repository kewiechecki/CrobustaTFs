ENSEMBL = KH-ENS.blast
ZIP = 'index.html?file=data%2Fcirobu%2FKH-ENS.blast.zip'
URL = https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fcirobu%2FKH-ENS.blast.zip

install: data
	R CMD INSTALL .

data: data-raw/cisbp_orthology.txt
	Rscript --vanilla data-raw/DATASET.R

data-raw/cisbp_orthology.txt: data-raw/$(ENSEMBL)
	Rscript --vanilla data-raw/getOrthologs.R

data-raw/$(ENSEMBL): $(ZIP)
	unzip -o $(ZIP)
	rm -f $(ZIP).zip
	mv $(ENSEMBL) data-raw

$(ZIP):
	wget $(URL)
