ENSEMBL = KH-ENS.blast
ZIP = 'index.html?file=data%2Fcirobu%2FKH-ENS.blast.zip'
URL = https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fcirobu%2FKH-ENS.blast.zip

install: build
	R CMD build build
	R CMD INSTALL build

build: data
	mkdir -p build
	cp -r R build
	cp -r data build
	cp -r man build
	cp DESCRIPTION build
	cp NAMESPACE build

data: data-raw/cisbp_orthology.txt
	Rscript --vanilla data-raw/DATASET.R

data-raw/cisbp_orthology.txt: data-raw/$(ENSEMBL)
	Rscript --vanilla data-raw/getOrthologs.R

data-raw/$(ENSEMBL):
	wget $(URL)
	unzip -o $(ZIP)
	mv $(ENSEMBL) data-raw

clean:
	rm -f $(ZIP)
	rm -rf build
