ENSEMBL = KH-ENS.blast
ZIP = 'index.html?file=data%2Fcirobu%2FKH-ENS.blast.zip'
URL = https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fcirobu%2FKH-ENS.blast.zip

install: CrobustaTFs_0.0.3.tar.gz
	R CMD INSTALL build

CrobustaTFs_0.0.3.tar.gz: build
	R CMD build build

build: data
	mkdir -p build
	cp -r R build
	cp -r data build
	cp -r man build
	cp DESCRIPTION build
	cp NAMESPACE build

data: data-raw/cisbp_orthologs.txt
	Rscript --vanilla data-raw/DATASET.R

data-raw/cisbp_orthologs.txt: data-raw/$(ENSEMBL)
	Rscript --vanilla data-raw/getOrthologs.R

$(ENSEMBL).zip:
	curl --insecure --output data-raw/$(ENSEMBL).zip $(URL)
	unzip -o data-raw/$(ENSEMBL).zip

clean:
	rm -f data-raw/$(ENSEMBL).zip*
	rm -rf build
