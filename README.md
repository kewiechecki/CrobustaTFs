# CrobustaTFs
An R package containing candidate *Ciona robusta* motifs developed for use with [ATAC-seq data](https://doi.org/10.7554/eLife.49921). It integrates motifs identified using [SELEX-seq](https://doi.org/10.1007/978-1-4939-9624-7_23) with known orthologs from [HOMER](http://homer.ucsd.edu/homer/index.html) and inferred orthologs from [Cis-BP](doi: 10.1016/j.cell.2014.08.009). The install script attempts to retrieve all source files and runs a Gale-Shapley-like algorithm to map candidate motifs to genes with minimal ambiguity

## Rationale
Motifs in *C. robusta* are much less well characterized than *H. sapiens* or *M. musculus*. The most common techniques to fill in the gaps use either orthologs from other chordates or to infer orthologs from protein DNA binding domains. 
These methods naturally raise the question of which motif is "best." The answer often comes down to experimenter judgment, and different sources may yield more useful motifs for different TFs. 
The problem is particularly pronounced with inferred orthologs, as there are often many-to-many relationships between motifs and TFs.
It is often desirable to be able to unambiguously assign a motif to a singe TF. This is done by iteratively removing associations of ambiguous motifs to TFs that are already unambiguously associated to a motif.
This package aims to simplify mapping motifs to TFs in *C. robusta* by unifying possible sources into a single database.

## Installation
Depends: 
    R (>= 2.10),
    TFBSTools

Imports:
    DBI,
    motifmatchr
```bash
git clone https://github.com/kewiechecki/CrobustaTFs
cd CrobustaMotifs
make
```

## Usage
```{R}
library(CrobustaMotifs)

# accessing tags(CrobustaMotifs) requires loading the TFBSTools library separately
library(TFBSTools)

# function to compress motif metadata into a list of data.frames
tfTags <- function(motifs,tag){
	res <- setNames(do.call(data.frame,
		lapply(tag,
		       function(x) sapply(tags(motifs),'[[',x))),
		tag)
	return(res)
}

# extract fields from tags(CrobustaMotifs)
meta <- tfTags(CrobustaMotifs, 
	       names(tags(CrobustaMotifs[[1]])))

# split metadata by TF
meta <- split(meta,meta$KYID)

# select lowest mean shannon entropy motif for each TF
sel <- sapply(meta,function(x) row.names(x)[which.min(x$meanEntropy)])
motifsub <- CrobustaMotifs[sel]

```
