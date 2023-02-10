library(CrobustaTFs)
library(tfenrichr)
library(Biostrings)
library(BSgenome.Crobusta.HT.KY)
library(TFBSTools)
library(purrr)

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

enhancer <- readDNAStringSet("depdc.fa")

motifs <- matchMotifs(motifsub,enhancer,
		      Crobusta,"genome","positions")
moitifs <- Filter(compose(partial(`!=`,0),length,
			  partial(`[[`,i=1)),motifs)
motifs <- Filter(compose(partial(`!=`,0),
	       length,partial(`[[`,i=1)),motifs)
motifs <- lapply(names(motifs),
		 function(x){
			 y <- motifs[[x]][[1]]
			 names(y) <- rep(x,length(y))
			 y
		 })
motifs <- unlist(do.call(IRangesList,motifs))
meta <- mcols(motifs)
mcols(motifs) <- NULL
motifs <- GRanges("Depdc.enhancer",motifs,strand=meta$strand,score=meta$score)
export(motifs,'depdc_motifs.bed')
