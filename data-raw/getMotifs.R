# geneToMotif is a table with one row for each motif-to-gene association
# To test for enrichment of binding sites, we would ideally want to unambiguously
# pair each TF to a single motif.
# select least ambiguous motifs based on having fewest possible gene associations
# To do this, we begin by selecting all motifs unambiguously mapped to a single gene
# We add these association to the output and remove all TFs present in the output 
# from the input.
disambigMotifs <- function(x,y,res=NULL){
  motifs <- split(x$geneToMotif$KHID,x$geneToMotif$ID)
  tmp <- motifs[sapply(motifs,length)<=y]
  if(length(tmp)>0){
    res <- append(res,tmp)
    x$geneToMotif <- x$geneToMotif[!x$geneToMotif$KHID%in%unlist(res),]
    disambigMotifs(x,y,res)
  }else{
    x$motifs <- append(x$motifs,res)
    return(x)
  }
}

reduceMotifs <- function(motifs,geneToMotif){
  require(TFBSTools)
  require(motifmatchr)
  require(DBI)
  reduced <- Reduce(
    disambigMotifs,
    1:max(sapply(split(geneToMotif$KHID,geneToMotif$ID),length)),
    list(motifs=list(),geneToMotif=geneToMotif)
  )
  motifs <- motifs[names(reduced$motifs)]
  names(motifs) <- sapply(reduced$motifs,paste,collapse=';')
  return(motifs)
}

mergeGeneName <- function(x,y){
  expr <- sub('^([A-Za-z][A-Za-z0-9]*[A-Za-z])[0-9/]+$','\\1',x)
  sel <- sapply(paste0("^",expr,"[0-9/]+$"),grepl,y)
  if(any(sel)){
    nums <- as.numeric(unlist(strsplit(sub(expr,'',x),'/')))
    numsToAdd <- as.numeric(unlist(strsplit(sub(expr,'',y),'/')))
    nums <- unique(c(nums,numsToAdd))
    nums <- nums[order(nums)]
    y <- paste(as.character(nums),collapse='/')
    x[sel] <- paste0(expr,y)
  }else{
    x <- append(x,y)
  }
  return(x)
}

lowerGeneName <- function(x){
  if(
    grepl("[A-Z]{2}",x)&!(grepl("^KH\\.",x)|grepl("^SI:",x))
  ){
    tmp <- sub(".*?([A-Z])([A-Z]+)",'\\2',x)
    return(sub(tmp,tolower(tmp),x))
  } else {x}
}

khToName <- function(x,gene.names){
  x <- gene.names[x,"GeneName"]
  x <- gsub(
    '(;?)KH\\.[A-Z][0-9]+\\.[0-9]+_','\\1',gsub("KH2013:",'',x)
  )
  x <- sub(';$','',x)
  x <- sapply(
    strsplit(x,';'),
    function(y) paste(
      sapply(y,lowerGeneName),
      collapse=';')
  )
  return(x)
}

nameMotifs <- function(motifs,gene.names,khid.sub=F){
  require(TFBSTools)
  tf.family <- sapply(tags(motifs),'[[',"Family_Name")
  tf.family[grep("(Gata|Zn?F)",tf.family,T)] <- "Zinc finger"
  tf.family[tf.family=="NR"] <- "Nuclear receptor"
  tf.family[tf.family=="ETS"] <- "Ets"
  tf.family[tf.family=="MAD"] <- "SMAD"
  tf.family[tf.family=="MADS"] <- "MADS-box"
  tf.family[tf.family=="Paired"] <- "Homeodomain,POU"
  tf.family[tf.family=="Homeodomain,Paired box"] <- "Paired box"
  tf.family[tf.family=="Paired,Homeobox"] <- "Paired box"
  tf.family[tf.family=="EBF"] <- "HLH"
  tf.family[tf.family=="POU,Homeobox,HMG"] <- "HMG"
  tf.family[tf.family=="CTF,Forkhead"] <- "CTF"
  tf.family[tf.family=="ETS:IRF"] <- "IRF"
  tf.family[tf.family=="promoter"] <- "TBP"
  tf.family[tf.family=="Homeobox,bHLH"] <- "Homeodomain"
  tf.family[tf.family=="AP2"] <- "AP-2"
  tf.family[tf.family=="E2F/TDP"] <- "E2F"
  tf.family[tf.family=="Homeobox"] <- "Homeodomain"
  
  tf.tags <- mapply(function(x,y){
    x$Family_Name <- y
    return(x)
  },tags(motifs),tf.family)
  
  tf.name <- sapply(tags(motifs),'[[',"DBID.1")
  tf.kh <- names(motifs)
  
  tf.kh.gene <- strsplit(tf.kh,';')
  tf.kh.gene <- lapply(tf.kh.gene,khToName,gene.names)
  tf.kh.gene <- lapply(tf.kh.gene,function(x) x[!duplicated(x)])
  
  if(khid.sub){
    sel <- sapply(tf.kh.gene,function(x) all(grepl("^KH\\.[A-Z][0-9]+",x)))
    tf.kh.gene[sel] <- sapply(ID(motifs)[sel],list)
  }
  
  tf.kh.gene <- lapply(tf.kh.gene,sub,pattern='V\\$',replacement='')
  
  tf.kh.gene <- lapply(tf.kh.gene,Reduce,f=mergeGeneName)
  
  tf.kh.gene <- sapply(tf.kh.gene,paste,collapse=';')
  
  tf.family[tf.kh.gene=="Ctcf"] <- "Zinc finger"
  tf.family[tf.kh.gene=="Rbpj"] <- "CSL"
  
  motifs <- mapply(
    PWMatrix,
    names(motifs),
    tf.kh.gene,
    strand="*",
    tags=tf.tags,
    profileMatrix=Matrix(motifs)
  )
  names(motifs) <- make.unique(tf.kh.gene)
  
  motifs <- do.call(PWMatrixList,motifs)
  return(motifs)
}

mergeMotifs <- function(){

	geneid <- read.delim('gene_name.txt')
	names(geneid) <- c("KHID","GeneName")
	row.names(geneid) <- geneid$KHID

	# read PWMs
	selex.pwm <- getSelex()
	homer.pwm <- getHomerMotifs()
	cisbp.pwm <- getCisbpMotifs()

	# read motif-to-gene matches
	khToHomer <- read.delim('homer_orthologs.txt',stringsAsFactors=F,header=T)
	khToCisbp <- read.delim('cisbp_orthologs.txt',stringsAsFactors=F,header=T)

	# rename duplicated motifs in selex.pwm
	sel <- names(selex.pwm)%in%names(homer.pwm)
	names(selex.pwm)[sel] <- paste0(names(selex.pwm)[sel],"ANISEED")

	# append motifs
	comb.pwm <- Reduce(append,list(selex.pwm,homer.pwm,cisbp.pwm))

	# motif-to-KHID table
	khToMotif <- rbind(
	  cbind(ID=names(selex.pwm),KHID=name(selex.pwm)),
	  setNames(khToHomer[,c("ID","GeneID")],c("ID","KHID")),
	  setNames(khToCisbp[,c("Motif_ID","KHID")],c('ID',"KHID")),
	  stringsAsFactors=F
	)

	# remove motifs not in comb.pwm
	khToMotif <- khToMotif[khToMotif$ID%in%names(comb.pwm),]
	# remove duplicate rows
	khToMotif <- khToMotif[!duplicated(khToMotif),]

	# khToMotif is a table with one row for each motif-to-gene association
	# To test for enrichment of binding sites, we would ideally want to unambiguously
	# pair each TF to a single motif.
	# select least ambiguous motifs based on having fewest possible gene associations
	# To do this, we begin by selecting all motifs unambiguously mapped to a single gene
	# We add these association to the output and remove all TFs present in the output 
	# from the input.
	motifs <- reduceMotifs(comb.pwm,khToMotif)
	motifs <- nameMotifs(motifs,geneid)
	return(motifs)
}
