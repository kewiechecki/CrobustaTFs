#!/bin/r
## convert CisBP orthologs from ENSEMBL to KH2013

if(sub('.*\\/','',getwd())!='data-raw'){
	 setwd('data-raw')
}

library(biomaRt)

ensembl <- read.delim('KH-ENS.blast',stringsAsFactors = F,header=F)

mart <- useMart('ensembl','cintestinalis_gene_ensembl')
# listAttributes(mart)
transcriptToGene <- select(mart,ensembl[,1],c('ensembl_gene_id','ensembl_transcript_id'),"ensembl_transcript_id")

ensembl.gene <- data.frame(ensembl_transcript_id=ensembl[,1],GeneID=sub('.v.*','',ensembl[,2]))
ensembl.gene <- merge(ensembl.gene,transcriptToGene,'ensembl_transcript_id')
ensembl.gene <- ensembl.gene[!duplicated(ensembl.gene[,-1]),-1]

cisbp2 <- read.delim('CisBP/Cintestinalis2.00/TF_Information_all_motifs_plus.txt',stringsAsFactors = F)
cisbp1.02 <- read.delim('CisBP/Cintestinalis1.02/TF_Information_all_motifs_plus.txt',stringsAsFactors = F)
cisbp.geneid <- merge(cisbp2,cisbp1.02,all.x=T,all.y=T)
cisbp.geneid <- merge(cisbp.geneid,ensembl.gene,by.x='DBID',by.y='ensembl_gene_id')
cisbp.geneid$GeneID <- sub("KH2012","KH2013",cisbp.geneid$GeneID)

gene.names <- read.delim("gene_name.txt", row.names=1, stringsAsFactors=F)
cisbp.geneid$GeneName <- gene.names[cisbp.geneid$GeneID,]
# dir.tab(cisbp.geneid,'cisbpGeneID')

cisbpDat <- cisbp.geneid[,c("Motif_ID","DBID.1","GeneID","GeneName","Family_Name","Motif_Type")]
cisbpDat <- cisbpDat[!duplicated(cisbpDat)&cisbpDat$Motif_ID!='.',]
write.table(cisbpDat, "cisbp_orthologs.txt", quote=F, sep='\t')
