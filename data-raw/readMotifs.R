getSelex <- function(dname='Selex_seq_Cirobu_All'){
  # package for creating PWMs
  require(TFBSTools)

  # data locations inside dname
  selex <- paste0(dname,'/SELEX_project_TableS1_S3_Nitta_et_al_2019.csv')
  selexGene <- paste0(dname,'/SELEXgeneID.csv')
  pwms <- paste0(dname,'/pwms/')

  #read data
  selex <- read.csv(selex, stringsAsFactors = F)
  selexGene <- read.csv(selexGene,stringsAsFactors = F)
  # merge data.frames
  selex <- merge(selex,selexGene,by.x=1,by.y=1)
  # select desired field if a motif was found
  selex <- selex[selex$Selex.Motif.Found=="Yes", c(
    'Primary.gene.name', 'Best.transcript.model', "Family.Classification",
    "Best.Pfm.6mer","Best.Pfm.8mer"
  )]
  # get KHID
  selex$Best.transcript.model <- paste0(
	"KH2013:", sub('.v.*','',selex$Best.transcript.model)
  )
  # remove duplicates
  selex <- selex[!duplicated(selex[,-1]),]
  # rename motifs with the same Primary.gene.name
  selex[duplicated(selex$Primary.gene.name),"Primary.gene.name"] <- paste0(
    selex[duplicated(selex$Primary.gene.name),"Primary.gene.name"],'v2'
  )
  
  # read PFMs
  profileMatrix8mer <- apply(
    selex[,c("Best.Pfm.8mer","Best.Pfm.6mer")],
    1,
    # if Best.Pfm.8mer is missing, substitute Best.Pfm.6mer
    function(x) {
      if(x[1]==''){
        file <- x[2]
      }else{
        file <- x[1]
      }
      return(as.matrix(
        read.table(paste0(pwms,file),row.names=c("A","C","G","T"))
      ))
    }
  )
  
  # convert PFM to PWM
  profileMatrix8mer <- lapply(profileMatrix8mer,function(x) t(t(x)/apply(x,2,sum)))
  
  selex.pwm8mer <- mapply(
    PWMatrix,
    selex$Primary.gene.name,
    selex$Best.transcript.model,
    profileMatrix=profileMatrix8mer,
    tags=mapply(
      list,
      Family_Name=selex$Family.Classification,
      DBID.1=selex$Primary.gene.name,
      MoreArgs = list(Motif_Type="ANISEED",strand="*"),
      SIMPLIFY = F
    )
  )
  selex.pwm8mer <- do.call(PWMatrixList,selex.pwm8mer)
  return(selex.pwm8mer)
} 

getHomerMotifs <- function(homerdat='known.motifs'){
  motifs <- readLines(homerdat)
  motifs <- split(
    motifs,
    Reduce(function(x,y) if(y) c(x,x[length(x)]+1) else c(x,x[length(x)]),grepl('^>',motifs))
  )
  require(TFBSTools)
  motifs <- lapply(motifs,strsplit,'\t')
  homer.motifs <- do.call(PWMatrixList,lapply(
    motifs,
    function(x){
      profileMatrix <- sapply(x[-1],as.numeric)
      profileMatrix <- t(t(profileMatrix)/apply(profileMatrix,2,sum))
      row.names(profileMatrix) <- c("A","C","G","T")
      # ID <- sub('^>','',x[[1]][1])
      tags <- strsplit(x[[1]][2],'/')[[1]]
      tags <- c(strsplit(tags[1],'[()]')[[1]],tags[-1])
      DBID.1 <- tags[1]
      ID <- make.names(DBID.1,T)
      Family_Name <- if(
        is.na(tags[2])|grepl("\\?",tags[2])|grepl("\\.",tags[2])
      ){"NA"} else if(
        grepl("Bias",tags[2],ignore.case = T)|grepl("repeat",tags[2],ignore.case = T)
      ){"SeqBias"} else if(
        grepl("promoter",tags[2],ignore.case = T)
      ){"promoter"} else tags[2]
      return(PWMatrix(
        ID=ID,
        name=ID,
        profileMatrix = profileMatrix,
        tags = list(Family_Name=Family_Name,DBID.1=DBID.1,Motif_Type="HOMER")
      ))
    }
  ))
  names(homer.motifs) <- make.names(ID(homer.motifs),T)
  return(homer.motifs)
}

getCisbpMotifs <- function(cisbp="CISBP/Ciona_intestinalis_2017_11_01_3_41_pm/",promoters=c(
  "BREd","BREu","DCEI-DCEIII","DPE","DRE","E-box","Inr_fly","Inr_human",
  "ohler","Pause_button",'TATA-box',"TCT","XCPE","MTE"
)){
  require(TFBSTools)
  motifs <- lapply(
    list.files(
      paste0(cisbp,'/pwms_all_motifs/'),
      full.names = T),
    function(x){
      mat <- t(as.matrix(read.delim(x)[,c("A","C","G","T")]))
      if(dim(mat)[2]>0){
        sapply(1:ncol(mat),function(y){
          mat[which.max(mat[,y]),y] <<- mat[which.max(mat[,y]),y]-(sum(mat[,y])-1)
        })
        # mat <- mat/apply(mat,2,sum)
        return(list(ID=sub('.txt','',sub('.*\\/','',x)),profileMatrix = mat,strand='*'))
      }
    })
  motifs <- motifs[sapply(motifs,class)=="list"]
  
  motifid <- unlist(sapply(motifs,'[','ID'))
  
  motifdat <- read.delim(
    paste0(cisbp,'/TF_Information_all_motifs_plus.txt'),
    stringsAsFactors = F)
  addmotif <- motifid[!motifid%in%motifdat$Motif_ID]
  sapply(addmotif, function(x) {
    motifdat[x,c('Motif_ID',"DBID.1",'Family_Name')] <<- x
    motifdat[x,'Family_Name'] <<- sub('[0-9]+$','',motifdat[x,'Family_Name'])
    if(
      motifdat[x,'Family_Name']%in%promoters
    ) motifdat[x,'Family_Name'] <<- "promoter"
  })
  motifdat <- motifdat[!duplicated(motifdat$Motif_ID),]
  row.names(motifdat) <- motifdat$Motif_ID
  motifdat <- motifdat[motifid,]
  sapply(1:length(motifs),function(i) {
    motifs[[i]]$name <<- motifdat[i,'DBID.1']
    motifs[[i]]$tags <<- motifdat[i,]
  })
  motifs <- lapply(motifs,function(x) do.call(PWMatrix,x))
  motifs <- do.call(PWMatrixList,motifs)
  names(motifs) <- ID(motifs)
  return(motifs)
}

