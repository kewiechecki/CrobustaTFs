#!/bin/r
## write CrobustaMotifs to data/

if(sub('.*\\/','',getwd())!='data-raw'){
	 setwd('data-raw')
}
source("readMotifs.R")
source("getMotifs.R")
CrobustaMotifs <- mergeMotifs()
usethis::use_data(CrobustaMotifs, overwrite = TRUE)
devtools::document()
