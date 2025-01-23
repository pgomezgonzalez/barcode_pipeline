#!/usr/bin/Rscript
#make_dorado_samplesheet.R

##Transform metadata file into sample_sheet for dorado

args = commandArgs(trailingOnly=TRUE)

#args[1] = metadata file (should be in .txt tab delimited)
#columns to have
###experiment_id
###flow_cell_id
###kit
###sample_id
###barcode
###alias
###time_point
###replicate 

pkgs <- "dplyr"
installed_packages <- pkgs %in% rownames(installed.packages())
if(any(installed_packages==FALSE)){
  install.packages(packages[!installed_packages])
}

library(dplyr)
metadata <- read.table(args[1],sep="\t",header=T,stringsAsFactors=FALSE)

#For DORADO sample sheet we need a csv file comma-separated with 
##experiment_id
##flow_cell_id
##kit
##barcode
##alias

#Find those columns in metadata and create new data frame 

columns_sample_sheet <- c("experiment_id","flow_cell_id","kit","barcode","alias")

#substract columns that match the columns_sample_sheet
df <- metadata %>% dplyr::select(matches(columns_sample_sheet))

write.table(df,file=paste(args[2],".csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
