#! /usr/bin/env Rscript
## convert_barcode_names.R

args = commandArgs(trailingOnly=TRUE)

#suppressWarnings(
  metadata <- read.table(args[1],sep="\t",header=T,stringsAsFactors=F)

  ##add "barcode" or "barcode0" to the barcode column if not given 

  barcodes <- metadata$barcode
  new_barcodes <- c()

  for(i in 1:length(barcodes)){
    if(length(grep("barcode*",barcodes[i])) == 0){
     if(nchar(barcodes[i])==1){
        new_barcodes[i] <- paste("barcode0",barcodes[i],sep="") 
      }else{
        new_barcodes[i] <- paste("barcode",barcodes[i],sep="")
     }
   }else{
      new_barcodes[i] <- barcodes[i]
   }
  }

  metadata$barcode <- new_barcodes

  write.table(metadata,file=args[1],sep="\t",quote=F,row.names=F)
#)