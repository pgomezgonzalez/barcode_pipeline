#! /usr/bin/env Rscript
#barcodes_used.R

###To remove those barcodes that haven't been used 

args = commandArgs(trailingOnly=TRUE)

#suppressWarnings(
    metadata <- read.table(args[1],sep="\t",header=T,stringsAsFactors=F)
    barcodes_used <- metadata$barcode 
    write.table(barcodes_used,file="barcodes_used",sep="\t",quote=F,row.names=F,col.names=F)
#)