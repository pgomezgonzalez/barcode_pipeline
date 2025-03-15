#! /usr/bin/env Rscript
#coverage_plot.R

args = commandArgs(trailingOnly=TRUE)

library(dplyr, quietly=TRUE)
library(ggplot2, quietly=TRUE)

suppressWarnings(

    cov_files <- list.files(path="./coverage/",pattern="*.cov.bed",full.name=TRUE)

    list_cov <- lapply(cov_files,read.table)
    names(list_cov) <- sub(".cov.bed","",cov_files)
    names(list_cov) <- sub("./coverage//","",names(list_cov))

    metadata <- read.table(args[1],sep="\t",header=T,stringsAsFactors=F)

    for(i in 1:length(list_cov)){
        id <- which(metadata$barcode==names(list_cov)[[i]])
        subset <- list_cov[[i]]
        subset$barcode <- names(list_cov)[[i]]
        subset$sample_id <- metadata$sample_id[id]
        subset$concentration <- metadata$concentration[id]
        subset$replicate <- paste("Replicate ",metadata$replicate[id],sep="")
        list_cov[[i]] <- subset
    }

    ##Combine all data frames from list
    combined_list <- do.call(rbind, list_cov)

    ##Make data frame for region of barcode 
    data_barcode_region <- data.frame(matrix(nrow=length(unique(metadata$sample_id)),ncol=3))
    colnames(data_barcode_region) <- c("sample_id","start","end")
    data_barcode_region$sample_id <- unique(metadata$sample_id)
    region_data <- read.table(args[2],sep="\t",header=FALSE,stringsAsFactors = FALSE)
    data_barcode_region$start <- as.numeric(region_data$V2)
    data_barcode_region$end <- as.numeric(region_data$V3)

    max_y_vals <- combined_list %>% group_by(sample_id) %>% summarise(max_y = max(V3))
    merged <- left_join(max_y_vals,data_barcode_region,by="sample_id")

    combined_list$concentration <- as.factor(combined_list$concentration)

    ggplot(combined_list, aes(x=V2,y=V3,color=concentration)) + geom_line(linewidth=0.2) + theme_classic() + facet_grid(sample_id~replicate, scales="free") +
    geom_rect(data=merged,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=0,ymax=max_y),fill="orange",alpha=0.2) +
    xlab("gene position (bp)") + ylab("coverage")

    ggsave("coverage_plots.png")

)