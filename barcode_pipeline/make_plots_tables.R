#! /usr/bin/env Rscript
#make_plots_tables.R

args = commandArgs(trailingOnly=TRUE)

pkgs <- c("ggplot2","reshape2","plyr","dplyr","openxlsx","viridis")

installed_packages <- pkgs %in% rownames(installed.packages())
if(any(installed_packages==FALSE)){
  install.packages(pkgs[!installed_packages],repos="https://cran.uk.r-project.org/")
}

library(ggplot2, quietly=TRUE)
library(reshape2, quietly=TRUE)
library(plyr, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(openxlsx, quietly=TRUE)
library(viridis, quietly=TRUE)

#suppressWarnings(
  #read data of count reads (should be called count_reads_internal_barcodes and be tab delimited file with 3 or 4 columns (barcode,variant_id and read_count/read _count_missmatch))

  data <- read.table(args[1],sep="\t",header=F,stringsAsFactors=FALSE)

  #read table with barcode, total number of reads, mapped reads, and unmapped reads (table_number_reads)
  total <- read.table(args[2],sep="\t",header=F,stringsAsFactors=FALSE)

  #read metadata file 
  metadata <- read.table(args[3],sep="\t",header=TRUE,stringsAsFactors=FALSE)

  #read coverage file (should be 8 columns, but the 5 ones are the important ones: barcode ref start end mean_coverage) - it only covers the barcode region
  #coverage <- read.table(args[4],sep="\t",header=FALSE,stringsAsFactors=FALSE)

  ##Make a table of proportion of reads per barcode/internal barcode 

  NP_barcodes <- metadata$barcode
  internal_barcodes <- unique(data$V2)

  df <- data.frame(matrix(ncol=length(internal_barcodes),nrow=length(NP_barcodes)))
  colnames(df) <- internal_barcodes
  df$NP_barcode <- NP_barcodes
  df$sample_id <-metadata$sample_id
  df$concentration <- metadata$concentration

  for(i in 1:nrow(df)){
	  subset <- data[which(data$V1==df$NP_barcode[i]),]
	  if(length(data)==4){
      reads <- subset$V3 + subset$V4
    }else{
      reads <- subset$V3
    }
	  df[i,1:length(internal_barcodes)] <- reads
  }

  df$total_reads_barcodes <- rowSums(df[1:length(internal_barcodes)])
  for(i in 1:nrow(df)){
    idx <- which(total$V1==df$NP_barcode[i])
    df$total_reads[i] <- total$V2[idx]
    df$mapped_reads[i] <- total$V3[idx]
    df$unmapped_reads[i] <- total$V4[idx]
  }

  df$time_point <- "" 
  df$replicate <- "" 

  for(i in 1:nrow(df)){
	  idx <- which(metadata$barcode==df$NP_barcode[i])
	  df$time_point[i] <- metadata[idx,"time_point"]
	  df$replicate[i] <- metadata[idx,"replicate"]
  }



  #make a normalised column per each of them 
  mat1 <- as.matrix(df[1:nrow(df),1:length(internal_barcodes)])
  mat2 <- mat1 

  mat2 <- as.matrix(mat1 / df$total_reads_barcodes)
  mat3 <- as.matrix(mat2*100)
  mat3 <- as.data.frame(mat3)
  mat3$total_percentage <- rowSums(mat3)

  write.table(df,file="table_reads.txt",sep="\t",quote=F,row.names=F)

  df2 <- as.data.frame(mat2)
  df2$NP_barcode <- NP_barcodes
  df2$sample_id <- df$sample_id
  df2$concentration <- df$concentration
  df2$total_reads_barcodes <- df$total_reads_barcodes
  df2$total_reads <- df$total_reads
  df2$time_point <- "" 
  df2$replicate <- "" 

  for(i in 1:nrow(df2)){
	  idx <- which(metadata$barcode==df2$NP_barcode[i])
	  df2$time_point[i] <- metadata[idx,"time_point"]
	  df2$replicate[i] <- metadata[idx,"replicate"]

  }

  write.table(df2,file="table_proportions.txt",sep="\t",quote=F,row.names=F)

  df3 <- as.data.frame(mat3)
  df3$NP_barcode <- NP_barcodes
  df3$sample_id <- df$sample_id
  df3$concentration <- df$concentration
  df3$total_reads_barcodes <- df$total_reads_barcodes
  df3$total_reads <- df$total_reads
  df3$time_point <- "" 
  df3$replicate <- "" 

  for(i in 1:nrow(df3)){
	  idx <- which(metadata$barcode==df3$NP_barcode[i])
	  df3$time_point[i] <- metadata[idx,"time_point"]
	  df3$replicate[i] <- metadata[idx,"replicate"]

  }

  write.table(df3,file="table_percentages.txt",sep="\t",quote=F,row.names=F)

  #generate averages by timepoint and by antibody used using the replicates

  df3$time_point <- "" 
  df3$replicate <- "" 

  for(i in 1:nrow(df3)){
	  idx <- which(metadata$barcode==df$NP_barcode[i])
	  df3$time_point[i] <- metadata[idx,"time_point"]
	  df3$replicate[i] <- metadata[idx,"replicate"]

  }

  ##Make line_plot
  df3_nototal <- subset(df3,select=-c(total_reads_barcodes,total_reads,total_percentage))

  ##Remove total_reads and total percentage 
  df3_nototal$time_point <- as.numeric(df3_nototal$time_point)

  df3_melt <- melt(df3_nototal,id=c("NP_barcode","time_point","replicate","sample_id","concentration"))



  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

      # New version of length which can handle NA's: if na.rm==T, don't count them
      length2 <- function (x, na.rm=FALSE) {
          if (na.rm) sum(!is.na(x))
          else       length(x)
      }

      # This does the summary. For each group's data frame, return a vector with
      # N, mean, and sd
      datac <- ddply(data, groupvars, .drop=.drop,
        .fun = function(xx, col) {
          c(N    = length2(xx[[col]], na.rm=na.rm),
            mean = mean   (xx[[col]], na.rm=na.rm),
            sd   = sd     (xx[[col]], na.rm=na.rm)
          )
        },
        measurevar
      )

      # Rename the "mean" column    
      # datac <- rename(datac, c("mean" = measurevar))

      datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

      # Confidence interval multiplier for standard error
      # Calculate t-statistic for confidence interval: 
      # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
      ciMult <- qt(conf.interval/2 + .5, datac$N-1)
      datac$ci <- datac$se * ciMult

      return(datac)
  }

  summary_df3 <- summarySE(df3_melt,measurevar="value",groupvars=c("time_point","sample_id","concentration","variable"))

  timepoints <- unique(df3_nototal$time_point)
  a <- min(timepoints)
  b <- max(timepoints)

  #Make barplots for proportion of reads 

  #calculate the number of rows and columns we want based on the number of variables (variables = antibodies/control)
  #We can make it be between 5 and 8 

  ab <- length(unique(summary_df3$sample_id)) ## number of antibodies and control used 

  if(ab %% 6 == 0){
    plot_rows <- ab/6
    plot_cols <- 6
  }else{
    plot_rows <- ceiling(ab/6)
    plot_cols <- 6
  }

  ##Make bar plot 

  ggplot(summary_df3,aes(x=concentration,y=mean,fill=variable)) + geom_bar(position="stack",stat="identity") + 
  theme_classic() + scale_fill_viridis(discrete=TRUE) + facet_wrap(.~sample_id,ncol=plot_cols,nrow=plot_rows) +
  labs(x="Concentration of antibody",y="Mean proportion of reads",fill="Variant")

  plot_height <- 3.5*plot_rows

  ggsave("barplot_barcodes.png",width=12,height=plot_height)

  summary_df3.2 <- subset(summary_df3,select=-c(ci))
  write.table(summary_df3.2,file="summary_proportions.txt",sep="\t",quote=F,row.names=F)

  list_results <- list("reads" = df, "proportions" = df2, "percentages" = df3, "percentages_long" = df3_melt, "summary" = summary_df3.2)
  write.xlsx(list_results,file="results.xlsx")


  ##Make a line plot 
  #The time point 0 is only for the controls with concentration 0, then time point 1 is for controls 
  # plus everything else with different concentrations 
  ##curate the table for it
  #add the time_point 0 - control for all sample_id 
  #summary_df3.3 <- summary_df3 
  #zero_tp <- summary_df3[which(summary_df3$time_point==0),]
  #sample_ids <- unique(summary_df3$sample_id)

  #list_zero_tps <- list()
  #for(i in 1:length(sample_ids)){
  #  zero_tp_i <- zero_tp
  #  zero_tp_i$sample_id <- sample_ids[i]
  #  list_zero_tps[[i]] <- zero_tp_i  
  #}

  #toremove <- which(summary_df3$time_point==0)
  #summary_df3.3 <- summary_df3.3[-toremove,]
  #list_zero_tps_df <- do.call("rbind",list_zero_tps)
  #summary_df3.3 <- rbind(summary_df3.3,list_zero_tps_df)

  #ggplot(summary_df3.3,aes(x=time_point,y=mean,col=variable)) + geom_point(size=2) + geom_line() + 
  #geom_errorbar(aes(ymin=mean-sd,ymax=mean-sd),colour="black",width=.2) + theme_classic() + scale_color_viridis(discrete=TRUE) +
  #labs(x="Time point",y="Proportion of reads",colour="Variant") + scale_x_continuous(breaks=seq(a,b,by=1)) +
  #facet_wrap(.~interaction(sample_id,concentration))


  #ggplot(summary_df3,aes(x=time_point,y=mean,col=variable)) + geom_point(size=2) + geom_line() + 
  #geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),colour="black",width=.2) + theme_classic() +
  #labs(x="Time point", y="Proportion of reads",colour="Internal barcode",title="Mean percentage of reads in each pool mapping to each internal barcode") + 
  #scale_x_continuous(breaks=seq(a,b,by=1))

  #ggsave("lineplot_barcodes.png")


  ##Make barplot
  #ggplot(summary_df3,aes(x=time_point,y=mean,fill=variable)) + geom_bar(position="stack",stat="identity") +
  #theme_classic() + labs(x="Time point",y="Proportion of reads",fill="Internal barcode") 

  #ggsave("barplot_barcodes.png")

#)


