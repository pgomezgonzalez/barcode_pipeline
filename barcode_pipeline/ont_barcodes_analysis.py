#! /usr/bin/env python
"""Basecalling of ONT pod5 data with dorado and demultiplex"""


import argparse
import os
import glob
import pandas as pd
import sys
import subprocess as sp 
import shutil
import logging

script_path =  os.path.dirname(os.path.realpath(__file__))

def cli():

	argparser = argparse.ArgumentParser(description="Basecalling and demultiplexing of ONT pod5 data with dorado followed by internal barcodes analysis")

	argparser.add_argument("metadata",help="metadata with experiment_id,flow_cell_id,kit_id,barcode,alias,time_point,replicate")
	argparser.add_argument("sample_sheet",help="name of sample sheet for dorado basecalling")
	argparser.add_argument("kit_name", help="name of the kit used")
	argparser.add_argument("data", help="the data directory (folder with all pod5 files)")
	argparser.add_argument("output", help="name of output bam file")
	argparser.add_argument("ref_fasta", help="reference gene")
	argparser.add_argument("internal_barcodes", help="file with internal barcodes and variant_id")
	argparser.add_argument("region_bed", help="bed file with region where internal barcodes are")
	argparser.add_argument("prefix", help="prefix for all the results files")
	argparser.add_argument("--skip-basecalling",action="store_true", help="Skip basecalling if the flag is provided")
	argparser.add_argument("--only-basecalling",action="store_true", help="Only perform basecalling")
	argparser.add_argument("--allow-missmatch", action="store_true", help="on reads where no exact match for internal barcodes has been found, allow for a missmatch")
	argparser.add_argument("--missmatch",type=int, help="number of missmatches allowed",default=1)

	args = argparser.parse_args()

	#function to convert xlsx file to txt tab del
	def convert_xlsx_to_txt(file_name):
		#check if it's xlsx 
		if file_name.lower().endswith('.xlsx'):
			print(f"Converting {file_name} to tab-delimited .txt file...")
			df = pd.read_excel(file_name, engine='openpyxl')
			output_file = file_name.replace('xlsx','txt')
			df.to_csv(output_file, sep="\t", index=False)
			return output_file
		else:
			print(f'{file_name} already a txt file.')
			output_file=file_name
			return output_file 


	def create_output_directory(directory_name):
		#check if demux directory already exists
		if not os.path.exists(directory_name):
			os.makedirs(directory_name) #create demux directory 
			print(f"Directory '{directory_name}' has been created")
		else:
			print(f"Directory '{directory_name}' already exists")



	output_file = convert_xlsx_to_txt(args.metadata)

	sp.run(f'Rscript {script_path}/convert_barcode_names.R {output_file}', shell=True)

	if not args.skip_basecalling:
		##Make sample_sheet from metadata 
		print("*****Generating sample sheet*****")
		sp.run(f'Rscript {script_path}/make_dorado_samplesheet.R {output_file} {args.sample_sheet}', shell=True)

	############################################################################################################################################################
	############################################################################################################################################################
	################################################################--------BASECALLING------###################################################################

		#Run dorado basecaller 
		print("...starting basecalling...")
		sp.run(f'dorado basecaller --min-qscore 10 --kit-name {args.kit_name} --sample-sheet {args.sample_sheet}.csv -r sup {args.data} > {args.output}.bam', shell=True)

		if args.only_basecalling:
			sys.exit("...basecalling finished...EXITING...")
	else:
		print("...SKIPPING BASECALLING...")


	create_output_directory("demux")


	############################################################################################################################################################
	############################################################################################################################################################
	##############################################################--------DEMULTIPLEXING------##################################################################

	#Demultiplex into NP barcodes. It should create 1 bam file per barcode in demux directory
	print("...demultiplexing...")
	sp.run(f'dorado demux --output-dir demux --no-classify --emit-summary {args.output}.bam', shell=True)


	############################################################################################################################################################
	############################################################################################################################################################
	################################################################--------FILTERING------#####################################################################

	#Filter the reads based on quality score >=8
	create_output_directory("fastqs")

	if not args.skip_basecalling:

		sp.run(r"""ls ./demux/*barcode*.bam | sed 's/.*\(barcode[0-9]*\)\.bam/\1/' > list_bams""", shell=True)
		sp.run(r'ls ./demux/*barcode*.bam > list_demux_bams2',shell=True)
		sp.run(r'paste list_demux_bams2 list_bams > list_bams_final',shell=True)
		sp.run(r'''cat list_bams_final | parallel -j 1 --col-sep "\t" "samtools view -b -e '[qs]>=8' {1} | samtools fastq - | pigz -c > ./fastqs/{2}.fastq.gz"''',shell=True)
	
	else:
		print("removing not used barcodes")
		sp.run(r"""ls ./demux/*barcode*.bam | sed 's/.*\(barcode[0-9]*\)\.bam/\1/' > list_bams""", shell=True)
		sp.run(r'ls ./demux/*barcode*.bam > list_demux_bams2',shell=True)
		sp.run(r'paste list_demux_bams2 list_bams > list_all_bams_final',shell=True)
		sp.run(f'Rscript {script_path}/barcodes_used.R {output_file}',shell=True)
		sp.run(r'cat barcodes_used | parallel -j 1 "grep {} list_all_bams_final" > list_bams_final', shell=True)
		#change name of barcodes_used to list_bams to be the same as in running basecalling
		sp.run(r'mv barcodes_used list_bams', shell=True)
		print("...filtering and creating fastqs...")
		sp.run(r'''cat list_bams_final | parallel -j 1 --col-sep "\t" "samtools view -b -e '[qs]>=8' {1} | samtools fastq - | pigz -c > ./fastqs/{2}.fastq.gz"''',shell=True)


	############################################################################################################################################################
	############################################################################################################################################################
	#################################################################--------MAPPING------######################################################################

	#Map de reads 
	#we need a reference gene = fasta 
	create_output_directory("mapping")

	print("...aligning barcodes to reference...")
	sp.run(f'dorado aligner {args.ref_fasta} ./fastqs/ --output-dir ./mapping --emit-summary', shell=True)


	##Calculate the number of total reads, mapped reads and unmapped reads 
	sp.run(r'''cat list_bams | parallel -j 1 "samtools view ./mapping/{}.bam | wc -l >> total_reads"''', shell=True)
	sp.run(r'''cat list_bams | parallel -j 1 "samtools view -F4 ./mapping/{}.bam | wc -l >> mapped_reads"''', shell=True)
	sp.run(r'''cat list_bams | parallel -j 1 "samtools view -f4 ./mapping/{}.bam | wc -l >> unmapped_reads"''', shell=True)
	sp.run(r'paste list_bams total_reads mapped_reads unammaped_reads > table_number_reads', shell=True)


	############################################################################################################################################################
	############################################################################################################################################################
	################################################################--------ANALYSIS------######################################################################


	##One mapped, we can extract the number of reads per barcode for each internal barcode 
	##We need a file with internal barcodes and variant_id called internal_barcodes.txt

	shutil.copyfile(args.internal_barcodes,"internal_barcodes")

	##create bam files excluding unmapped reads 
	print("...creating mapped bam files...")
	sp.run(r'''cat list_bams | parallel -j 1 "samtools view -F4 ./mapping/{}.bam > ./mapping/{}.mapped.bam"''',shell=True)
	
	print("...counting reads...")

	sp.run(r'''while read line; do echo $line; cat internal_barcodes | parallel -j 1 --col-sep "\t" "samtools view ./mapping/$line.mapped.bam | grep {2} | wc -l >> barcode_read_count"; done < list_bams''', shell=True)
	sp.run(r'''while read line; do echo $line; cat internal_barcodes | parallel -j 1 --col-sep "\t" "samtools view ./mapping/$line.mapped.bam | grep {2} | awk '{print \$1}' >> read_ids_with_barcode"; done < list_bams''', shell=True)
	sp.run(r'''while read line; do cat internal_barcodes | parallel -j 1 --col-sep "\t" "echo $line '\t'{1} >> barcodes_variants"; done < list_bams''',shell=True)

	##CALCULATE COVERAGE???? 
	#calculate coverage
	#print("***calculating coverage***")
	#sp.run(f'''cat list_bams | parallel -j 1 "bedtools coverage -a {args.region_bed} -b ./mapping/{}.bam >> coverage.bed"''',shell=True)
	#sp.run(r'paste list_bams coverage.bed > coverage2.bed', shell=True)

	if args.allow_missmatch:
		#remove the reads that have picked up an internal barcode with exact match from bam file (into a different bam file called barcodeXX_rest.bam)
		sp.run(r'cat list_bams | parallel -j 1 "samtools view -h ./mapping/{}.bam | grep -vf read_ids_to_remove | samtools view -bS -o ./mapping/{}_rest.bam"', shell=True)
		#count reads allowing for a missmatch 
		sp.run(f'''while read line; do cat internal_barcodes | parallel -j 1 --col-sep "\t" "samtools view ./mapping/$line_rest.bam | agrep -n {args.missmatch} {2} | wc -l >> count_reads_missmatch"; done < list_bams''',shell=True)
		sp.run(r'''while read line; do cat internal_barcodes | parallel -j 1 --col-sep "\t" "echo $line_rest'\t'{1} >> barcodes_variants_missmatch"; done < list_bams''', shell=True)
		#calculate coverage
		#print("***calculating coverage***")
		#sp.run(f'cat list_bams | parallel -j 1 "bedtools coverage -a {args.region_bed} -b ./mapping/{}_rest.bam >> rest_coverage.bed"',shell=True)

		sp.run(r'paste barcodes_variants barcode_read_count count_reads_missmatch > count_reads_internal_barcodes',shell=True)
		#	sp.run(r'paste coverage2.bed rest_coverage.bed > barcodes_coverage.bed', shell=True)

	else:
		sp.run(r'paste barcodes_variants barcode_read_count > count_reads_internal_barcodes',shell=True)
		#sp.run(r'cp coverage2.bed barcodes_coverage.bed', shell=True)



	#Run Rscript to create plots and tables 
	print("...generating tables and plots...")
	sp.run(f'Rscript {script_path}/make_plots_tables.R count_reads_internal_barcodes.csv total_reads {output_file} barcodes_coverage.bed', shell=True)

	###Remove temporary files 
	#sp.run(r'rm list_bams list_demux_bams2 list_all_bams_final list_bams_final total_reads count_reads barcodes_variants count_reads_internal_barcodes.csv', shell=True)

	###Rename files with prefix of the run to differentiate 
	sp.run(f'mv lineplot_barcodes.png {args.prefix}.lineplot_barcodes.png', shell=True)
	sp.run(f'mv barplot_barcodes.png {args.prefix}.barplot_barcodes.png', shell=True)
	sp.run(f'mv results.xlsx {args.prefix}.results.xlsx',shell=True)
	sp.run(f'mv table_reads.txt {args.prefix}.table_reads.txt', shell=True)
	sp.run(f'mv table_proportions.txt {args.prefix}.table_proportions.txt', shell=True)
	sp.run(f'mv table_percentages.txt {args.prefix}.table_percentages.txt', shell=True)
	sp.run(f'mv internal_barcodes {args.prefix}.internal_barcodes', shell=True)
	sp.run(f'mv sample_sheet.csv {args.prefix}.sample_sheet.csv', shell=True)
	sp.run(f'mv summary_proportions.txt {args.prefix}.summary_proportions.txt',shell=True)	
	sp.run(f'mv coverage.bed {args.prefix}.coverage.bed',shell=True)
	sp.run(f'mv rest_coverage.bed {args.prefix}.rest_coverage.bed', shell=True)

	print("ALL DONE!")
