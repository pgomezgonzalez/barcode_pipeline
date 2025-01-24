# barcode_pipeline

Package for analysis pipeline of Oxford Nanopore sequencing data with different barcodes. 
Counts the relative abundance of different mutants based on internal barcodes at different time points defined by the NP barcodes. 

## Installation

Before installing the package, create an adequate conda environment with all dependencies needed. This can be done through the `.yml` file available. To do so, open the terminal and make a new directory for the yml file. Can be called `ont_barcodes`.
The file with all the dependencies and packages needed is `ont_barcodes.yml`.

Copy the file into your profile in the server. In the terminal, go to the directory where you have the `ont_barcodes.yml` file (Downloads for example). 

```
cd ~/Downloads
scp ont_barcodes.yml username@10.18.0.25
```
It will ask you for your password (same as your `username`), and then copy the file into your profile in the server. 

Now go to the server:
```
ssh -X username@10.18.0.25
```
And create the conda enviroment by typing 
```
mamba env create -f ont_barcodes.yml
```
It will take some time to download all the dependencies. Once done it will show a message like 
```
done
#
# To activate this environment, use
#
#     $ conda activate ont_barcodes
#
# To deactivate an active environment, use
#
#     $ conda deactivate
```

The conda environment will now be created. 

Activate it to install the analyis package:
```
conda activate ont_barcodes
```

And we can install the package by typing:
``` 
pip install git+https://github.com/pgomezgonzalez/barcode_pipeline
```

This should get it installed within the `ont_environment` and ready to use. 

## Usage

As it is installed in a conda environment, remember to activate always before using it the `ont_barcodes` environment. 
```
conda activate ont_barcodes
```

To run the program, you just need to type `barcodes-pipeline` and a message like this should appear:
```
usage: barcodes-pipeline [-h] metadata sample_sheet kit_name data output ref_fasta internal_barcodes
barcodes-pipeline: error: the following arguments are required: metadata, sample_sheet, kit_name, data, output, ref_fasta, internal_barcodes
```

If you want to see the help from the command just type `barcodes-pipeline -h` and you should see something like:
```
usage: barcodes-pipeline [-h] [--skip-basecalling] metadata sample_sheet kit_name data output ref_fasta internal_barcodes

Basecalling and demultiplexing of ONT pod5 data with dorado followed by internal barcodes analysis

positional arguments:
  metadata            metadata with experiment_id,flow_cell_id,kit_id,barcode,alias,time_point,replicate
  sample_sheet        name of sample sheet for dorado basecalling
  kit_name            name of the kit used
  data                the data directory (folder with all pod5 files)
  output              name of output bam file
  ref_fasta           reference gene
  internal_barcodes   file with internal barcodes and variant_id

options:
  -h, --help          show this help message and exit
  --skip-basecalling  Skip basecalling if the flag is provided
```

The arguments that need to be given to the command are:

- **metadata**: The file with metadata, that can be a xlsx or txt (tab-delimited) file with the following data in columns:
    - `experiment_id` = whatever name/ID you want to give to this experiment
    - `nanopore_run_id` = identifier for the nanopore run (i.e. AI_NP004)
    - `flow_cell_id` = flow cell identification used (i.e. FAZ01857)
    - `kit` = name of the kit used for library preparation (i.e. SQK-NBD114-96)
    - `sample_id` = name given to the samples run (each sample is 1 barcode)
    - `barcode` = list of barcodes used (i.e. barcode01)
    - `alias` = an alias name you want to give to each barcode (has to be different than barcode name, can be sample_id)
    - `time_point` = time point associated with each barcode
    - `replicate` = number of replicates of each time point

It should look like this:

```
experiment_id | nanopore_run_id | flow_cell_id |      kit      | sample_id |  barcode  |   alias   | time_point | replicate 
---------------------------------------------------------------------------------------------------------------------------
kk_GIA_1      |    KK_NP005     |   FAZ01857   | SQK-NBD114-96 | PvDBPa    | barcode01 | time1rep1 |      1     |     1     
---------------------------------------------------------------------------------------------------------------------------
kk_GIA_1      |    KK_NP005     |   FAZ01857   | SQK-NBD114-96 | PvDBPa    | barcode02 | time1rep2 |      1     |     2
----------------------------------------------------------------------------------------------------------------------------
kk_GIA_1      |    KK_NP005     |   FAZ01857   | SQK-NBD114-96 | PvDBPa    | barcode04 | time2rep1 |      2     |     1
...
...
```
<mark><p style ="text-algin: center;">If you make it in excel, save it as .xlsx or tab-delimited (.txt) file!</mark></p>

- **sample_sheet**: this is just the name of the sample sheet that the basecaller will take to only allocate reads to the barcodes used. Recommend to just leave it as `sample_sheet`.

- **kit_name**: name of the kit used. For example, SQK-NBD114-96.

- **data**: this is the folder containing all the pod5 files, the nanopore data for basecalling (usually will be called `pod5`). 

- **output**: this is the name you want to give to the file containing the reads after basecalling (i.e. kk_GIA_20240125)

- **ref_fasta**: this is the reference gene/amplicon in fasta format to which the sequencing reads will be aligned. 

- **internal_barcodes**: a file containing the sequence of internal barcodes (second column) and variant they represent (first column). This should be a tab delimited file (.txt) with 2 columns and no header:
```
---------------------------------
DBP_P	|   TCGAGAAAGAATAGCAGTAA
---------------------------------
DBP_O	|  AGTCGCAAGAATTCGTCTAAT
---------------------------------
DBP_AH	|  TCGAGAAAGAACAGCTCAAAT
...
...
```


The nanopore sequencing data should be stored in a folder named by the run name and then another one with the date. 
There should be there a folder called `pod5`, together with more files. 

Enter the directory of the run where you can see the `pod5` folder and other files.

To have all the files/data needed together, make there an `analysis` folder, and transfer here the required files: 
 - metadata file
 - ref_fasta file
 - internal_bacodes file

 ```
 mkdir analysis
 ```
To copy the files we can use `scp` 
In your local terminal, go to the directory where you have the files.
Once there, you can copy them to your profile in the server with:
```
scp metadata.txt internal_barcodes.txt amplicon.fa username@10.18.0.25:<run_folder>/<run_date_folder>/analysis/
```

The pipeline runs the whole analysis from basecalling of the pod5 files and outputs a results spreadsheet and 2 plots (line plot and barplot). 
But if the basecalling process has already been done, we can skip basecalling, as it takes few hours, and go directly to demultiplex the reads into barcodes.

- To run the whole pipeline (basecalling + demultiplex + analysis)
```
barcodes-pipeline <metadata> <sample_sheet> <kist_name> <data> <output> <ref_fasta> <internal_barcodes>
```
For example:
```
barcodes-pipeline metadata.txt sample_sheet SQK-NBD114-96 ../pod5 AI_NP004ab_merged_basecalled_bam_dorado_sup.allfiles  dbp_template.fa internal_barcodes.txt
```

Keep in mind that the `pod5` folder, if we have created an analysis one and we are running the command from the analysis one, will be in `../pod5` (we need to give the directory path as well). 


- To skip the basecalling process we need to use the flag `--skip-basecalling` and use the same arguments
```
barcodes-pipeline --skip-basecalling <metadata> <sample_sheet> <kist_name> <data> <output> <ref_fasta> <internal_barcodes>
```
When skipping basecalling, just make sure that the `<output>` argument is set to the bam file with all the reads (output of the previos basecalling process), but without the extension (`.bam`)!

### outputs

The outputs of the pipeline will be:
- all reads bam file (named as your `<output>` argument `.bam`)
- a `demux` folder containing the demultiplexed reads into barcodes 
- a `sample_sheet.csv` with the infor for dorado basecalling process 
- a `fastqs` folder with fastq files filtered by quality score 
- a `mapping` folder with the mapped bam files 
- `results.xlsx` file with 5 tabs: table with reads per barcode/internal barcode; same with proportions; same with percentages; a long version of the same percentages table used for R plots; and a summary table with means and standard deviations 
- a `line_plot_barcodes.png`
- a `barplot_barcodes.png`

You can transfer your results to your computer from your computer terminal. Go to the directory 


## Upgrading package

If there has been any upgrades in the package, you can reinstall the new version with 
```
conda activate ont_barcodes     #Make sure you are in the right environment where it is installed
pip install --upgrade --force-reinstall git+https://github.com/pgomezgonzalez/barcode_pipeline
```
