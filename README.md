# barcode_pipeline

Package for analysis pipeline of Oxford Nanopore sequencing data with different barcodes. 
Counts the relative abundance of different mutants based on internal barcodes at different time points defined by the NP barcodes. 

## Install

Before installing the package, create an adequate conda environment with all dependencies needed. This can be done through the .yml file available. To do so, open the terminal and make a new directory for the yml file. Can be called ont_barcodes
The file with all the dependencies and packages needed is `ont_barcodes.yml`.

Copy the file into your profile in the server. In the terminal, go to the directory where you have the `ont_barcodes.yml` file (Downloads for example). 

```
cd ~/Downloads
scp ont_barcodes.yml username@10.18.0.25
```
It will ask you for your password, and then copy the file into your profile in the server. 

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

Now activate it to install the analyis package
```
conda activate ont_barcodes
```

And we can install the package by
``` 
pip install git+https://github.com/pgomezgonzalez/barcode_pipeline
```

## Usage
Once install we can run the pipeline. 
To do it you just need to type `barcodes-pipeline` and a message like this should appear
```
usage: barcodes-pipeline [-h] metadata sample_sheet kit_name data output ref_fasta internal_barcodes
barcodes-pipeline: error: the following arguments are required: metadata, sample_sheet, kit_name, data, output, ref_fasta, internal_barcodes
```

The nanopore sequencing data should be stored in a folder named by the run name and then another one with the date. 
There should there a folder called `pod5`, together with more files. 

Make an `analysis` folder 
