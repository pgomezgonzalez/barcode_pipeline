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
It will ask you 
The conda environment will now be created 
pip install balblabl
```
