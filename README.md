# MicroC_RCMC_analysis
A processing pipeline for alignment, processing, and generation of contact maps for Micro-C and RCMC data.

## Setup
### Conda setup (can be skipped if you already have conda installed)
The pipeline uses conda environments to run all of the necessary tools for processing the data. If you do not have conda installed, install conda and mamba as follows (or by any other method appropriate for your system), answering `yes` to all asked questions:
```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```
Restart your terminal and conda will be activated.

### Environment setup
In order to run snakemake, a conda environment containing snakemake is needed. First, ensure the base environment is activated:
```
conda activate base
```
Next, create the required environment (conda can also be used here if mamba is not installed):
```
mamba env create --name snakemake-env -c conda-forge -c bioconda python=3.12 snakemake=9.12
```
Once the environment is made, activate it:
```
conda activate snakemake-env
```
### Preparing files to run
Create and move to a directory for your analysis:
```
mkdir /path/to/your/analysis/directory
cd /path/to/your/analysis/directory
```
Download the repository to the directory:
```
clone git https://github.com/ahansenlab/MicroC_RCMC_analysis
```
### Running demo analysis - Micro-C
The repository includes files for running a quick demo analysis on heavily downsampled data. The fastq files for this data contain reads mapped to the human chromosome 16 from the hg19 build. All of the files included are already configured to run on the demo data, but a genome file is needed. This can be downloaded with the following command, or you can use your own file (see below):
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
```
Run the pipeline, setting an appropriate number of cores for your system - this number serves as the upper limit for what your processes will use. If setting this lower than 8, also change the value set for the `threads` parameter in `config/config.yaml` using a text editor. This will first index the genome file that was downloaded - this may take ~2 hours. If you already have a `bwa-mem2` indexed genome to use, you can use this instead by modifying the `genome` parameter in `config/config.yaml` to have the full path to the location of the genome file and its associated 
```
snakemake --cores 8 --sdm conda
```
### Running demo analysis - RCMC
RCMC analysis can be run in much the same way as Micro-C analysis. To run in RCMC mode, comment out the blank `regions` entry in `config/config.yaml` and uncomment `regions: demo/capture_locus.bed`. This provides a `.bed` format file which indicates the position of a captured region (no capture was performed on the demo data, but the same procedure can be performed to run the pipeline this way). Once `config/config.yaml` has been modified, you can again run the demo analysis to produce a new `.mcool` file which only contains reads in the indicated region:
```
snakemake --cores 8 --sdm conda
```
