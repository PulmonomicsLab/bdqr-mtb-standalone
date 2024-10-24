# BDQR-MTB-standalone

Prediction of Bedaquiline (BDQ) resistance in MTB clinical isolates using XAI

* [Workflows](#workflows)
* [Pre-requisites](#pre-requisites)
* [Installation](#installation)
* [Usage](#usage)
  * [1. Prediction of BDQ resistance class of one or many MTB paired-end FASTQ files](#1-prediction-of-bdq-resistance-class-of-one-or-many-mtb-paired-end-fastq-files)
  * [2. Prediction of BDQ resistance class of a single MTB VCF file](#2-prediction-of-bdq-resistance-class-of-a-single-mtb-vcf-file)
  * [3. Prediction of BDQ resistance class of multiple MTB VCF files](#3-prediction-of-bdq-resistance-class-of-multiple-mtb-vcf-files)
  * [4. Prediction of BDQ resistance class of merged MTB VCF file](#4-prediction-of-bdq-resistance-class-of-merged-mtb-vcf-file)
* [Demo runs](#demo-runs)
  * [1. Prediction from FASTQ (Single isolate)](#1-prediction-from-fastq-single-isolate)
  * [2. Prediction from FASTQ (Multiple isolates)](#2-prediction-from-fastq-multiple-isolates)
  * [3. Prediction from single VCF](#3-prediction-from-single-vcf)
  * [4. Prediction from multiple VCFs](#4-prediction-from-multiple-vcfs)
  * [5. Prediction from merged VCFs](#5-prediction-from-merged-vcfs)
* [Team](#team)
* [Disclaimer](#disclaimer)

## Workflows

There are 4 different scripts for prediction:
* ***bdqr-WGS-predict.sh:*** Prediction of BDQ resistance of single or multiple MTB WGS data (paired-end *fastq* files).
* ***bdqr-VCF-predict.sh:*** Prediction of BDQ resistance of single MTB VCF file.
* ***bdqr-multi-VCF-predict.sh:*** Prediction of BDQ resistance of multiple MTB VCF files.
* ***bdqr-merge-VCF-predict.sh:*** Prediction of BDQ resistance of MTB merged VCF file.

## Pre-requisites
* trim-galore (version `0.6.7`) - quality check and trimming of read sequences
* bwa (version `0.7.17-r1188`) - reference based alignment
* samtools (version `1.13`)- processing the BAM files
* freebayes (version `v1.3.6`) - variant calling
* libvcflib-tools (version `1.0.7`) - processing the VCF files
* libvcflib-dev (version `1.0.7`)	- processing the VCF files
* bgzip (version `1.13+ds`) - zipping files
* R (version `4.2.2 Patched (2022-11-10 r83330)`) - data operation and machine learning model usage


## Installation
### Step 1: Install dependent packages/tools
*(For ubuntu)*

    sudo apt-get install trim-galore bwa samtools freebayes libvcflib-tools libvcflib-dev bgzip

*(For other distributions)*

The installation steps for the different packages/tools are given in the following links:

* trim-galore - https://github.com/FelixKrueger/TrimGalore
* bwa - https://github.com/lh3/bwa
* samtools, bcftools, bgzip(htstools) - http://www.htslib.org/download/
* freebayes - https://github.com/freebayes/freebayes
* vcflib - https://github.com/vcflib/vcflib

R should be installed in the user system/PC. R installation steps are given in https://cran.r-project.org/.
The following R packages are required for the prediction of BDQ resistance class:

    install.packages("caret")
    install.packages("reshape2")
    install.packages("ggplot2")
    install.packages("DALEX")
    install.packages("RSNNS")

----
### Step 2: Install bdqr-mtb-standalone
  #### I. Download the software from GitHub repository
   Create a clone of the repository:

      git clone https://github.com/PulmonomicsLab/bdqr-mtb-standalone

   **Note:** Creating a clone of the repository requires `git` to be installed.

   The `git` can be installed using

      sudo apt-get install git
  **OR**

   Download using `wget`:

      wget https://github.com/PulmonomicsLab/bdqr-mtb-standalone/archive/refs/heads/main.zip
      unzip main.zip
  **Note:** wget can be installed using

      sudo apt-get install wget

  <br>

  #### II. Make the shell scripts executable

    chmod +x INSTALLATION_DIR/bdqr-mtb-standalone/config.sh
    chmod +x INSTALLATION_DIR/bdqr-mtb-standalone/bdqr-WGS-predict.sh
    chmod +x INSTALLATION_DIR/bdqr-mtb-standalone/bdqr-VCF-predict.sh
    chmod +x INSTALLATION_DIR/bdqr-mtb-standalone/bdqr-multi-VCF-predict.sh
    chmod +x INSTALLATION_DIR/bdqr-mtb-standalone/bdqr-merge-VCF-predict.sh

  `INSTALLATION_DIR` = Directory where bdqr-mtb-standalone is installed

  <br>

  #### III. update the paths in config.sh (optional)

  The `config.sh` looks like

	freebayes_path=/usr/bin/freebayes
	samtools_path=/usr/bin/samtools
	bwa_path=/usr/bin/bwa
	trim_galore_path=/usr/bin/trim_galore
	vcflib_path=/usr/bin/vcflib
	bgzip_path=/usr/bin/bgzip
	bcftools_path=/usr/bin/bcftools
	trim_galore_cores=4
	bwa_mem_cores=4
	samtools_cores=4

Note: It shows the default paths of the executables files for `freebayes`, `samtools`, `bwa`, `trim galore!`, `vcflib`, `bgzip` and `bcftools`. The users need to update the paths of the executables, in case these tools were installed in ways other than the `apt-get install` command.

## Usage

Initially change the directory to the directory where bdqr-mtb-standalone is installed.

    cd INSTALLATION_DIR/bdqr-mtb-standalone

Different operations can be performed by calling the appropriate scripts with two command-line arguments: `INPUT_DIR` and `OUTPUT_DIR`.

`INPUT_DIR` = the path (absolute or relative) of the folder containing the input files.

`OUTPUT_DIR` = the path (absolute or relative) of the folder in which bdqr-mtb-standalone will store the outputs.

The executable script, and contents of `INPUT_DIR` and `OUTPUT_DIR` depends on the choice of operations. The different operations are explained below.

### 1. Prediction of BDQ resistance class of one or many MTB paired-end FASTQ files
**(*.fastq* to BDQ resistance class (*Susceptible* or *Resistance*))**

    ./bdqr-WGS-predict.sh INPUT_DIR OUTPUT_DIR

`INPUT_DIR` must contain paired end FASTQ (*<sample_id>*_1.fastq.gz & *<sample_id>*_2.fastq.gz) of 1 or more isolates.

`OUTPUT_DIR` will contain
* one folder for each ISOLATE ID. Each folder will contain
  * the VCF file (*<sample_id>*.vcf)
  * the intermediate BAM files (*<sample_id>*.bam, *<sample_id>*_sorted.bam)
* the MERGED.vcf file (only in case of multiple isolates)
* the intermediate TSV file (*<sample_id>*.tsv/MERGED.tsv)
* the prediction result performed with the full model (prediction.tsv)
* the SHAP result for each isolate with the 50-feature model (shap_result_50_features_*<sample_id>*.tsv)
* the SHAP result plot for each isolate with the 50-feature model (shap_plot_50_features_*<sample_id>*.svg)

**Note:** The merged.vcf file will not be present if there was only one isolate in the INPUT_DIR. <br/>
The prediction and the SHAP results output will also be displayed on the terminal.

----

### 2. Prediction of BDQ resistance class of a single MTB VCF file
**(.vcf to BDQ resistance class (*Susceptible* or *Resistance*))**

    ./bdqr-VCF-predict.sh INPUT_DIR OUTPUT_DIR

`INPUT_DIR` must contain VCF file (*<sample_id>*.vcf)

`OUTPUT_DIR` will contain
* intermediate .tsv file (*<sample_id>*.tsv)
* the prediction result performed with the full model (prediction.tsv)
* the SHAP result for the isolate with the 50-feature model (shap_result_50_features_*<sample_id>*.tsv)
* the SHAP result plot for the isolate with the 50-feature model (shap_plot_50_features_*<sample_id>*.svg)

The prediction and the SHAP results output will also be displayed on the terminal.

----

### 3. Prediction of BDQ resistance class of multiple MTB VCF files
**(multiple .vcf to BDQ resistance class (*Susceptible* or *Resistance*))**

    ./bdqr-multi-VCF-predict.sh INPUT_DIR OUTPUT_DIR

`INPUT_DIR` must contain more than 1 VCF file (*<sample_id_1>*.vcf, *<sample_id_2>*.vcf, … , *<sample_id_n>*.vcf)

`OUTPUT_DIR` will contain
* the MERGED.vcf file
* the intermediate .tsv file (MERGED.tsv)
* the prediction result performed with the full model (prediction.tsv)
* the SHAP result for each isolate with the 50-feature model (shap_result_50_features_*<sample_id>*.tsv)
* the SHAP result plot for each isolate with the 50-feature model (shap_plot_50_features_*<sample_id>*.svg)

The prediction and the SHAP results output will also be displayed on the terminal.

----

### 4. Prediction of BDQ resistance class of merged MTB VCF file
**(merged.vcf to BDQ resistance class (*Susceptible* or *Resistance*))**

    ./bdqr-merge-VCF-predict.sh INPUT_DIR OUTPUT_DIR

`INPUT_DIR` must contain a merged VCF file (MERGED.vcf) with information of one or more isolates.

`OUTPUT_DIR` will contain
* the intermediate .tsv file (MERGED.tsv)
* the prediction result performed with the full model (prediction.tsv)
* the SHAP result for each isolate with the 50-feature model (shap_result_50_features_*<sample_id>*.tsv)
* the SHAP result plot for each isolate with the 50-feature model (shap_plot_50_features_*<sample_id>*.svg)

The prediction and the SHAP results output will also be displayed on the terminal.

----


## Demo runs

### 1. Prediction from FASTQ (Single isolate)

  Step 1. Create an Input directory

   	   mkdir /home/username/Input_Dir1

  Step 2. Get Data

  Download the whole genome sequencing FASTQ files of a MTB isolate run, SRR28888046 (SRR28888046_1.fastq & SRR28888046.fastq) from https://www.ebi.ac.uk/ena/browser/view/SRR28888046

  Step 3. Store these files in `Input_Dir1`

  Step 4. Create an Output directory

   	   mkdir /home/username/Output_Dir1

  Step 5. Go to the bdqr-mtb-standalone installation directory

   	   cd /home/username/bdqr-mtb-standalone/

  Step 6. Run bdqr-WGS-predict.sh

   	   ./bdqr-WGS-predict.sh /home/username/Input_Dir1/ /home/username/Output_Dir1/

`Input_Dir1` contains `SRR28888046_1.fastq`, `SRR28888046_2.fastq`

`Output_Dir1` contains -
* A folder - `SRR28888046` - which contains:
	* reference folder - reference genome and index files
	* Trim galore outputs - `SRR28888046_1_val_1.fq.gz`, `SRR28888046_2_val_2.fq.gz`, `SRR28888046_1_trimming_report.txt`, `SRR28888046_2_trimming_report.txt`
	* Bwa-mem output - `SRR28888046.bam`
	* Intermediate BAM files - `SRR28888046_fix.bam`, `SRR28888046_namesort.bam`, `SRR28888046_positionsort.bam`, `SRR28888046_markdup.bam`
	* BAM index - `SRR28888046.bam.bai`
	* Freebayes output - `SRR28888046.vcf`
	* VCF compressed - `SRR28888046.vcf.gz`
	* VCF index - `SRR28888046.vcf.gz.csi`
* `SRR28888046.tsv` - the intermediate TSV file
* `prediction.tsv` - the prediction result performed with the full model
* `shap_result_50_features_SRR28888046.tsv` - the SHAP result for SRR28888046 with the 50-feature model
* `shap_plot_50_features_SRR28888046.svg` - the SHAP result plot for SRR28888046 with the 50-feature model

----
### 2. Prediction from FASTQ (Multiple isolates)

  Step 1. Create an Input directory

    mkdir /home/username/Input_Dir2

  Step 2. Get Data

Download the whole genome sequencing FASTQ files of MTB ISOLATE runs, SRR28888046 (SRR28888046_1.fastq & SRR28888046_2.fastq) and SRR1103491 (SRR1103491_1.fastq & SRR1103491_2.fastq) from https://www.ebi.ac.uk/ena/browser/view/SRR28888046 and https://www.ebi.ac.uk/ena/browser/view/SRR1103491

  Step 3. Store these files in `Input_Dir2`

  Step 4. Create an Output directory

    mkdir /home/username/Output_Dir2

  Step 5. Go to the bdqr-mtb-standalone installation directory

    cd /home/username/bdqr-mtb-standalone/

  Step 6. Run bdqr-WGS-predict.sh

    ./bdqr-WGS-predict.sh /home/username/Input_Dir1/ /home/username/Output_Dir2/

`Input_Dir2` contains `SRR28888046_1.fastq`, `SRR28888046_2.fastq`, `SRR5341464_1.fastq`, `SRR5341464_2.fastq`

`Output_Dir2` contains -
* Two folders - `SRR28888046`, `SRR5341464` - each containing:
	* reference folder - reference genome and index files
	* Trim galore outputs - `ISOLATENAME_1_val_1.fq.gz`, `ISOLATENAME_2_val_2.fq.gz`, `ISOLATENAME_1_trimming_report.txt`, `ISOLATENAME_2_trimming_report.txt`
	* Bwa-mem output - `ISOLATENAME.bam`
	* Intermediate BAM files - `ISOLATENAME_fix.bam`, `ISOLATENAME_namesort.bam`, `ISOLATENAME_positionsort.bam`, `ISOLATENAME_markdup.bam`
	* BAM index - `ISOLATENAME.bam.bai`
	* Freebayes output - `ISOLATENAME.vcf`
	* VCF compressed - `ISOLATENAME.vcf.gz`
	* VCF index - `ISOLATENAME.vcf.gz.csi`
* `merged.vcf`
* `merged.tsv` - the intermediate TSV file
* `prediction.tsv` - the prediction result performed with the full model
* `shap_result_50_features_SRR28888046.tsv` - the SHAP result for SRR28888046 with the 50-feature model
* `shap_plot_50_features_SRR28888046.svg` - the SHAP result plot for SRR28888046 with the 50-feature model
* `shap_result_50_features_SRR5341464.tsv` - the SHAP result for SRR5341464 with the 50-feature model
* `shap_plot_50_features_SRR5341464.svg` - the SHAP result plot for SRR5341464 with the 50-feature model

----
### 3. Prediction from single VCF

  Step 1. Create an Input directory

    mkdir /home/username/Input_Dir3

  Step 2. Store a VCF file generated from variant calling of a MTB isolate based on MTB H37Rv reference genome  (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz)

  Step 3. Create an Output directory

    mkdir /home/username/Output_Dir3

  Step 4. Go to the bdqr-mtb-standalone installation directory

    cd /home/username/bdqr-mtb-standalone/

  Step 5. Run bdqr-VCF-predict.sh

    ./bdqr-VCF-predict.sh /home/username/Input_Dir3/ /home/username/Output_Dir3/

`Input_Dir3` contains `SRR28888046.vcf`

`Output_Dir3` contains -
* `SRR28888046.tsv` - the intermediate TSV file
* `prediction.tsv` - the prediction result performed with the full model
* `shap_result_50_features_SRR28888046.tsv` - the SHAP result for SRR28888046 with the 50-feature model
* `shap_plot_50_features_SRR28888046.svg` - the SHAP result plot for SRR28888046 with the 50-feature model

----
### 4. Prediction from multiple VCFs

  Step 1. Create an Input directory

    mkdir /home/username/Input_Dir4

  Step 2. Store multiple VCF files generated from variant calling of  MTB isolates based on MTB H37Rv reference genome  (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz)

  Step 3. Create an Output directory

    mkdir /home/username/Output_Dir4

  Step 4. Go to the bdqr-mtb-standalone installation directory

    cd /home/username/bdqr-mtb-standalone/

  Step 5. Run bdqr-multi-VCF-predict.sh

    ./bdqr-multi-VCF-predict.sh /home/username/Input_Dir4/ /home/username/Output_Dir4/

`Input_Dir4` contains `SRR28888046.vcf`, `SRR5341464.vcf`

`Output_Dir4` contains -
* Compressed VCFs - `SRR28888046.vcf.gz` and `SRR5341464.vcf.gz`
* VCF indices - `SRR28888046.vcf.gz.csi` and `SRR5341464.vcf.gz.csi`
* `merged.vcf`
* `merged.tsv` - the intermediate TSV file
* `prediction.tsv` - the prediction result performed with the full model
* `shap_result_50_features_SRR28888046.tsv` - the SHAP result for SRR28888046 with the 50-feature model
* `shap_plot_50_features_SRR28888046.svg` - the SHAP result plot for SRR28888046 with the 50-feature model
* `shap_result_50_features_SRR5341464.tsv` - the SHAP result for SRR5341464 with the 50-feature model
* `shap_plot_50_features_SRR5341464.svg` - the SHAP result plot for SRR5341464 with the 50-feature model

----
### 5. Prediction from merged VCFs

  Step 1. Create an Input directory

    mkdir /home/username/Input_Dir5

  Step 2. Store a merged VCF file generated by merging multiple VCF files from variant calling of  MTB isolates based on MTB H37Rv reference genome  (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz)

  Step 3. Create an Output directory

    mkdir /home/username/Output_Dir5

  Step 4. Go to the bdqr-mtb-standalone installation directory

    cd /home/username/bdqr-mtb-standalone/

  Step 5. Run bdqr-merge-predict.sh

    ./bdqr-merge-predict.sh /home/username/Input_Dir5/ /home/username/Output_Dir5/

`Input_Dir5` contains `merged.vcf`

`Output_Dir5` contains -
* `merged.tsv` - the intermediate TSV file
* `prediction.tsv` - the prediction result performed with the full model
* `shap_result_50_features_SRR28888046.tsv` - the SHAP result for SRR28888046 with the 50-feature model
* `shap_plot_50_features_SRR28888046.svg` - the SHAP result plot for SRR28888046 with the 50-feature model
* `shap_result_50_features_SRR5341464.tsv` - the SHAP result for SRR5341464 with the 50-feature model
* `shap_plot_50_features_SRR5341464.svg` - the SHAP result plot for SRR5341464 with the 50-feature model

## Team

**Stuti Ghosh, Sudipto Bhattacharjee, and Sudipto Saha**

## Disclaimer

The scripts were tried and tested on the Ubuntu Operating system.

This tool is strictly for Research Use Only. By using this tool the user acknowledges no intended medical purpose such as patient diagnosis.

