# Variant_callng_colon_cancer
This code of review by [variant calling pipeline using GATK4 and nextflow](https://github.com/gencorefacility/variant-calling-pipeline-gatk4) and download [whole exome sequencing for metastatic colorectal cancer](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA726023).

## Part 0: Setup
1. [Anaconda](https://www.anaconda.com/products/distribution) or [Colab-pro](https://colab.research.google.com/?hl=zh-tw) write automatically wxs analysis pipeline code by jupyter notebook.
2. [Sra-toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) download and fastq-dump by whole exome sequencing for metastatic colorectal cancer.
3. [Cutadapt](https://anaconda.org/bioconda/cutadapt) trimmed wxs sequence data .
4. [Docker](https://desktop.docker.com/mac/main/amd64/Docker.dmg?utm_source=docker&utm_medium=webreferral&utm_campaign=dd-smartbutton&utm_location=module)(>=19.03) or [Singularity(by HPC)](https://docs.sylabs.io/guides/latest/user-guide/) automatically run bioinformation container .
5. [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) automatically  bulid sequence anaylsis pipeline .
6. [GATK](https://github.com/broadinstitute/gatk/releases)(=4.2.3) package variant call .

## Part 1: Download whole exome sequencing for metastatic colorectal cancer(PRJNA726023) and trimmed by fold of [Jupyter notebook](https://github.com/twobrassiere/Variant_callng_colon_cancer/tree/main/Jupyter%20notebook) 
Step 1 : Download accession list by [SRA run selector](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=5&WebEnv=MCID_63a51fbfe936a5469741cbf0&o=acc_s%3Aa) 

Step 2 : Automatically download sequencing code by [Automatic download sequence by sra toolkit.ipynb](https://github.com/twobrassiere/Variant_callng_colon_cancer/blob/main/Jupyter%20notebook/Automatic%20download%20sequence%20by%20sra%20toolkit.ipynb)

Step 3 : Search adapter by [Illumina Adapter Sequences](https://support-docs.illumina.com/SHARE/AdapterSeq/illumina-adapter-sequences.pdf)

Step 4 : Automatically trimmed wxs sequence data by [Use pair-end sequnce triming by cutapat.ipynb](https://github.com/twobrassiere/Variant_callng_colon_cancer/blob/main/Jupyter%20notebook/Use%20pair-end%20sequnce%20%20triming%20by%20cutapat.ipynb)
   
## We bulid workflow of  variant calling pipeline using nextflow![iamge](https://github.com/twobrassiere/Variant_callng_colon_cancer/blob/main/workflow.jpg)

## Part 2 : variant calling using nexflow
```sh
 $ ./nextflow run main.nf --help
=========================================
neoflow => WXS anylsis
=========================================
Usage:
nextflow run main.nf
Arguments:
  --reads                     Reads data in fastq.gz or fastq format. For example, "*_{1,2}.fastq.gz"
  --ref_dir                   Reference  sequence folder
  --seqtype                   Read type, dna or rna. Default is dna.
  --singleEnd                 Single end or not, default is false (pair end reads)
  --cpu                       The number of CPUs, default is 4.
  --vcf_dir                   Folder of variant file , default is "./"
  --help                      Print help message
 ```
The  output of `main.nf` is a txt format file containing vcf  for a sample and zip file . This file is generated by [vcf file](https://github.com/twobrassiere/Variant_callng_colon_cancer/tree/main/vcf)
   
#### Example by servse
```sh
 ./nextflow run  main.nf
 --reads "./fastq_trimmed/SRR14463457_pass_{1,2}_trimmed.fastq.gz" \ 
 --ref_dir ./reference_genome/GRCh37 \
 --vcf_dir ./vcf \
 --cpu 14
```

#### Example by Taiwania 3(HPC)
```sh
#!/usr/bin/sh
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J Job_name         # Job name
#SBATCH -p ngs186G        # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 28              # 使用的數 請參考Queue資源設定
#SBATCH --mem=186g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log          # Path to the standard output file
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=elephantliu@nycu.edu.tw    # email
#SBATCH --mail-type=BEGIN,END              # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, ALL

module load  biology/Samtools/1.15.1
module load  biology/OpenJDK/17.0.2+8
module load biology/BWA/0.7.17
module load  biology/GATK/4.2.3.0
 ./nextflow run main.nf 
 --reads "./fastq_trimmed/SRR14463457_pass_{1,2}_trimmed.fastq.gz" \
 --ref_dir ./reference_genome/GRCh37 \
 --vcf_dir ./vcf \
 --cpu 14
 ```
# Acknowledgements 
We appreciate Taiwania 3 by TWCC to variant call of colon cancer 
