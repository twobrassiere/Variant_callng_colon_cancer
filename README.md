# Variant_callng_colon_cancer
This code of review by [variant calling pipeline using GATK4 and nextflow](https://github.com/gencorefacility/variant-calling-pipeline-gatk4) and download [whole exome sequencing for metastatic colorectal cancer](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA726023).

# Part 0: Setup
1. [Anaconda](https://www.anaconda.com/products/distribution) or [Colab-pro](https://colab.research.google.com/?hl=zh-tw) write automatically wxs analysis pipeline code by jupyter notebook.
2. [Sra-toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) download and fastq-dump by whole exome sequencing for metastatic colorectal cancer.
3. [Cutadapt](https://anaconda.org/bioconda/cutadapt) trimmed wxs sequence data .
4. [Docker](https://desktop.docker.com/mac/main/amd64/Docker.dmg?utm_source=docker&utm_medium=webreferral&utm_campaign=dd-smartbutton&utm_location=module)(>=19.03) or [Singularity(by HPC)](https://docs.sylabs.io/guides/latest/user-guide/) automatically run bioinformation container .
5. [Neoflow](https://github.com/bzhanglab/neoflow) and [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) automatically  bulid sequence anaylsis pipeline .
6. [GATK](https://github.com/broadinstitute/gatk/releases)(=4.2.3) package variant call .

# Part 1: Download whole exome sequencing for metastatic colorectal cancer(PRJNA726023) and trimmed by fold of [Jupyter notebook](https://github.com/twobrassiere/Automatically-analysis-whole-exome-sequencing-of-colon-cancer-dataset/tree/main/Jupyter%20notebook) .
   Step 1 : Download accession list by [SRA run selector](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=5&WebEnv=MCID_63a51fbfe936a5469741cbf0&o=acc_s%3Aa) 
   
   Step 2 : Automatically download sequencing code by [Automatic download sequence by sra toolkit.ipynb](# Part 1: Download whole exome sequencing for metastatic colorectal cancer(PRJNA726023) and trimmed by fold of [Jupyter notebook](https://github.com/twobrassiere/Automatically-analysis-whole-exome-sequencing-of-colon-cancer-dataset/tree/main/Jupyter%20notebook) .
   Step 1 : Download accession list by [SRA run selector](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=5&WebEnv=MCID_63a51fbfe936a5469741cbf0&o=acc_s%3Aa) 
   
   Step 2 : Automatically download sequencing code by [Automatic download sequence by sra toolkit.ipynb].
   
   Step 3 : Search adapter by [Illumina Adapter Sequences](https://support-docs.illumina.com/SHARE/AdapterSeq/illumina-adapter-sequences.pdf)
   
   Step 4 : Automatically trimmed wxs sequence data by [Use pair-end sequnce triming by cutapat.ipynb](https://github.com/twobrassiere/Automatically-analysis-whole-exome-sequencing-of-colon-cancer-dataset/blob/main/Jupyter%20notebook/Use%20pair-end%20sequnce%20%20triming%20by%20cutapat.ipynb)).
   
   Step 3 : Search adapter by [Illumina Adapter Sequences](https://support-docs.illumina.com/SHARE/AdapterSeq/illumina-adapter-sequences.pdf)
   
   Step 4 : Automatically trimmed wxs sequence data by [Use pair-end sequnce triming by cutapat.ipynb](https://github.com/twobrassiere/Automatically-analysis-whole-exome-sequencing-of-colon-cancer-dataset/blob/main/Jupyter%20notebook/Use%20pair-end%20sequnce%20%20triming%20by%20cutapat.ipynb)
