# BONITA

BONITA- Boolean Omics Network Invariant-Time Analysis is a package for the inference of Boolean rules and pathway analysis on omics data. It can be applied to help uncover underlying relationships in biological data. Please see our publication for more information: [insert eventual publication here]. 

**Authors**: _Rohith Palli, Mukta G. Palshikar, and Juilee Thakar_\
**Maintainer**: Please contact Rohith Palli at rohith_palli@urmc.rochester.edu

# Installation
BONITA is designed for use with distributed computing systems. A desktop version installation and use guide are provided in the wiki (see bar across top of github page). Necessary SLURM commands are included. If users are having trouble translating to PBS or other queueing standards for their computing environment, please contact maintainer Rohith Palli at rohith_palli@urmc.rochester.edu. 

## Install Python
Please create a virtual environment with the latest version of Intel Python 2 and install the following packages (all accessible by pip):
* networkx==1.11
* pandas
* requests
* deap
* lxml
* bs4
* seaborn\
The software was tested with Intel Python 2.7.12 and the following package versions: 
networkx==1.11, pandas==0.19.0, requests==2.11.1, deap==1.0.2, 
bs4==0.0.1, lxml==4.1.1, seaborn==0.8.1. In our experience, only the networkx version affects BONITA functionality. 

## Install BONITA
You can download and use BONITA in one of two ways:
1. Download a zipped folder containing all the files you need (github download link in green box above and to the right)\
1. Clone this git repository in the folder of your choice using the command `git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY`\

Next, the C code must be compiled using the make file. Simply type make while in the BONITA folder. \
`make`\
Now you have a fully functional distribution of BONITA! Time to gather your data and get started. 

# Usage

You will need the following files to run BONITA:
* omics data as a plaintext table (csv, tsv, or similar) with the first row containing a holder for gene symbol column then sample names and subsequent rows containing gene symbol in first column and column-normalized (rpm or rpkm in transcriptomics) abundance measures in other columns. 
* gmt file with list of KEGG pathways to be considered (can be downloaded from msigdb)
* matrix of conditions with each line representing a sample and the first column containing the names of the samples and subsequent columns describing 1/0 if the sample is part of that condition or not. 
* list of contrasts you would like to run with each contrast on a single line

There are three main steps in BONITA: prepare pathways for rule inference, rule inference, and pathway analysis. All necessary files for an example run are provided in the pathway_analysis folder within experiments folder. The preparation step requires internet access to access the KEGG API. 

**Step 1: Pathway prepearation**

**This step requires internet access.** BONITA needs omics data, gmt file, and an indication of what character is used to separate columns in the file. For example, a traditional comma separated value file (csv) would need BONITA input "-sep ,". Since tab can't be passed in as easily, a -t command will automatically flag tab as the separator. The commands are below:
comma separated: `python pathway_analysis_setup.py Your_omics_data Your_gmt_file -sep ,`
tab separated: `python pathway_analysis_setup.py -t Your_omics_data Your_gmt_file`

**Step 2: Rule inference**

Simply run the script find_rules_pathway_analysis.sh which will automatically submit appropriate jobs to SLURM queue.

**Step 3: Pathway Analysis**
To accomplish this, the proper inputs must be provided to pathway_analysis_score_pathways.py

`python pathway_analysis_score_pathways.py Your_omics_data Your_condition_matrix Your_desired_contrasts -sep Separator_used_in_gmt_and_omics_data`\
If your files are tab seperated, then the following command can be used: `python pathway_analysis_score_pathways.py -t Your_omics_data Your_condition_matrix Your_desired_contrasts`
