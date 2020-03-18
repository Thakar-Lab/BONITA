# BONITA

BONITA- Boolean Omics Network Invariant-Time Analysis is a package for the inference of Boolean rules and pathway analysis on omics data. It can be applied to help uncover underlying relationships in biological data. Please see our [publication](https://doi.org/10.1371/journal.pcbi.1007317) for more information. 

**Authors**: _Rohith Palli, Mukta G. Palshikar, and Juilee Thakar_\
**Maintainer**: Please contact Rohith Palli at rohith_palli@urmc.rochester.edu

# Citation
We would appreciate the citation of our manuscript describing BONITA, below, for any use of our code. 

Palli R, Palshikar MG, Thakar J (2019) Executable pathway analysis using ensemble discrete-state modeling for large-scale data. PLoS Comput Biol 15(9): e1007317. (https://doi.org/10.1371/journal.pcbi.1007317)

# Installation
BONITA is designed for use with distributed computing systems. A desktop version installation and use guide are provided in the wiki (https://github.com/Thakar-Lab/BONITA/wiki). Necessary SLURM commands are included. If users are having trouble translating to PBS or other queueing standards for their computing environment, please contact maintainer Rohith Palli at rohith_palli@urmc.rochester.edu. 

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

**Step 1: Pathway preparation**

**This step requires internet access.** 

**There are three ways to complete this process: 1) on a gmt of human pathways, 2) on all KEGG pathways for any organism, or 3) on a list of KEGG pathways for any organism**

**Only Option 1 was used and tested in our manuscript. Caution should be exercised in interpreting results of other two methods. At a minimum, graphmls with impact scores and relative abundance should be examined before drawing conclusions about pathway differences.**

**Option 1: On a gmt of human pathways** 
BONITA needs omics data, gmt file, and an indication of what character is used to separate columns in the file. For example, a traditional comma separated value file (csv) would need BONITA input "-sep ,". Since tab can't be passed in as easily, a -t command will automatically flag tab as the separator. The commands are below:

comma separated: `python pathway_analysis_setup.py -gmt Your_gmt_file -sep , Your_omics_data `

tab separated: `python pathway_analysis_setup.py -t  -gmt Your_gmt_file Your_omics_data`

**Option 2: On all KEGG pathways for any organism** 
BONITA needs omics data, organism code, and an indication of what character is used to separate columns in the file. For example, a traditional comma separated value file (csv) would need BONITA input "-sep ,". Since tab can't be passed in as easily, a -t command will automatically flag tab as the separator. A three letter organism code from KEGG must be provided (lower case). Example codes include mmu for mouse and hsa for human. The commands are below:
comma separated: `python pathway_analysis_setup.py -org Your_org_code -sep , Your_omics_data `

comma separated, human: `python pathway_analysis_setup.py -org hsa -sep , Your_omics_data `

comma separated, mouse: `python pathway_analysis_setup.py -org mmu -sep , Your_omics_data `

tab separated: `python pathway_analysis_setup.py -t  -gmt Your_org_code Your_omics_data`

**Option 3: On a list of KEGG pathways for any organism** 
BONITA needs omics data, organism code, the list of pathways, and an indication of what character is used to separate columns in the file. For example, a traditional comma separated value file (csv) would need BONITA input "-sep ,". Since tab can't be passed in as easily, a -t command will automatically flag tab as the separator. A three letter organism code from KEGG must be provided (lower case). Example codes include mmu for mouse and hsa for human. The list of pathways must include the 5 digit pathway identifier, must be seperated by commas, and must not include any other numbers. An example paths.txt is included in the inputData folder. The commands are below:
comma separated: `python pathway_analysis_setup.py -org Your_org_code -sep , -paths Your_pathway_list Your_omics_data `

comma separated, human: `python pathway_analysis_setup.py -org hsa -sep , -paths Your_pathway_list Your_omics_data `

comma separated, mouse: `python pathway_analysis_setup.py -org mmu -sep , -paths Your_pathway_list Your_omics_data `

tab separated: `python pathway_analysis_setup.py -t  -gmt Your_org_code -paths Your_pathway_list Your_omics_data`

**Step 2: Rule inference**

Simply run the script find_rules_pathway_analysis.sh which will automatically submit appropriate jobs to SLURM queue 

**Step 3: Pathway Analysis**
To accomplish this, the proper inputs must be provided to pathway_analysis_score_pathways.py

`python pathway_analysis_score_pathways.py Your_omics_data Your_condition_matrix Your_desired_contrasts -sep Separator_used_in_gmt_and_omics_data`\
If your files are tab separated, then the following command can be used: `python pathway_analysis.py -t Your_omics_data Your_condition_matrix Your_desired_contrasts`
