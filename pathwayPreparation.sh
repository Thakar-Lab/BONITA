#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J bonitaStep1
#SBATCH -o bonitaStep1.log
#SBATCH -t 1:00:00

module load intelpython/2.7.12

#Option 1: On a gmt of human pathways BONITA needs omics data, gmt file, and an indication of what character is used to separate columns in the file
#comma separated
#python pathway_analysis_setup.py -gmt Your_gmt_file -sep , Your_omics_data
#tab separated
python pathway_analysis_setup.py -t -gmt Your_gmt_file Your_omics_data

#Option 2: On all KEGG pathways for any organism BONITA needs omics data, organism code, and an indication of what character is used to separate columns in the file.
#comma separated, human: 
#python pathway_analysis_setup.py -org hsa -sep , Your_omics_data #MOST COMMON USAGE
#comma separated, mouse: 
#python pathway_analysis_setup.py -org mmu -sep , Your_omics_data
#tab separated: 
#python pathway_analysis_setup.py -t -org Your_org_code Your_omics_data

#Option 3: On a list of KEGG pathways for any organism BONITA needs omics data, organism code, the list of pathways, and an indication of what character is used to separate columns in the file. 
#comma separated, human
#python pathway_analysis_setup.py -org hsa -sep , -paths Your_pathway_list Your_omics_data
#comma separated, mouse
#python pathway_analysis_setup.py -org mmu -sep , -paths Your_pathway_list Your_omics_data
#tab separated
#python pathway_analysis_setup.py -t -org Your_org_code -paths Your_pathway_list Your_omics_data
