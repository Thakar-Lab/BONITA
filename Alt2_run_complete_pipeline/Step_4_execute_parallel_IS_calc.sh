#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J localSearch
#SBATCH -o localSearch.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=1G

sbatch submit_IS_prallel_local_search.sh IS_1_setup.pickle IS_individuals_1_setup.pickle
sbatch submit_IS_prallel_local_search.sh IS_2_setup.pickle IS_individuals_2_setup.pickle
sbatch submit_IS_prallel_local_search.sh IS_3_setup.pickle IS_individuals_3_setup.pickle
sbatch submit_IS_prallel_local_search.sh IS_4_setup.pickle IS_individuals_4_setup.pickle
sbatch submit_IS_prallel_local_search.sh IS_5_setup.pickle IS_individuals_5_setup.pickle