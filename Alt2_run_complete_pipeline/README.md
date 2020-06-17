#Parallelized version of BONITA's local search

## Input:
One network, as a graphml
Expression dataset
USER MUST PASTE THESE FILE NAMES INTO THE STEP 1 SCRIPT
ALSO NEED TO CHANGE 'temp_series1_net_sss' name in the parallel_local_search_find_node_importance to match file name of input graphml
Finally, paste the file names for data, matrix of sample groups, and contrasts in text file into Step_5_compile_final_results.sh in place of "dataName.csv" "matrix.name" "diff.file"
## Step 1: Set up data and network
Runs pathway analysis setup then starts jobs for local search
Need to change names of  "dataName.csv" and "pathway.graphml" in the file to the names of the files to be used in the run
"sbatch Step_1_create_GA_model_and_set_up_parallel_local_search.sh"

## Step 2: Run parallel local search
"sbatch Step_2_parallel_local_search_find_node_importance.sh"


## Step 3: Clean up results from local search and setup IS calculations
Cleans up the results of step 1
 "sbatch Step_3_cleanup_local_search.sh"

## Step 4: Run Impact Score calculations
Runs Impact score caluclations
 "sbatch Step_4_execute_parallel_IS_calc.sh"

## Step 5: Integrate findings with pathway analysis scripts
Compile the fineal results together
"sbatch Step_5_compile_final_results.sh"
