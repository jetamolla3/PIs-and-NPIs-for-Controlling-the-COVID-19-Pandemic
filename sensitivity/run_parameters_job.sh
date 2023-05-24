#!/bin/bash
#SBATCH --time=32:00:00
#SBATCH --mem=96
#SBATCH --array=429,459,564,644,951

parameters_file="waning.csv"
column_num=$((${SLURM_ARRAY_TASK_ID}+1))
entry=$(sed "${column_num}q;d" $parameters_file)

IFS="," read -r -a array <<< $entry
python sensitivity.py "${array[@]}"

printf "%s," "${array[@]}" >> "output.txt"



