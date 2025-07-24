#!/bin/sh

#SBATCH --job-name=run_rhime_2014_weekly_%a
#SBATCH --output=run_rhime_2014_weekly_%a.out
#SBATCH --time=72:00:00
#SBATCH --mem=100gb
#SBATCH --array=0-3
#SBATCH --partition=cpu-long

year="2014"
year_end="2015"

start=("${year}-03-01" "${year}-03-08" "${year}-03-15" "${year}-03-22" "${year}-04-01" "${year}-04-08" "${year}-04-15" "${year}-04-22")
end=("${year}-03-08" "${year}-03-15" "${year}-03-22" "${year}-04-01" "${year}-04-08" "${year}-04-15" "${year}-04-22" "${year}-05-01")

echo "Job ID number:" $SLURM_ARRAY_TASK_ID
echo "start date:" "${start[$SLURM_ARRAY_TASK_ID]}"

python /home/users/alice.ramsden/rhime-co2-inversions/run_rhime_201401.py "${start[$SLURM_ARRAY_TASK_ID]}" "${end[$SLURM_ARRAY_TASK_ID]}"
