#!/bin/sh

#SBATCH --job-name=run_rhime_2014_monthly_%a
#SBATCH --output=run_rhime_2014_monthly_%a.out
#SBATCH --time=72:00:00
#SBATCH --mem=100gb
#SBATCH --array=3,4,5
#SBATCH --partition=cpu-long

year="2014"
year_end="2015"

start=("${year}-01-01" "${year}-02-01" "${year}-03-01" "${year}-04-01" "${year}-05-01" "${year}-06-01" "${year}-07-01" "${year}-08-01" "${year}-09-01" "${year}-10-01" "${year}-11-01" "${year}-12-01")
end=("${year}-02-01" "${year}-03-01" "${year}-04-01" "${year}-05-01" "${year}-06-01" "${year}-07-01" "${year}-08-01" "${year}-09-01" "${year}-10-01" "${year}-11-01" "${year}-12-01" "${year_end}-01-01")

echo "Job ID number:" $SLURM_ARRAY_TASK_ID
echo "start date:" "${start[$SLURM_ARRAY_TASK_ID]}"

python /home/users/alice.ramsden/rhime-co2-inversions/run_rhime_2014_monthly.py "${start[$SLURM_ARRAY_TASK_ID]}" "${end[$SLURM_ARRAY_TASK_ID]}"
