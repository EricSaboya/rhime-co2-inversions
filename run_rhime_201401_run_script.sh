#!/bin/sh

#SBATCH --job-name=run_rhime_201401
#SBATCH --output=run_rhime_201401.out
#SBATCH --time=6:00:00
#SBATCH --mem=30gb

python /home/users/alice.ramsden/rhime-co2-inversions/run_rhime_201401.py
