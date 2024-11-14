#!/bin/bash

#SBATCH -N 1
#SBATCH -n 5
#SBATCH --mem=20GB
#SBATCH -t 1-
#SBATCH --output=Fit_EndoPRS/Refit.txt

module load r

Rscript 02.Refit_endoPRS.R endoPRS_model_results/

