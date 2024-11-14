#!/bin/bash

#SBATCH -N 1
#SBATCH -n 5
#SBATCH --mem=20GB
#SBATCH --array=1-25
#SBATCH -t 1-
#SBATCH --output=Fit_EndoPRS/%A_%a.txt

module load r

Rscript 01.endoPRS_gridsearch.R 0.01 $SLURM_ARRAY_TASK_ID endoPRS_model_results/
Rscript 01.endoPRS_gridsearch.R 1e-4 $SLURM_ARRAY_TASK_ID endoPRS_model_results/
Rscript 01.endoPRS_gridsearch.R 1e-6 $SLURM_ARRAY_TASK_ID endoPRS_model_results/
