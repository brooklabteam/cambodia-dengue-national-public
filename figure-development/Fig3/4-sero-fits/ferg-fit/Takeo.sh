#!/bin/bash
#SBATCH --job-name=Takeo
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --time=36:00:00


module load R

Rscript fig2-model-fitting-province-Takeo.R