#!/bin/bash
#SBATCH --job-name=hyp2-19-lambda
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=2000
#SBATCH --time=36:00:00


module load R

Rscript hyp2-lambda-2019.R