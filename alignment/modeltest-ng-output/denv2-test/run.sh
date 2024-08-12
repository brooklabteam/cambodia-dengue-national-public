#!/bin/bash
#SBATCH --job-name=denv2-modtest
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2000

module load flex/2.6.4
module load vim/8.1  
module load openmpi/3.1.2
module load cmake/3.15 
module load python/cpython-3.7.0
module load gcc/10.2.0
module load emacs/26
module load java/1.8


~/modeltest-ng/modeltest-ng-static -i ~/modeltest-ng/denv2-test/DENV2aligned.fasta -t ml -p 1
