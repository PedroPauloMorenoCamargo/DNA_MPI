#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --output=saidas/preprocess.txt
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=1024
#SBATCH --time=00:05:00
#SBATCH --partition=normal

mpirun executaveis/preprocess