#!/bin/bash
#SBATCH --job-name=sintese_proteica
#SBATCH --output=saidas/sintese_proteica.txt
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024
#SBATCH --time=00:05:00
#SBATCH --partition=normal

mpirun executaveis/codons