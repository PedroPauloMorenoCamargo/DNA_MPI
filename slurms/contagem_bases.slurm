#!/bin/bash
#SBATCH --job-name=contagem_bases
#SBATCH --output=saidas/contagem_bases.txt
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024
#SBATCH --time=00:05:00
#SBATCH --partition=normal

mpirun executaveis/contagem_bases