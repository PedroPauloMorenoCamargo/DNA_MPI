# DNA_MPI

This project processes large DNA sequence files, performing the following tasks using MPI:
1. Filters the original DNA sequence files
2. Counts each type of DNA base (adenine, cytosine, guanine, and thymine).
3. Transcribes the DNA sequence to RNA.
4. Counts the number of proteins (polypeptides) in the RNA sequence.
5. Extracts and outputs the codons from the original DNA sequence.

## How to Run

Follow these steps to run the project. All commands should be executed from the root directory of the repository.

### 1. Download the Files

Use the following command to obtain the DNA sequence files:

```bash
/entradas/download.sh
```

### 2. Compile

Use the following command to compile the program:

```bash
mpic++ -fopenmp -std=c++17 -O3 codigos/XXXXX.cpp -o executaveis/XXXXX
```

### 3.Run

Use the following command to run the program:

```bash
sbatch slurms/XXXXXX.slurm
```