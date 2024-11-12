#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>
#include <filesystem>

namespace fs = std::filesystem;
using namespace std;

size_t count_complete_proteins_in_rna(const string& rna_file) {
    // Abre arquivo de RNA
    ifstream input(rna_file, ios::binary);
    if (!input) {
        cerr << "Error opening RNA file: " << rna_file << endl;
        return 0;
    }

    // Cria buffer
    const size_t buffer_size = 8192;
    vector<char> buffer(buffer_size);

    // Conta proteinas
    size_t protein_count = 0;

    // Quando não for de 3 em 3 pega os caracteres passados
    string prev_chars;

    bool in_protein = false;

    // Codon Parcial
    string partial_codon;

    while (input.read(buffer.data(), buffer_size) || input.gcount()) {
        // Numeros de bytes Lidos
        size_t bytes_read = input.gcount();
        string data = partial_codon + string(buffer.data(), bytes_read);

        size_t data_length = bytes_read+partial_codon.length();

        size_t num_codons = data_length / 3;

        // Processa sequencialmente os codons
        for (size_t i = 0; i < num_codons; ++i) {
            size_t idx = i * 3;
            string codon = data.substr(idx, 3);

            if (!in_protein) {
                //Codon de inicio
                if (codon == "aug") {
                    in_protein = true;
                }
            } else {
                //Codon de parada
                if (codon == "uaa" || codon == "uag" || codon == "uga") {
                    in_protein = false;
                    protein_count++;
                }
            }
        }
        size_t remaining_chars = data_length % 3;
        partial_codon = data.substr(data_length - remaining_chars, remaining_chars);
    }
    return protein_count;
}

int main(int argc, char *argv[]) {
    // Inicializa MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time = MPI_Wtime();

    // Inputs
    vector<string> rna_files = {
        "rna/chr1.rna", "rna/chr2.rna", "rna/chr3.rna",
        "rna/chr4.rna", "rna/chr5.rna", "rna/chr6.rna",
        "rna/chr7.rna", "rna/chr8.rna", "rna/chr9.rna",
        "rna/chr10.rna", "rna/chr11.rna", "rna/chr12.rna",
        "rna/chr13.rna", "rna/chr14.rna", "rna/chr15.rna",
        "rna/chr16.rna", "rna/chr17.rna", "rna/chr18.rna",
        "rna/chr19.rna", "rna/chr20.rna", "rna/chr21.rna",
        "rna/chr22.rna"
    };

    // Divide arquivos processos
    int num_files = rna_files.size();
    int files_per_proc = num_files / size;
    int remainder = num_files % size;
    int start_idx, end_idx;

    if (rank < remainder) {
        files_per_proc++;
        start_idx = rank * files_per_proc;
    } else {
        start_idx = rank * files_per_proc + remainder;
    }
    end_idx = start_idx + files_per_proc;

    // Variavel de contagem total de cada processo
    size_t total_protein_count = 0;

    // Processa os arquivos
    for (int i = start_idx; i < end_idx && i < num_files; ++i) {
        string rna_file = rna_files[i];

        size_t protein_count = count_complete_proteins_in_rna(rna_file);

        cout << "Process " << rank << " found " << protein_count << " complete proteins in file: " << rna_file << endl;

        total_protein_count += protein_count;
    }

    // Pega o valor de todos os processos
    size_t global_protein_count = 0;
    MPI_Reduce(&total_protein_count, &global_protein_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    double max_time;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Total complete proteins found: " << global_protein_count << endl;
        cout << "Total execution time: " << max_time << " seconds." << endl;
    }

    MPI_Finalize();
    return 0;
}
