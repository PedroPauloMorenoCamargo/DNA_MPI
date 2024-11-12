#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>
#include <filesystem>
#include <cctype>

namespace fs = std::filesystem;
using namespace std;

void count_bases(const string& input_name, array<long long, 4>& local_counts) {
    // Abre o arquivo binário
    ifstream input(input_name, ios::binary);
    if (!input) {
        cerr << "Error opening input file: " << input_name << endl;
        return;
    }

    // Cria o buffer de leitura
    const size_t buffer_size = 8192;
    vector<char> buffer(buffer_size);

    while (input.read(buffer.data(), buffer_size) || input.gcount()) {
        //Bytes Lidos
        size_t bytes_read = input.gcount();

        // Contagem de bases
        long long count_a = 0, count_t = 0, count_c = 0, count_g = 0;

        //Conta Bases paralelamente
        #pragma omp parallel for reduction(+:count_a, count_t, count_c, count_g)
        for (size_t i = 0; i < bytes_read; ++i) {
            char ch = buffer[i];
            switch (ch) {
                case 'a': count_a++; break;
                case 't': count_t++; break;
                case 'c': count_c++; break;
                case 'g': count_g++; break;
                default: break; 
            }
        }

        //Updata a contagem local
        local_counts[0] += count_a;
        local_counts[1] += count_t;
        local_counts[2] += count_c;
        local_counts[3] += count_g;
    }
}

int main(int argc, char *argv[]) {
    //Começa MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Começa Tempo
    double start_time = MPI_Wtime();

    // Lista Inputs
    vector<string> input_files = {
        "entradas_filtradas/chr1.subst.filtered", "entradas_filtradas/chr2.subst.filtered", "entradas_filtradas/chr3.subst.filtered",
        "entradas_filtradas/chr4.subst.filtered", "entradas_filtradas/chr5.subst.filtered", "entradas_filtradas/chr6.subst.filtered",
        "entradas_filtradas/chr7.subst.filtered", "entradas_filtradas/chr8.subst.filtered", "entradas_filtradas/chr9.subst.filtered",
        "entradas_filtradas/chr10.subst.filtered", "entradas_filtradas/chr11.subst.filtered", "entradas_filtradas/chr12.subst.filtered",
        "entradas_filtradas/chr13.subst.filtered", "entradas_filtradas/chr14.subst.filtered", "entradas_filtradas/chr15.subst.filtered",
        "entradas_filtradas/chr16.subst.filtered", "entradas_filtradas/chr17.subst.filtered", "entradas_filtradas/chr18.subst.filtered",
        "entradas_filtradas/chr19.subst.filtered", "entradas_filtradas/chr20.subst.filtered", "entradas_filtradas/chr21.subst.filtered",
        "entradas_filtradas/chr22.subst.filtered"
    };

    // Divide os arquivos por processos
    int num_files = input_files.size();
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

    //Contagem Local de Cada Base
    array<long long, 4> local_counts = {0, 0, 0, 0};

    // Conta as Bases dos Arquivos
    for (int i = start_idx; i < end_idx && i < num_files; ++i) {
        string input_name = input_files[i];
        cout << "Process " << rank << " processing file: " << input_name << endl;
        count_bases(input_name, local_counts);
    }

    //Reduz as contas locais
    array<long long, 4> global_counts = {0, 0, 0, 0};
    MPI_Reduce(local_counts.data(), global_counts.data(), 4, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    //Tempo Percorrido
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    // Printa resultados
    if (rank == 0) {
        cout << "Total counts:" << endl;
        cout << "A: " << global_counts[0] << endl;
        cout << "T: " << global_counts[1] << endl;
        cout << "C: " << global_counts[2] << endl;
        cout << "G: " << global_counts[3] << endl;
        cout << "Total execution time: " << elapsed_time << " seconds." << endl;
    }

    // Finaliza MPI
    MPI_Finalize();
    return 0;
}
