#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>
#include <filesystem>
#include <array>


namespace fs = std::filesystem;
using namespace std;

//Função lambda para inicializar a lookup table de transcrição
const array<char, 256> DNA_TO_RNA = []() -> std::array<char, 256> {
    array<char, 256> table = {};
    table.fill('n');
    
    table['a'] = 'u';
    table['t'] = 'a';
    table['c'] = 'g';
    table['g'] = 'c';
    
    return table;
}();


void transcribe_dna_to_rna(const string& input_name, const string& output_name) {
    // I/O
    ifstream input(input_name, ios::binary);
    ofstream output(output_name, ios::binary);
    if (!input) {
        cerr << "Error opening input file: " << input_name << endl;
        return;
    }
    if (!output) {
        cerr << "Error opening output file: " << output_name << endl;
        return;
    }

    // Buffer
    const size_t buffer_size = 1 << 20;
    vector<char> buffer(buffer_size);

    while (input.read(buffer.data(), buffer_size) || input.gcount()) {
        size_t bytes_read = input.gcount();

        // Transcricao
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < bytes_read; ++i) {
            buffer[i] = DNA_TO_RNA[static_cast<unsigned char>(buffer[i])];
        }

        // Escreve a file
        output.write(buffer.data(), bytes_read);
    }
}

int main(int argc, char *argv[]) {
    // Inicializa MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time = MPI_Wtime();

    // Inputs
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

    // Distribui os arquivos entre proocessos
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

    for (int i = start_idx; i < end_idx && i < num_files; ++i) {
        string input_name = input_files[i];
        fs::path input_path(input_name);
        string base_name = input_path.stem().string();
        size_t pos = base_name.find(".subst");
        if (pos != string::npos) {
            base_name = base_name.substr(0, pos);
        }

        fs::path output_path = fs::path("rna") / (base_name + ".rna");
        string output_name = output_path.string();

        fs::create_directories(output_path.parent_path());

        transcribe_dna_to_rna(input_name, output_name);

        cout << "Process " << rank << " processed file: " << input_name 
             << " and wrote output to: " << output_name << endl;
    }

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    double max_time;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Total execution time: " << max_time << " seconds." << endl;
    }

    MPI_Finalize();
    return 0;
}