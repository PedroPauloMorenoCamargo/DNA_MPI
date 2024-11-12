#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cctype>
#include <string>
#include <mpi.h>
#include <omp.h>
#include <filesystem>

namespace fs = std::filesystem;
using namespace std;

void preprocess(const string& input_name, const string& output_name, const array<char, 256>& table, const size_t& buffer_size) {
    // Abre as files no modo binário
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

    //Lê a primeira linha do arquivo
    string first_line;
    getline(input, first_line);

    //Cria buffer de input e output
    vector<char> input_buffer(buffer_size);
    vector<char> output_buffer;
    output_buffer.reserve(buffer_size);

    //Lê o buffer
    while (input.read(input_buffer.data(), buffer_size) || input.gcount()) {
        // Numero de Bytes
        streamsize bytes_read = input.gcount();
        output_buffer.clear();

        // Preenche o Buffer
        for (streamsize i = 0; i < bytes_read; ++i) {
            unsigned char ch = static_cast<unsigned char>(input_buffer[i]);
            char mapped_char = table[ch];
            if (mapped_char) {
                output_buffer.push_back(mapped_char);
            }
        }
        // Escreve para o output
        if (!output_buffer.empty()) {
            output.write(output_buffer.data(), output_buffer.size());
            if (!output) {
                cerr << "Error writing to output file." << endl;
                return;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    // Inicializa MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Começa o Tempo
    double program_start_time = MPI_Wtime();

    // Tabela de LookUP
    array<char, 256> table = {};
    for (int i = 0; i < 256; ++i) {
        table[i] = 0;
    }
    // Nucleotideos
    const string nucleotides = "ATCG";
    for (char ch : nucleotides) {
        unsigned char uc = static_cast<unsigned char>(ch);
        table[uc] = tolower(uc);
        table[tolower(uc)] = tolower(uc);
    }

    //Define o Tamanho do Buffer
    const size_t buffer_size = 8192;

    // Lista dos inputs
    vector<string> input_files = {
        "entradas/chr1.subst.fa", "entradas/chr2.subst.fa", "entradas/chr3.subst.fa",
        "entradas/chr4.subst.fa", "entradas/chr5.subst.fa", "entradas/chr6.subst.fa",
        "entradas/chr7.subst.fa", "entradas/chr8.subst.fa", "entradas/chr9.subst.fa",
        "entradas/chr10.subst.fa", "entradas/chr11.subst.fa", "entradas/chr12.subst.fa",
        "entradas/chr13.subst.fa", "entradas/chr14.subst.fa", "entradas/chr15.subst.fa",
        "entradas/chr16.subst.fa", "entradas/chr17.subst.fa", "entradas/chr18.subst.fa",
        "entradas/chr19.subst.fa", "entradas/chr20.subst.fa", "entradas/chr21.subst.fa",
        "entradas/chr22.subst.fa"
    };

    // Cada processo lê arquivos diferentes
    int files_per_proc = input_files.size() / size;
    int start_idx = rank * files_per_proc;
    int end_idx = (rank == size - 1) ? input_files.size() : start_idx + files_per_proc;

    for (int i = start_idx; i < end_idx; ++i) {
        string input_name = input_files[i];

        // Extrai o nome do arquivo
        size_t last_slash = input_name.find_last_of("/\\");
        if (last_slash == string::npos) {
            last_slash = -1;
        }
        size_t dot_pos = input_name.find(".fa", last_slash);
        if (dot_pos == string::npos) {
            cerr << "Invalid input file name: " << input_name << endl;
            continue;
        }
        string base_name = input_name.substr(last_slash + 1, dot_pos - last_slash - 1);
        string output_name = (output_dir / (base_name + ".filtered")).string();

        cout << "Process " << rank << " processing file: " << input_name << " -> " << output_name << endl;
        preprocess(input_name, output_name, table, buffer_size);
    }

    // Termina o Tempo
    double program_end_time = MPI_Wtime();
    double program_elapsed_time = program_end_time - program_start_time;

    // Soma todos os Tempos
    double max_program_time;
    MPI_Reduce(&program_elapsed_time, &max_program_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Total execution time: " << max_program_time << " seconds." << endl;
    }

    // Finaliza MPI
    MPI_Finalize();
    return 0;
}
