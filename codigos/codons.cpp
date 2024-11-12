#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>
#include <filesystem>
namespace fs = std::filesystem;
using namespace std;
// Tabela de char para bits
unsigned char char_to_bits[256];

void init_char_to_bits() {
    // Mapeia char para bits
    for (int i = 0; i < 256; ++i) {
        char_to_bits[i] = 0xFF; 
    }
    // 00
    char_to_bits['a'] = 0;
    // 01 
    char_to_bits['c'] = 1;
    // 10
    char_to_bits['g'] = 2;
    // 11
    char_to_bits['u'] = 3;
}

int codon_to_index(const char* codon) {
    unsigned char b0 = char_to_bits[static_cast<unsigned char>(codon[0])];
    unsigned char b1 = char_to_bits[static_cast<unsigned char>(codon[1])];
    unsigned char b2 = char_to_bits[static_cast<unsigned char>(codon[2])];
    // Checa caso de caracter invalido
    if (b0 == 0xFF || b1 == 0xFF || b2 == 0xFF) {
        return -1;
    }
    // Mapeia bits de 0 a 63
    return (b0 << 4) | (b1 << 2) | b2;
}

void create_codon_table(int amino_acid_table[64]) {
    // Inicializa Tabela
    for (int i = 0; i < 64; ++i) {
        amino_acid_table[i] = 0;
    }

    // Methionine (Start codon)
    amino_acid_table[codon_to_index("aug")] = 1;

    // Phenylalanine
    amino_acid_table[codon_to_index("uuu")] = 2;
    amino_acid_table[codon_to_index("uuc")] = 2;

    // Leucine
    amino_acid_table[codon_to_index("uua")] = 3;
    amino_acid_table[codon_to_index("uug")] = 3;
    amino_acid_table[codon_to_index("cuu")] = 3;
    amino_acid_table[codon_to_index("cuc")] = 3;
    amino_acid_table[codon_to_index("cua")] = 3;
    amino_acid_table[codon_to_index("cug")] = 3;

    // Isoleucine
    amino_acid_table[codon_to_index("auu")] = 4;
    amino_acid_table[codon_to_index("auc")] = 4;
    amino_acid_table[codon_to_index("aua")] = 4;

    // Valine
    amino_acid_table[codon_to_index("guu")] = 5;
    amino_acid_table[codon_to_index("guc")] = 5;
    amino_acid_table[codon_to_index("gua")] = 5;
    amino_acid_table[codon_to_index("gug")] = 5;

    // Serine
    amino_acid_table[codon_to_index("ucu")] = 6;
    amino_acid_table[codon_to_index("ucc")] = 6;
    amino_acid_table[codon_to_index("uca")] = 6;
    amino_acid_table[codon_to_index("ucg")] = 6;
    amino_acid_table[codon_to_index("agu")] = 6;
    amino_acid_table[codon_to_index("agc")] = 6;

    // Proline
    amino_acid_table[codon_to_index("ccu")] = 7;
    amino_acid_table[codon_to_index("ccc")] = 7;
    amino_acid_table[codon_to_index("cca")] = 7;
    amino_acid_table[codon_to_index("ccg")] = 7;

    // Threonine
    amino_acid_table[codon_to_index("acu")] = 8;
    amino_acid_table[codon_to_index("acc")] = 8;
    amino_acid_table[codon_to_index("aca")] = 8;
    amino_acid_table[codon_to_index("acg")] = 8;

    // Alanine
    amino_acid_table[codon_to_index("gcu")] = 9;
    amino_acid_table[codon_to_index("gcc")] = 9;
    amino_acid_table[codon_to_index("gca")] = 9;
    amino_acid_table[codon_to_index("gcg")] = 9;

    // Tyrosine
    amino_acid_table[codon_to_index("uau")] = 10;
    amino_acid_table[codon_to_index("uac")] = 10;

    // Histidine
    amino_acid_table[codon_to_index("cau")] = 11;
    amino_acid_table[codon_to_index("cac")] = 11;

    // Glutamine
    amino_acid_table[codon_to_index("caa")] = 12;
    amino_acid_table[codon_to_index("cag")] = 12;

    // Asparagine
    amino_acid_table[codon_to_index("aau")] = 13;
    amino_acid_table[codon_to_index("aac")] = 13;

    // Lysine
    amino_acid_table[codon_to_index("aaa")] = 14;
    amino_acid_table[codon_to_index("aag")] = 14;

    // Aspartic Acid
    amino_acid_table[codon_to_index("gau")] = 15;
    amino_acid_table[codon_to_index("gac")] = 15;

    // Glutamic Acid
    amino_acid_table[codon_to_index("gaa")] = 16;
    amino_acid_table[codon_to_index("gag")] = 16;

    // Cysteine
    amino_acid_table[codon_to_index("ugu")] = 17;
    amino_acid_table[codon_to_index("ugc")] = 17;

    // Tryptophan
    amino_acid_table[codon_to_index("ugg")] = 18;

    // Arginine
    amino_acid_table[codon_to_index("cgu")] = 19;
    amino_acid_table[codon_to_index("cgc")] = 19;
    amino_acid_table[codon_to_index("cga")] = 19;
    amino_acid_table[codon_to_index("cgg")] = 19;
    amino_acid_table[codon_to_index("aga")] = 19;
    amino_acid_table[codon_to_index("agg")] = 19;

    // Glycine
    amino_acid_table[codon_to_index("ggu")] = 20;
    amino_acid_table[codon_to_index("ggc")] = 20;
    amino_acid_table[codon_to_index("gga")] = 20;
    amino_acid_table[codon_to_index("ggg")] = 20;

    // Stop codons
    amino_acid_table[codon_to_index("uaa")] = -1;
    amino_acid_table[codon_to_index("uag")] = -1;
    amino_acid_table[codon_to_index("uga")] = -1;

}

void translate_rna_to_amino_acids(const std::string& input_name, const std::string& output_name, const int amino_acid_table[64]) {
    // I/O 
    ifstream input(input_name, ios::binary);
    ofstream output(output_name, ios::binary);

    if (!input) {
        std::cerr << "Error opening input file: " << input_name << endl;
        return;
    }
    if (!output) {
        std::cerr << "Error opening output file: " << output_name << endl;
        return;
    }

    // Cria buffer
    const size_t buffer_size = 8192;
    vector<char> buffer(buffer_size);

    string partial_codon;

    while (input.read(buffer.data(), buffer_size) || input.gcount()) {
        size_t bytes_read = input.gcount();

        // Combina Codong parcial com a data
        string data = partial_codon + string(buffer.data(), bytes_read);

        size_t data_length = bytes_read + partial_codon.length();
        size_t num_codons = data_length / 3;


        //Vetor que guarda os amino_acidos
        vector<int> amino_acids(num_codons);

        // Processa codons em paralelo
        #pragma omp parallel for
        for (size_t i = 0; i < num_codons; ++i) {
            size_t idx = i * 3;
            const char* codon_ptr = data.data() + idx;

            int index = codon_to_index(codon_ptr);

            if (index != -1) {
                amino_acids[i] = amino_acid_table[index];
            } 
            else {
                amino_acids[i] = 0;
            }
        }

        //Escreve para o arquivo
        for (size_t i = 0; i < num_codons; ++i) {
            int aa = amino_acids[i];
            if (aa > 0) {
                output << aa << ' ';
            }
        }

        // Tratativa codons parciais
        size_t remaining_chars = data_length % 3;
        partial_codon = data.substr(data_length - remaining_chars, remaining_chars);
    }

    output.close();
    input.close();
}

int main(int argc, char* argv[]) {
    // Inicializa MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Mapeiamento de caracteres
    init_char_to_bits();

    double start_time = MPI_Wtime();

    // Cria tabela de codons
    int amino_acid_table[64];
    create_codon_table(amino_acid_table);

    // Input files
    std::vector<std::string> input_files = {
        "rna/chr1.rna", "rna/chr2.rna", "rna/chr3.rna",
        "rna/chr4.rna", "rna/chr5.rna", "rna/chr6.rna",
        "rna/chr7.rna", "rna/chr8.rna", "rna/chr9.rna",
        "rna/chr10.rna", "rna/chr11.rna", "rna/chr12.rna",
        "rna/chr13.rna", "rna/chr14.rna", "rna/chr15.rna",
        "rna/chr16.rna", "rna/chr17.rna", "rna/chr18.rna",
        "rna/chr19.rna", "rna/chr20.rna", "rna/chr21.rna",
        "rna/chr22.rna"
    };

    // Distribui os arquivos entre os processos
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

    // Roda a função
    for (int i = start_idx; i < end_idx && i < num_files; ++i) {
        string input_name = input_files[i];
        fs::path input_path(input_name);
        string base_name = input_path.stem().string();

        fs::path output_path = fs::path("sintese_proteica") / (base_name + ".protein");
        string output_name = output_path.string();

        fs::create_directories(output_path.parent_path());

        translate_rna_to_amino_acids(input_name, output_name, amino_acid_table);

        cout << "Process " << rank << " processed file: " << input_name << " and wrote output to: " << output_name << endl;
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
