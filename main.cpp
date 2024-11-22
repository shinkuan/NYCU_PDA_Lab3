#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <omp.h>
#include "legalizer.h"

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <lg_file> <opt_file> <output_file>" << std::endl;
        return 1;
    }
    omp_set_num_threads(8); // Limit OpenMP to 8 threads

    auto start = std::chrono::high_resolution_clock::now();

    Legalizer legalizer;
    legalizer.load_lg(argv[1]);
    legalizer.load_opt(argv[2]);
    legalizer.legalize(argv[3]);
    // legalizer.dump_loaded_row("final_row.txt");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Elapsed time: " << elapsed.count() << "s" << std::endl;

    return 0;
}