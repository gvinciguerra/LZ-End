#include <iostream>
#include <chrono>
#include "lzend_parser.hpp"

int main() {
    // Print the LZ-End parsing of a string
    std::cout << lzend_parse("alabar_a_la_alabarda") << std::endl;

    // Print the LZ-End parsing of an integer vector
    std::cout << lzend_parse(sdsl::int_vector<16>{500, 1, 2, 1, 1, 1, 2, 2, 1}) << std::endl;


    // Parse a text file
    sdsl::int_vector<8> v;
    sdsl::load_vector_from_file(v, "proteins.100MB"); // download from http://pizzachili.dcc.uchile.cl/texts/protein/

    std::cout << "Input length:  " << v.size() << std::endl;

    auto t0 = std::chrono::high_resolution_clock::now();
    auto parse = lzend_parse(v); // shorthand for:  LZEnd<unsigned char, LZEndParsingStrategyKN> parse(v);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Parsing time:  " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms"
              << std::endl << "Phrases count: " << parse.size() << std::endl;

    auto decoded = parse.decode();
    for (size_t i = 0; i < v.size(); ++i) {
        if (decoded[i] != v[i]) {
            std::cout << "Error at " << i << std::endl;
            return 1;
        }
    }
    std::cout << "Decompression test passed" << std::endl;

    return 0;
}
