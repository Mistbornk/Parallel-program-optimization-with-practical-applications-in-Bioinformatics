#ifndef SW_H
#define SW_H
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <xsimd/xsimd.hpp>
#define MATCH 2
#define MISTMATH 1
#define GAPOPEN 1
#define GAPEXTEND 1

using value_type = int16_t;
std::string read_fasta_sequence(const std::string& filename);

struct SmithWaterman {
    int score;
    std::string aligned_seq1;
    std::string aligned_seq2;
    std::string match_line;
    size_t start1, end1, start2, end2;
};
struct Position {
    value_type score;
    std::size_t start_i, start_j;
    std::size_t end_i, end_j;
};

SmithWaterman striped_smith_waterman(std::string_view ref, std::string_view query,
    int16_t match = MATCH, int16_t mismatch = MISTMATH, int16_t gap_open = GAPOPEN, int16_t gap_extend = GAPEXTEND);

SmithWaterman naive_sw(const std::string& s1, const std::string& s2,
    int match = MATCH, int mismatch = MISTMATH, int gap_open = GAPOPEN, int gap_extend = GAPEXTEND);

#endif
