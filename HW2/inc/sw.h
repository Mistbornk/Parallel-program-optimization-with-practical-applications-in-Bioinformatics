#ifndef SW_H
#define SW_H
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <xsimd/xsimd.hpp>
#define BANDWIDTH 25
#define MATCH 2
#define MISTMATH 2
#define GAPOPEN 3
#define GAPEXTEND 1

std::string read_fasta_sequence(const std::string& filename);

struct SmithWaterman {
    int score;
    std::string aligned_seq1;
    std::string aligned_seq2;
    std::string match_line;
    int start1, end1, start2, end2;
};
struct Position {
    int end1, end2;
};

SmithWaterman banded_smith_waterman(const std::string& s1, const std::string& s2,
    int band_width = BANDWIDTH, int match = MATCH, int mismatch = MISTMATH, int gap_open = GAPOPEN, int gap_extend = GAPEXTEND);

Position  striped_smith_waterman(const std::vector<int16_t> &ref, const std::vector<int16_t> &query,
    int match = MATCH, int mismatch = MISTMATH, int gap_open = GAPOPEN, int gap_extend = GAPEXTEND);

int naive_sw(const std::string& s1, const std::string& s2,
    int match = MATCH, int mismatch = MISTMATH, int gap_open = GAPOPEN, int gap_extend = GAPEXTEND);

SmithWaterman naive_sw_traceback(const std::string& s1, const std::string& s2, int end_i, int end_j,
    int match = MATCH, int mismatch = MISTMATH, int gap_open = GAPOPEN, int gap_extend = GAPEXTEND);

std::vector<int16_t> seq2vec(const std::string &s);

#endif
