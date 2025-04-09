#ifndef SW_H
#define SW_H

#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <xsimd/xsimd.hpp>

using xsimd::batch;
using xsimd::batch_cast;
using simd_t = xsimd::batch<uint8_t, xsimd::default_arch>;

inline simd_t saturated_add(const simd_t& a, const simd_t& b) {
    return xsimd::min(a + b, simd_t(255));
}

inline simd_t saturated_sub(const simd_t& a, const simd_t& b) {
    return xsimd::max(a - b, simd_t(0));
}

std::string read_fasta_sequence(const std::string& filename);

struct SmithWaterman {
    int score;
    std::string aligned_seq1;
    std::string aligned_seq2;
    std::string match_line;
    int start1, end1, start2, end2;
};

SmithWaterman banded_smith_waterman(const std::string& seq1, const std::string& seq2,
    int band_width = 20, int match = 2, int mismatch = -2, int gap_open = -3, int gap_extend = -1);

SmithWaterman xsimd_smith_waterman(const std::string& s1, const std::string& s2,
    int match = 2, int mismatch = -2, int gap_open = -3, int gap_extend = -1);

std::vector<simd_t> build_query_profile(const std::vector<uint8_t>& read,
    int match, int mismatch, uint8_t bias);

#endif
