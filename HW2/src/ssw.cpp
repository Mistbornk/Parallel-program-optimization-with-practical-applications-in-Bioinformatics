#include "sw.h"

std::vector<simd_t> build_query_profile(const std::vector<uint8_t>& read, int match, int mismatch, uint8_t bias) {
	const std::vector<uint8_t> bases = {'A', 'C', 'G', 'T'};
	const int bases_size = bases.size();
	int readLen = read.size();
	int segLen = (readLen + 15) / 16;
	std::vector<simd_t> profile(bases_size * segLen);

    for (int b = 0; b < bases_size; b++) {
        uint8_t base = bases[b];
        for (int i = 0; i < segLen; i++) {
            alignas(16) uint8_t buf[16] = {};
            for (int seg = 0; seg < 16; seg++) {
                int j = i + seg * segLen;
                buf[seg] = j >= readLen ? bias
                    : ((base == read[j]) ? match : mismatch) + bias;
            }
            profile[b * segLen + i] = simd_t::load_aligned(buf);
        }
    }
    return profile;
}

SmithWaterman xsimd_smith_waterman(const std::string& s1, const std::string& s2,
    int match, int mismatch, int gap_open, int gap_extend) {
    int m = s1.size(), n2 = s2.size();
    std::vector<uint8_t> query(s1.begin(), s1.end());
    std::vector<uint8_t> ref(s2.begin(), s2.end());

    uint8_t bias = std::abs(mismatch);
    int segLen = (m + 15) / 16;
    auto profile = build_query_profile(query, match, mismatch, bias);

    std::vector<simd_t> H(segLen, simd_t(0)), E(segLen, simd_t(0)), Hmax(segLen, simd_t(0));
    simd_t vGapO(gap_open), vGapE(gap_extend), vBias(bias);

    uint8_t max_score = 0;
    int max_i = 0, max_j = 0;

    for (int i = 0; i < n2; ++i) {
        simd_t F(0), vMaxColumn(0);
        simd_t vH = H[segLen - 1];
        vH = xsimd::slide_left<1>(vH);

        int base_idx = 0;
        switch (ref[i]) {
            case 'A': base_idx = 0; break;
            case 'C': base_idx = 1; break;
            case 'G': base_idx = 2; break;
            case 'T': base_idx = 3; break;
            default: base_idx = 0; break;
        }

        auto vP = &profile[base_idx * segLen];
        std::swap(H, Hmax);

        for (int j = 0; j < segLen; ++j) {
            vH = saturated_add(vH, vP[j]);
            vH = saturated_sub(vH, vBias);

            simd_t e = E[j];
            vH = xsimd::max(vH, e);
            vH = xsimd::max(vH, F);

            vMaxColumn = xsimd::max(vMaxColumn, vH);
            H[j] = vH;

            vH = saturated_sub(vH, vGapO);
            e = saturated_sub(e, vGapE);
            E[j] = xsimd::max(e, vH);

            F = saturated_sub(F, vGapE);
            F = xsimd::max(F, vH);

            vH = Hmax[j];
        }

        auto reduced = vMaxColumn;
        reduced = xsimd::max(reduced, xsimd::slide_right<8>(reduced));
        reduced = xsimd::max(reduced, xsimd::slide_right<4>(reduced));
        reduced = xsimd::max(reduced, xsimd::slide_right<2>(reduced));
        reduced = xsimd::max(reduced, xsimd::slide_right<1>(reduced));
        uint8_t curr_max = reduced.get(0);

        if (curr_max > max_score) {
            max_score = curr_max;
            max_i = segLen - 1; // approximate last row
            max_j = i;
        }
    }

    std::string align1, align2, matchline;
    int i = max_i * 16, j = max_j; // approximate index in original sequence

    while (i > 0 && j > 0 && s1[i - 1] && s2[j - 1]) {
        if (s1[i - 1] == s2[j - 1]) {
            align1 += s1[i - 1];
            align2 += s2[j - 1];
            matchline += '|';
            --i; --j;
        } else {
            align1 += s1[i - 1];
            align2 += s2[j - 1];
            matchline += '*';
            --i; --j;
        }
    }

    std::reverse(align1.begin(), align1.end());
    std::reverse(align2.begin(), align2.end());
    std::reverse(matchline.begin(), matchline.end());

    return {static_cast<int>(max_score), align1, align2, matchline, i + 1, max_i * 16, j + 1, max_j};
}