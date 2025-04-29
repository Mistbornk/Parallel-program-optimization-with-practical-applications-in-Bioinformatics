#include "sw.h"
#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <xsimd/xsimd.hpp>

constexpr auto kBases = uint8_t{4};
using sub_t = std::array<std::array<int16_t, kBases>, kBases>;
using vector_unaligned = std::vector<value_type>;
using vector_aligned = std::vector<value_type, xsimd::default_allocator<value_type>>;

sub_t substitution_matrix(int16_t m, int16_t x) {
    return sub_t {{
        {m, x, x, x},
        {x, m, x, x},
        {x, x, m, x},
        {x, x, x, m}
    }};
}

uint8_t encode(char c) {
	switch (c) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		default:  return 0; // fallback
	}
}

SmithWaterman striped_smith_waterman(std::string_view ref, std::string_view query, 
                                 int16_t match, int16_t mismatch, int16_t gap_open, int16_t gap_extend) {

    constexpr std::size_t batchSize = xsimd::batch<value_type>::size;
    size_t segLen = (size(query) + batchSize) / batchSize;
    size_t profile_len = segLen * batchSize;
	sub_t sub_mat = substitution_matrix(match, -mismatch);

    std::array<vector_aligned, kBases> profile;
	for (uint8_t base = 0; base < kBases; ++base) {
		profile[base].resize(profile_len, 0);
		for (size_t j = 0; j < segLen; ++j) {
			for (size_t k = (j==0) ? 1 : 0; k < batchSize; ++k) {
				size_t idx = k * segLen + j - 1;
				if (idx < size(query)) {
					uint8_t code = encode(query[idx]);
					profile[base][j * batchSize + k] = sub_mat[base][code];
				}
			}
		}
	}
	//std::cout << "=== Profile Debug ===" << std::endl;
	//for (uint8_t base = 0; base < kBases; ++base) {
	//	std::cout << "Base: " << "ACGT"[base] << std::endl;
	//	for (size_t j = 0; j < segLen; ++j) {
	//		for (size_t k = 0; k < batchSize; ++k) {
	//			size_t idx = j * batchSize + k;
	//			if (idx < profile[base].size()) {
	//				std::cout << std::setw(4) << int(profile[base][idx]) << " ";
	//			}
	//		}
	//		std::cout << std::endl;
	//	}
	//	std::cout << "-------------------------" << std::endl;
	//}
	
    auto H = std::vector(size(ref) + 1, vector_aligned(profile_len));
    auto E = std::vector(size(ref) + 1, vector_aligned(profile_len));
    auto F = std::vector(size(ref) + 1, vector_aligned(profile_len));
	auto GAP_O = xsimd::batch<value_type>(gap_open);
	auto GAP_E = xsimd::batch<value_type>(gap_extend);
	auto ZERO = xsimd::batch<value_type>{0};

    size_t max_i  = 0, max_j = 0, max_k = 0;
    value_type max_score = 0;

    for (size_t i = 1; i <= size(ref); ++i) {
        auto vH = std::array{
            xsimd::slide_left<sizeof(value_type)>(
                xsimd::load_aligned(&H[i - 1][profile_len - batchSize])),
            ZERO
        };
        auto vF = ZERO;
        auto base = encode(ref[i - 1]);
        for (size_t j = 0; j < segLen; ++j) {

			bool even = (j & 1) == 0;
			auto& vH_prev = vH[even];	// (j%2==0) ? vH[1] : vH[0] -> i-1, j-1 H buffer
			auto& vH_curr = vH[!even];	// (j%2==0) ? vH[0] : vH[1] -> i-1, j H buffer
			
			// 更新 vF（vertical gap extension）
			vF = xsimd::max(vF - GAP_E, vH_prev - GAP_O);
			
			// 載入 H[i-1][j]
			vH_prev = xsimd::load_aligned(&H[i - 1][j * batchSize]);
			
			// 更新 vE（horizontal gap extension）
			auto vE = xsimd::max(
				xsimd::load_aligned(&E[i - 1][j * batchSize]) - GAP_E,
				vH_prev - GAP_O
			);
			
			// 更新 vH（diagonal match/mismatch or extend from E/F）
			vH_curr = xsimd::max(
				max(
					vH_curr + xsimd::load_aligned(&profile[base][j * batchSize]),
					ZERO
				),
				max(vE, vF)
			);
			
			// 更新目前 max_score
			if (value_type new_score = xsimd::reduce_max(vH_curr); new_score > max_score) {
				max_i = i;
				max_j = j;
				max_score = new_score;
			}
			
			// 存回 H、E、F
			xsimd::store_aligned(&H[i][j * batchSize], vH_curr);
			xsimd::store_aligned(&E[i][j * batchSize], vE);
			xsimd::store_aligned(&F[i][j * batchSize], vF);			
		}
		// ========== Lazy F Loop ========== //
		size_t last_idx = ~segLen & 1; // if (segLen % 2 ==0) ? vH[1] : vH[0] 從 vH buffer 取得前一個 stripe 最後的資料
        vF = xsimd::slide_left<sizeof(value_type)>(
			xsimd::max(
				vF - GAP_E,
				vH[last_idx] - GAP_O
        ));

		size_t j = 0, pass = 0;
        while (true) {
            auto vh = xsimd::load_aligned(&H[i][j * batchSize]);
            if (xsimd::count(vF > vh - GAP_O) == 0) break;
            xsimd::store_aligned(&H[i][j * batchSize], max(vh, vF));
            xsimd::store_aligned(&F[i][j * batchSize],
                    max(xsimd::load_aligned(&F[i][j * batchSize]), vF));
            if (j + 1 == segLen) {
                vF = xsimd::slide_left<sizeof(value_type)>(vF - GAP_E);
                j = 0;
				if (++pass > 2) break; 
            } else {
                vF -= GAP_E;
				j++;
            }
        }
    }

    // ========== Traceback ========== //
	for (size_t k = 0; k < batchSize; ++k) {
		if (H[max_i][max_j * batchSize + k] == max_score) {
			max_k = k;
			break;
		}
	}
	size_t i = max_i, j = max_j, k = max_k;
	value_type best_score = 0;
	// if j==0: 跳回上一列最後一個 stripe
	auto move_left = [&]() {
        if (j == 0) {
            --k; // slide right
            j = segLen - 1;
        } else {
            --j;
        }
    };

	std::string align1, align2, matchline;
	while (H[i][j * batchSize + k] > 0) {
		auto q_pos = k * segLen + j - 1;
		//auto score = (ref[i - 1] == query[q_pos]) ? match : -mismatch;
		if (H[i][j * batchSize + k] == E[i][j * batchSize + k]) { // horizonal gap
			align1 += ref[i - 1];
            align2 += '-';
            matchline += ' ';

            while (E[i][j * batchSize + k] == E[i - 1][j * batchSize + k] - gap_extend) { // if gap extend from prev row
                --i;
				align1 += ref[i - 1];
				align2 += '-';
				matchline += ' ';
            }
			--i;
		}
		else if (H[i][j * batchSize + k] == F[i][j * batchSize + k]) {
			align1 += '-';
            align2 += query[q_pos];
            matchline += ' ';

			auto prev_idx = (j == 0) ? profile_len - batchSize + k - 1 : (j - 1) * batchSize + k;
			while (F[i][j] == F[i][prev_idx] - gap_extend) {
				move_left();
				prev_idx = (j == 0) ? profile_len - batchSize + k - 1 : (j - 1) * batchSize + k;

				align1 += '-';
				align2 += query[q_pos];
				matchline += ' ';
			}
			move_left();
		}
		else { // if (H[i][j * batchSize + k] == H[i - 1][prev_idx] + score)
			align1 += ref[i - 1];
			align2 += query[q_pos];
			matchline += (ref[i - 1] == query[q_pos]) ? '|' : '*';
			--i;
			move_left();
		}
	}

    std::reverse(align1.begin(), align1.end());
    std::reverse(align2.begin(), align2.end());
    std::reverse(matchline.begin(), matchline.end());

    return SmithWaterman{
        .score = max_score,
        .aligned_seq1 = std::move(align1),
        .aligned_seq2 = std::move(align2),
        .match_line = std::move(matchline),
        .start1 = i,
        .end1 = max_i,
        .start2 = k * segLen + j,
        .end2 = max_k * segLen + max_j
    };
}
