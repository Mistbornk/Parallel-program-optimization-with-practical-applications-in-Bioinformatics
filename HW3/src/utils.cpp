#include "sw.h"
#include <iostream>

std::string read_fasta_sequence(const std::string& filename) {
    std::ifstream file(filename);
    std::string line, sequence;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        sequence += line;
    }
    return sequence;
}

SmithWaterman naive_sw(const std::string& s1, const std::string& s2,
    int match, int mismatch, int gap_open, int gap_extend) {
    int len1 = s1.size();
    int len2 = s2.size();

    std::vector<std::vector<int>> H(len1 + 1, std::vector<int>(len2 + 1, 0));
    std::vector<std::vector<int>> E(len1 + 1, std::vector<int>(len2 + 1, 0));
    std::vector<std::vector<int>> F(len1 + 1, std::vector<int>(len2 + 1, 0));
    size_t max_i = 0, max_j = 0;
	value_type max_score = 0;

    for (int i = 1; i <= len1; ++i) {
        for (int j = 1; j <= len2; ++j) {
			E[i][j] = std::max(H[i-1][j] - gap_open, std::max(0, E[i-1][j] - gap_extend));
			F[i][j] = std::max(H[i][j-1] - gap_open, std::max(0, F[i][j-1] - gap_extend));
			int h = H[i-1][j-1] + (s1[i-1] == s2[j-1] ? match : -mismatch);
			H[i][j] = std::max(0, std::max({h, E[i][j], F[i][j]}));
            if (H[i][j] > max_score) {
                max_score = H[i][j];
				max_i = i;
				max_j = j;
            }
        }
    }
	// trace back
	std::string align1, align2, matchline;
	size_t i = max_i, j = max_j;
    //while (i > 0 && j > 0 && H[i][j] > 0) {
    //    int score = (s1[i - 1] == s2[j - 1]) ? match : mismatch;
    //    // 先判斷是否來自斜對角 match/mismatch
    //    if (H[i][j] == H[i - 1][j - 1] + score) {
    //        align1 += s1[i - 1];
    //        align2 += s2[j - 1];
    //        matchline += (s1[i - 1] == s2[j - 1]) ? '|' : '*';
    //        --i; --j;
    //    }
    //    // gap in s2（s1 有字元，s2 補 -）
    //    else if (H[i][j] == E[i][j]) {
    //        align1 += s1[i - 1];
    //        align2 += '-';
    //        matchline += ' ';
    //        --i;
    //    }
    //    // gap in s1（s2 有字元，s1 補 -）
    //    else  { // if (H[i][j] == F[i][j])
    //        align1 += '-';
    //        align2 += s2[j - 1];
    //        matchline += ' ';
    //        --j;
    //    }
    //}  
	while (i > 0 && j > 0 && H[i][j] > 0) {
		if (H[i][j] == E[i][j]) {
			align1 += s1[i-1];
			align2 += '-';
			matchline += ' ';
			i--;
		}else if (H[i][j] == F[i][j]) {
			align1 += '-';
			align2 += s2[j-1];
			matchline += ' ';
			j--;
		} else {
			align1 += s1[i-1];
			align2 += s2[j-1];
			matchline += (s1[i-1] == s2[j-1] ? '|' : '*');
			i--, j--;
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
        .start2 = j,
        .end2 = max_j
    };
}