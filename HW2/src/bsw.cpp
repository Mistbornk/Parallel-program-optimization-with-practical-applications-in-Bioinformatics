#include "sw.h"

std::string read_fasta_sequence(const std::string& filename) {
    std::ifstream file(filename);
    std::string line, sequence;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        sequence += line;
    }
    return sequence;
}

SmithWaterman banded_smith_waterman (const std::string& s1, const std::string& s2,
    int band_width, int match, int mismatch, int gap_open, int gap_extend) {
	// initialize dp matrix
	int len1 = s1.size();
    int len2 = s2.size();
	// Affine gap matrix, H:High score matrix, E:Extend gap vertically, F:Extend gap horizontally
    std::vector<std::vector<int>> H(len1 + 1, std::vector<int>(len2 + 1, 0));
    std::vector<std::vector<int>> E(len1 + 1, std::vector<int>(len2 + 1, 0));
    std::vector<std::vector<int>> F(len1 + 1, std::vector<int>(len2 + 1, 0));

	// initialize score and matrix position
	int max_score = 0, max_i = 0, max_j = 0;
	// dp score board
	for(int i=1; i<=len1; i++) {
		int j_start = std::max(1, i - band_width);
		int j_end = std::min(len2, i + band_width);
		for (int j = j_start; j<=j_end; j++) {
			E[i][j] = std::max(H[i-1][j] - gap_open, E[i-1][j] - gap_extend);
			F[i][j] = std::max(H[i][j-1] - gap_open, F[i][j-1] - gap_extend);
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
	int i = max_i, j = max_j;
	while (i > 0 && j > 0 && H[i][j] > 0) {
		if (H[i][j] == H[i-1][j-1] + (s1[i-1] == s2[j-1] ? match : -mismatch)) {
			align1 += s1[i-1];
			align2 += s2[j-1];
			matchline += (s1[i-1] == s2[j-1] ? '|' : '*');
			i--, j--;
		}else if (H[i][j] == E[i][j]) {
			align1 += s1[i-1];
			align2 += '-';
			matchline += ' ';
			i--;
		}else {
			align1 += '-';
			align2 += s2[j-1];
			matchline += ' ';
			j--;
		}
	}
	std::reverse(align1.begin(), align1.end());
	std::reverse(align2.begin(), align2.end());
	std::reverse(matchline.begin(), matchline.end());

	return {max_score, align1, align2, matchline, i + 1, max_i, j + 1, max_j};
}


int naive_sw(const std::string& s1, const std::string& s2,
    int match, int mismatch, int gap_open, int gap_extend) {
    int len1 = s1.size();
    int len2 = s2.size();

    std::vector<std::vector<int>> H(len1 + 1, std::vector<int>(len2 + 1, 0));
    std::vector<std::vector<int>> E(len1 + 1, std::vector<int>(len2 + 1, 0));
    std::vector<std::vector<int>> F(len1 + 1, std::vector<int>(len2 + 1, 0));
    int max_score = 0;

    for (int i = 1; i <= len1; ++i) {
        for (int j = 1; j <= len2; ++j) {
			E[i][j] = std::max(H[i-1][j] - gap_open, E[i-1][j] - gap_extend);
			F[i][j] = std::max(H[i][j-1] - gap_open, F[i][j-1] - gap_extend);
			int h = H[i-1][j-1] + (s1[i-1] == s2[j-1] ? match : -mismatch);
			H[i][j] = std::max(0, std::max({h, E[i][j], F[i][j]}));
            if (H[i][j] > max_score) {
                max_score = H[i][j];
            }
        }
    }
    return max_score;
}