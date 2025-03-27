#include "FM_Index.hpp"

void FMIndex::buildBWT(const std::string& T) {
	// check if '$' at the end
	if (T.find('$') == std::string::npos) {
		throw std::invalid_argument("Must have \"$\" at the end. ");
	}
	// store the rotation case for string
	const int n = T.size();
	std::vector<std::string> rotations;
	for (int i=0; i<n; i ++) {
		rotations.push_back(T.substr(i) + T.substr(0, i));
	}
	// sort the rotation list
	std::sort(rotations.begin(), rotations.end());
	// initialize
	bwt = "";
	suffix_array.clear();
	// store the bwt
	for (int i=0; i<n; i++) {
		bwt += rotations[i][n-1];
		suffix_array.push_back((n - rotations[i].find('$')) % n);
	}
}

void FMIndex::buildC() {
	std::string sorted_bwt = bwt;
	sort(sorted_bwt.begin(), sorted_bwt.end());
	
	std::unordered_map<char, int> count;
	for (char c : sorted_bwt) {
		count[c]++;
	}
	int total = 0;
	for (char c : sorted_bwt) {
		if (C.count(c) == 0) {
			C[c] = total;
			total += count[c];
		}
	}
}

void FMIndex::buildOcc() {
	std::unordered_map<char, int> freq;
	for (int i = 0; i < bwt.size(); ++i) {
		for (auto& p : C) {
			char c = p.first;
			(i==0) ? Occ[c].push_back(bwt[i]==c ? 1 : 0)
				   : Occ[c].push_back(Occ[c][i-1] + (bwt[i]==c ? 1 : 0));
		}
	}
}

void FMIndex::print() {
	std::cout << "BWT: " << bwt << std::endl;
	std::cout << "C table:\n";
	for (auto& p : C) {
		std::cout << "  C(" << p.first << ") = " << p.second << " ";
	}
	std::cout << std::endl;

	std::cout << "Suffix Array:\n";
	for (int i = 0; i < suffix_array.size(); ++i) {
		std::cout << "  SA[" << i << "] = " << suffix_array[i] << std::endl;
	}
}

std::vector<int> FMIndex::query(const std::string& pattern) {
	int m = pattern.size();
	int sp = 0, ep = bwt.size() - 1;
	
	for (int i=m-1; i>=0; i--) {
		char c = pattern[i];
		if (C.count(c) == 0) return {}; // not exist
		int c_occ_sp = sp > 0 ? Occ[c][sp - 1] : 0;
		int c_occ_ep = Occ[c][ep];

		sp = C[c] + c_occ_sp;
		ep = C[c] + c_occ_ep - 1;
		if (sp > ep) return {};
	}

	// return suffix array
	std::vector<int> result;
	for (int i=sp; i<=ep; i++) {
		result.push_back(suffix_array[i]);
	}
	return result;
}