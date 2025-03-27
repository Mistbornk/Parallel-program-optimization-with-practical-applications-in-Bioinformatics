#ifndef FM_INDEX_H
#define FM_INDEX_H
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>

class FMIndex {
public:
	FMIndex(){ };
	FMIndex(const std::string& T) {
		buildBWT(T);
		buildC();
		buildOcc();
		print();
	}
	~FMIndex(){ };
	void buildBWT(const std::string& T);
	void buildC();
	void buildOcc();
	void print();
	std::vector<int> query(const std::string& pattern);

private:
	std::string bwt;
	std::vector<int> suffix_array;
	std::unordered_map<char, int> C;
	std::unordered_map<char, std::vector<int>> Occ;
};

#endif
