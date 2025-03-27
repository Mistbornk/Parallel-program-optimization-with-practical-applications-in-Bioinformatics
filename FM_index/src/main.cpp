#include "FM_Index.hpp"
#include <iostream>
using namespace std;

int main() {
    string T = "ACATNCCGTCATGGATTACGTACAG$"; // 注意結尾要加 $
    FMIndex fm(T);

    string pattern = "TTA";

    auto result = fm.query(pattern);
	cout << endl;
	if (result.empty()) cout << "No results\n";
	else  {
		cout << "Pattern \"" << pattern << "\" found at positions: ";
		for (int pos : result) cout << pos << " ";
		cout << "\nMatched substrings: \n";
		for (int pos : result) {
			int n = T.size();
			cout << (T.substr(pos-1, n)) <<endl;
		}
	}
    return 0;
}
