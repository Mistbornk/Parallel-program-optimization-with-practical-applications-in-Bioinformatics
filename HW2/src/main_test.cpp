#include "sw.h"
#include <iomanip>
#include <iostream>
#include <chrono>
using namespace std;

void PrintOutcome(SmithWaterman& result, chrono::duration<double> duration, bool use_simd) {
    cout << "Mode: " << (use_simd ? "XSIMD" : "Baseline") << endl;
    cout << "Seq1:  " << setw(5) << result.start1 << "   " << result.aligned_seq1 << "   " << result.end1 << "\n";
    cout << "       " << "        " << result.match_line << "\n";
    cout << "Seq2:  " << setw(5) << result.start2 << "   " << result.aligned_seq2 << "   " << result.end2 << "\n";

    cout << fixed << setprecision(6);
    cout << "Execution time: " << duration.count() << " seconds\n";

    int alignment_length = result.aligned_seq1.size();
    int match = count(result.match_line.begin(), result.match_line.end(), '|');
    float match_ratio = match / float(result.match_line.size());
    int mismatch = count(result.match_line.begin(), result.match_line.end(), '*');
    int gaps = count(result.aligned_seq1.begin(), result.aligned_seq1.end(), '-') +
               count(result.aligned_seq2.begin(), result.aligned_seq2.end(), '-');
    float quality = match_ratio - 0.1 * (mismatch + gaps) / alignment_length;

    cout << fixed << setprecision(4);
    cout << "\n[Alignment Quality Report]\n";
    cout << "Length     : " << alignment_length << "\n";
    cout << "Matches    : " << match << "\n";
    cout << "Match Ratio: " << match_ratio * 100 << "%\n";
    cout << "Mismatches : " << mismatch << "\n";
    cout << "Gaps       : " << gaps << "\n";
    cout << "Quality    : " << quality * 100<< "%\n\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " seq1.fa seq2.fa" << endl;
        return 1;
    }

    string seq1 = read_fasta_sequence(argv[1]);
    string seq2 = read_fasta_sequence(argv[2]);

    SmithWaterman result;
    // banded smith waterman
    cout << "==================================" << endl;
    auto start = chrono::high_resolution_clock::now();
    result = banded_smith_waterman(seq1, seq2, 25);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration_base = end - start;
    PrintOutcome(result, duration_base, false);
    // striped smith waterman 
    cout << "==================================" << endl;
    start = chrono::high_resolution_clock::now();
    result = xsimd_smith_waterman(seq1, seq2);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration_xsimd = end - start;
    PrintOutcome(result, duration_xsimd, true);
    // speed up
    cout << "==================================" << endl;
    double speedup = duration_base.count() / duration_xsimd.count();
    cout << "Speedup       : " << speedup << "x\n\n";

    return 0;
}
