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
    cout << "Alignment Score: " << result.score << endl;
    cout << "Execution time: " << duration.count() << " seconds\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " seq1.fa seq2.fa" << endl;
        return 1;
    }

    string seq1 = read_fasta_sequence(argv[1]);
    string seq2 = read_fasta_sequence(argv[2]);

    SmithWaterman result;
    // Ground truth
    cout << "==================================" << endl;
    cout << "Match score truth: " << naive_sw(seq1, seq2) << endl;
    // banded smith waterman
    cout << "==================================" << endl;
    auto start = chrono::high_resolution_clock::now();
    result = banded_smith_waterman(seq1, seq2);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration_base = end - start;
    PrintOutcome(result, duration_base, false);

    // striped smith waterman 
    cout << "==================================" << endl;
    start = chrono::high_resolution_clock::now();
    Position  pos = striped_smith_waterman(seq2vec(seq1), seq2vec(seq2));
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration_xsimd = end - start;
    result = naive_sw_traceback(seq1, seq2, pos.end1, pos.end2);
    PrintOutcome(result, duration_xsimd, true);
    // speed up
    cout << "==================================" << endl;
    double speedup = duration_base.count() / duration_xsimd.count();
    cout << "Speedup       : " << speedup << "x\n\n";

    return 0;
}
