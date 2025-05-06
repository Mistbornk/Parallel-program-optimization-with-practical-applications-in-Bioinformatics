#include "sw.hpp"
#include <iomanip>
#include <iostream>
#include <chrono>
#include <random>
#include <unistd.h>
#define TRY 100
using namespace std;
using namespace chrono;
extern "C" void cuda_warmup();

void striped_warmup() {
    std::string dummy1(128, 'A');
    std::string dummy2(128, 'A');
	striped_smith_waterman(dummy1, dummy2);
}

void randomTest() {
    std::cout << "==================================" << std::endl;
    std::cout << "Random test comparison ("<<TRY<<" runs)" << std::endl;

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> len_dist(1000, 10000);
    std::uniform_int_distribution<int> base_dist(0, 3);
    const char bases[] = {'A', 'C', 'G', 'T'};

    int correct = 0;
    for (int test = 0; test < TRY; ++test) {
        int len1 = len_dist(rng);
        int len2 = len_dist(rng);
        if (len1 < len2) swap(len1, len2);
        std::string r, q;
        for (int i = 0; i < len1; ++i) r += bases[base_dist(rng)];
        for (int i = 0; i < len2; ++i) q += bases[base_dist(rng)];

        auto result1 = striped_smith_waterman(r, q);
        auto result2 = cuda_smith_waterman(r, q);
        if (result1.score == result2.score) ++correct;
        //tuple t1 = {result1.start1, result1.start2, result1.end1, result1.end2};
        //tuple t2 = {result2.start1, result2.start2, result2.end1, result2.end2};
        //if (t1 == t2) ++correct;
    }

    std::cout << "Matching scores: " << correct << "/"<<TRY<<" = "
              << std::fixed << std::setprecision(2)
              << (correct / (double)TRY) * 100 << "%\n";
}

void PrintOutcome(SmithWaterman& result, string method, double time_ms) {
    bool enable_color = isatty(fileno(stdout));
    auto green = [](char c) { return "\033[32m" + std::string(1, c) + "\033[0m"; };
    auto red   = [](char c) { return "\033[31m" + std::string(1, c) + "\033[0m"; };

    std::string colored_seq1, colored_seq2, colored_match;
    for (size_t i = 0; i < result.match_line.size(); ++i) {
        char a = result.aligned_seq1[i];
        char b = result.aligned_seq2[i];
        char m = result.match_line[i];

        if (m == '|') {
            colored_seq1 += enable_color ? green(a) : std::string(1, a);
            colored_seq2 += enable_color ? green(b) : std::string(1, b);
            colored_match += enable_color ? green('|') : "|";
        } else if (m == '*') {
            colored_seq1 += enable_color ? red(a) : std::string(1, a);
            colored_seq2 += enable_color ? red(b) : std::string(1, b);
            colored_match += enable_color ? red('*') : "*";
        } else {
            colored_seq1 += a;
            colored_seq2 += b;
            colored_match += ' ';
        }
    }

    std::cout << "Method: " << method << "\n\n";
    if (method != "naive") {
        std::cout << "Seq1:  " << std::setw(5) << result.start1 << "   " << colored_seq1 << "   " << result.end1 << "\n";
        std::cout << "       " << "        " << colored_match << "\n";
        std::cout << "Seq2:  " << std::setw(5) << result.start2 << "   " << colored_seq2 << "   " << result.end2 << "\n";
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Alignment Score: " << result.score << "\n";
    //std::cout << "End at (i, j): " << "(" << result.end1 << ", " << result.end2 << ")" << "\n";
    std::cout << "Time: " << time_ms << " ms" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " seq1.fa seq2.fa" << endl;
        return 1;
    }

    string seq1 = read_fasta_sequence(argv[1]);
    string seq2 = read_fasta_sequence(argv[2]);

    SmithWaterman result1, result2, result3, result4;
    //double time_naive = 0.0, time_banded = 0.0;
    double  time_striped = 0.0, time_cuda = 0.0;

    // === Ground truth === //
    //cout << "==================================" << endl;
    //auto start1 = high_resolution_clock::now();
    //result1 = naive_sw(seq1, seq2);
    //auto end1 = high_resolution_clock::now();
    //time_naive = duration_cast<duration<double, milli>>(end1 - start1).count();
    //PrintOutcome(result1, "Naive", time_naive);
 
    // === striped smith waterman === //
    cout << "==================================" << endl;
    striped_warmup();
    auto start2 = high_resolution_clock::now();
    result2 = striped_smith_waterman(seq1, seq2);
    auto end2 = high_resolution_clock::now();
    time_striped = duration_cast<duration<double, milli>>(end2 - start2).count();
    PrintOutcome(result2, "Striped", time_striped);

    
    // === cuda smith waterman === //
    cout << "==================================" << endl;
    cuda_warmup();
    auto start3 = high_resolution_clock::now();
    result3 = cuda_smith_waterman(seq1, seq2);
    auto end3 = high_resolution_clock::now();
    time_cuda = duration_cast<duration<double, milli>>(end3 - start3).count();
    PrintOutcome(result3, "CUDA", time_cuda);

    // === banded smith waterman === //
    //cout << "==================================" << endl;
    //auto start4 = high_resolution_clock::now();
    //result4 = banded_smith_waterman(seq1, seq2);
    //auto end4 = high_resolution_clock::now();
    //time_banded = duration_cast<duration<double, milli>>(end4 - start4).count();
    //PrintOutcome(result4, "Banded", time_banded);

    // === Speedup === //
    cout << "==================================" << endl;
    cout << std::fixed << std::setprecision(2);
    cout << "Speedup (Striped SW / Cuda SW): " << time_striped / time_cuda  << "x" << endl;

    //randomTest();
    return 0;
}
