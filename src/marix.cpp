#include <iostream>
#include <chrono>
#include "col_major.hpp"
#include "row_major.hpp"
using namespace std;

int main() {
	Column_Major_Matrix<int> cc1 (1000, 1000);
	Row_Major_Matrix<int> rr1( 1000, 1000);
	Column_Major_Matrix<int> cc2 (cc1);
	Row_Major_Matrix<int> rr2 = (rr1);
	Column_Major_Matrix<int> cc3 = std::move( cc2 );
	Row_Major_Matrix<int> rr3 = std::move( rr2 );

	Column_Major_Matrix<int> cc (55, 1000);
	Row_Major_Matrix<int> rr (1000, 66);

	Row_Major_Matrix<int> rr_single = cc * rr;
	
	int size = 1000;  // 設定矩陣大小
    Row_Major_Matrix<int> A(size, size);
    Column_Major_Matrix<int> B(size, size);
	Row_Major_Matrix<int> D(size, size);

    // **測試單執行緒版本**
    auto start = std::chrono::high_resolution_clock::now();
    Row_Major_Matrix<int> C1 = A * B;
    auto end = std::chrono::high_resolution_clock::now();
    double t_single = std::chrono::duration<double>(end - start).count();  

    // **測試多執行緒版本**
    start = std::chrono::high_resolution_clock::now();
    Row_Major_Matrix<int> C2 = A % B;
    end = std::chrono::high_resolution_clock::now();
    double t_multi = std::chrono::duration<double>(end - start).count();

    // **輸出結果**
    std::cout << "Single-thread time: " << t_single << " seconds\n";
    std::cout << "Multi-thread time: " << t_multi << " seconds\n";
    std::cout << "Speedup: " << t_single / t_multi << "x\n";

}