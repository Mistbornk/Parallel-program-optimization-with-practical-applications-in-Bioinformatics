#include <iostream>
#include <chrono>
#include "colMajor.hpp"
#include "rowMajor.hpp"
using namespace std;

int main() {
    /* test */
	Column_Major_Matrix<int> cc1 (1000, 1000);
	Row_Major_Matrix<int> rr1( 1000, 1000);
	Column_Major_Matrix<int> cc2 (cc1);
	Row_Major_Matrix<int> rr2 = (rr1);
	Column_Major_Matrix<int> cc3 = std::move( cc2 );
	Row_Major_Matrix<int> rr3 = std::move( rr2 );

	Column_Major_Matrix<int> cc (55, 1000);
	Row_Major_Matrix<int> rr (1000, 66);

	Row_Major_Matrix<int> rr_single = cc * rr;
	
    cout << "Basic test success\n";

    /* test overload % */
	int size = 1000;  // set matrix size
    Row_Major_Matrix<int> A(size, size);
    Column_Major_Matrix<int> B(size, size);

    /* Row major x  Column major */
    cout << "------Row major x Column major------\n";
    auto start = chrono::high_resolution_clock::now();
    Row_Major_Matrix<int> C1 = A * B;
    auto end = chrono::high_resolution_clock::now();
    double t_single = chrono::duration<double>(end - start).count();  

    start = chrono::high_resolution_clock::now();
    Row_Major_Matrix<int> C2 = A % B;
    end = chrono::high_resolution_clock::now();
    auto t_multi = chrono::duration<double>(end - start).count();

    /* Outcome */
    if (C1 == C2) cout << "True" << endl;
    else cout << "False" <<endl;
    cout << "original time: " << t_single << " seconds\n";
    cout << "Multi-thread time: " << t_multi << " seconds\n";
    cout << "Speedup: " << t_single / t_multi << "x\n";

    /* Column major x Row major */
    cout << "------Column major x Row major------\n";
    start = chrono::high_resolution_clock::now();
    Column_Major_Matrix<int> D1 = B * A;
    end = chrono::high_resolution_clock::now();
    t_single = chrono::duration<double>(end - start).count();  

    start = chrono::high_resolution_clock::now();
    Column_Major_Matrix<int> D2 = B % A;
    end = chrono::high_resolution_clock::now();
    t_multi = chrono::duration<double>(end - start).count();

    /* Outcome */
    if (D1 == D2) cout << "True" << endl;
    else cout << "False" <<endl;
    cout << "original time: " << t_single << " seconds\n";
    cout << "Multi-thread time: " << t_multi << " seconds\n";
    cout << "Speedup: " << t_single / t_multi << "x\n";
    return 0;
}