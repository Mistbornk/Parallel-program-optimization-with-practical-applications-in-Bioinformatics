# HW1 Overview
## Github Repository
https://github.com/Mistbornk/Parallel-program-optimization-with-practical-applications-in-Bioinformatics.git
## Implements:
1. **Matrix Operations**:
   - `Row_Major_Matrix`: Row-major storage.
   - `Column_Major_Matrix`: Column-major storage.
   - Supports **matrix multiplication (`*`)** and **multi-threaded multiplication (`%`)**.
2. **Thread Pool**:
   - A **thread pool** that manages **5 threads** and processes jobs asynchronously.

---

## Structure:
```plaintext
.
├── src/            	   	# Source code directory
│   ├── threadPool.cpp     	# Thread Pool implementation
│   ├── matrix_test.cpp    	# Matrix multiplication test executable
│   ├── threadpool_test.cpp # Thread Pool test executable
├── inc/           	# Header files
│   ├── rowMajor.hpp       	# Row-Major matrix class and implementation
│   ├── colMajor.hpp       	# Column-Major matrix class and implementation
│   ├── threadPool.hpp    	# Thread Pool class
├── Makefile                # Build script
├── README.md               # documentation
```

## Usage
### **Compile the HW**  
To build all executables, run:  
```
make
```
This will compile:  
- `./matrix_test`  
- `./threadpool_test`  
---

### **Clean build files**  
To remove compiled object files (`.o`) and executables, run:  
```
make clean
```
---

### **Matrix Multiplication Test**  
To execute matrix multiplication tests:
```
./matrix_test
```
This will first test:
``` C++
Column_Major_Matrix<int> cc1 (1000, 1000);
Row_Major_Matrix<int> rr1( 1000, 1000);
Column_Major_Matrix<int> cc2 (cc1);
Row_Major_Matrix<int> rr2 = (rr1);
Column_Major_Matrix<int> cc3 = std::move( cc2 );
Row_Major_Matrix<int> rr3 = std::move( rr2 );

Column_Major_Matrix<int> cc (55, 1000);
Row_Major_Matrix<int> rr (1000, 66);

Row_Major_Matrix<int> rr_single = cc * rr;
```
And than test the Overload of (`*`) and (`%`):
```C++
int size = 1000;
Row_Major_Matrix<int> A(size, size);
Column_Major_Matrix<int> B(size, size);
Row_Major_Matrix<int> C1 = A * B;
Row_Major_Matrix<int> C2 = A % B;
```
Output: 
```
Basic test success
------Row major x Column major------
True/False: Verify that `C1` (original result) and `C2` (multi-threaded result) are identical. 
original time (seconds)
Multi-thread time (seconds)
Speedup
------Column major x Row major------
True/False: Verify that `C1` (original result) and `C2` (multi-threaded result) are identical. 
original time (seconds)
Multi-thread time (seconds)
Speedup
```
---

### **Thread Pool Test**  
To execute the thread pool test:  
```
./threadpool_test
```
This will test:
```
First send 496 print_1 functions and then 4 Print_2 functors into the pool
for (int i = 0; i < 496; ++i) {
	pool.enqueueJobs(print_1);
}
for (int i = 0; i < 4; ++i) {
	pool.enqueueJobs(Print2());
}
```
Output:
- print_1 results
- Print_2 results
- [ThreadPool Summary] follow by 5 Thread INFO
	- Thread ID
	- The total running time
	- The total life time
