# HW1 Overview
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
├── inc/           			# Header files
│   ├── rowMajor.hpp       	# Row-Major matrix class and implementation
│   ├── colMajor.hpp       	# Column-Major matrix class and implementation
│   ├── threadPool.hpp    	# Thread Pool class
├── Makefile           # Build script
├── README.md          # Project documentation
```

## Usage
### ** Compile the HW **  
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
---

### **Thread Pool Test**  
To execute the thread pool test:  
```
./threadpool_test
```