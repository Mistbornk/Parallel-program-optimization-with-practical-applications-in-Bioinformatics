# Striped Smith-Waterman (XSIMD) and CUDA Alignment

github repository: [https://github.com/Mistbornk/Parallel-program-optimization-with-practical-applications-in-Bioinformatics](https://github.com/Mistbornk/Parallel-program-optimization-with-practical-applications-in-Bioinformatics)

This project implements the Smith-Waterman local alignment algorithm in two modes:

* **XSIMD-accelerated version** using `xsimd` for vectorized CPU performance
* **CUDA version** using `cuda` for parallel GPU performance

## Build Instructions

### Prerequisites

* C++20 compiler (e.g., `g++` version supporting `-std=c++20`)
* `xsimd` header-only library (make sure it’s available in `./inc`)
* NVIDIA CUDA Toolkit (for compiling the CUDA version)

### Build

```bash
make
```

### Run Alignment Test

You can test the implementation with:

```bash
make test
```

Or run manually with:

```bash
./main_test seq1.fa seq2.fa
```

You can also use command to generate sequence
```bash
./generate.sh <ref length> <query length>
```

### Clean

To remove compiled objects and executable:

```bash
make clean
```

## A Example Expected Output Format
With reference = 20kbp, query = 20kbp
```
Running test case:
./main_test seq1.fa seq2.fa
==================================
Method: Striped
Seq1:     ...

Seq2:     ...
Alignment Score: 15091
Time: 1972.122222 ms
==================================
Method: CUDA
Seq1:     ...

Seq2:     ...
Alignment Score: 15091
Time: 811.441891 ms
==================================
Speedup (Striped SW / Cuda SW): 2.43x (uaually 2x ~ 3x at this case)
```
---
