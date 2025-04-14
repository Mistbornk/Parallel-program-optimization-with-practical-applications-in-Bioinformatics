# Striped Smith-Waterman (XSIMD) Alignment
github repository: https://github.com/Mistbornk/Parallel-program-optimization-with-practical-applications-in-Bioinformatics

This project implements the Smith-Waterman local alignment algorithm in two modes:

- **Baseline (non-SIMD) version**
- **XSIMD-accelerated version** using `xsimd` for vectorized performance

## Build Instructions

### Prerequisites
- C++20 compiler (e.g., `g++` version supporting `-std=c++20`)
- `xsimd` header-only library (make sure it’s available in `./inc`)

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

### Clean
To remove compiled objects and executable:
```bash
make clean
```

## Expected Output Format

```
==================================
Match score truth: 184
==================================
Mode: Baseline
Seq1:     21   CGGATCCGGTAGCTAGCTACCGTTAGGCCTAGGTTAGGCTAGCTAGGTTACCTAGGATCGATCGTACGTTAGCTAGGGTACCGTAGCTA-GCTAGGTCGATCGTACG-TAG   129
               |||||||||*||||||||| |||||||||*||||||||  |||||||||||||||||||**|||||||||||||||||||||||||||| |*|*|*|||**||| || |||
Seq2:      1   CGGATCCGGCAGCTAGCTA-CGTTAGGCCAAGGTTAGG--AGCTAGGTTACCTAGGATCAGTCGTACGTTAGCTAGGGTACCGTAGCTAGGGTCGATCGTACGT-CGCTAG   107
Alignment Score: 158
Execution time: 0.000112 seconds
==================================
Mode: XSIMD
Seq1:     20   CGGATCCGGTAGCTAGCTACCGTTAGGCCTAGGTTAGGCTAGCTAGGTTACCTAGGATCGATCGTACGTTAGCTAGGGTACCGTAGCTAGCTAGGTCGATCGTACGTAGCTAGG   134
               |||||||||*||||||||| |||||||||*||||||||  |||||||||||||||||||**|||||||||||||||||||||||||||||   ||||||||||||||*||||||
Seq2:      0   CGGATCCGGCAGCTAGCTA-CGTTAGGCCAAGGTTAGG--AGCTAGGTTACCTAGGATCAGTCGTACGTTAGCTAGGGTACCGTAGCTAG---GGTCGATCGTACGTCGCTAGG   108
Alignment Score: 184
Execution time: 0.000022 seconds
==================================
Speedup       : 5.033744x
```

---