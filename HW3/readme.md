- build
```bash
# build
make
# use simple test
make run
# clean obj & exe
make clean
```
- usage
```
Usage: ./aligment [options] <seq1.fasta> <seq2.fasta>

Options:
  -m <int>     Match score (default: 2)
  -x <int>     Mismatch penalty (default: -1)
  -o <int>     Gap open penalty (default: -2)
  -e <int>     Gap extension penalty (default: -1)
  -b <int>     Band width (default: 10)
  -h, --help   Show this help message

# note:
# -b only modify the band width of baseline version
```
- Notes
  - The baseline implementation uses banded Smith-Waterman, which only computes scores near the diagonal. Therefore, when the two sequences differ significantly in length and the band width is small, only the beginning portion of the longer sequence may be computed.
  
  - The SIMD version first uses SIMD instructions to quickly identify the highest-scoring alignment end position, then traces back to find the start position, thereby determining an alignment region. The final alignment within this region is computed using the baseline banded Smith-Waterman. In this case, **the band width is automatically determined** and is set to **the length difference between the two sequences plus 50**.