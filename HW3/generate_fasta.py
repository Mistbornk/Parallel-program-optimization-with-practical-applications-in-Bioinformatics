import random
import argparse

def generate_random_dna(length):
    return ''.join(random.choices('ACGT', k=length))

def write_fasta(filename, seq_id, sequence):
    with open(filename, 'w') as f:
        f.write(f">{seq_id}\n")
        f.write(sequence + "\n")

def generate_two_fasta_files(len1, len2, file1="seq1.fa", file2="seq2.fa"):
    seq1 = generate_random_dna(len1)
    seq2 = generate_random_dna(len2)

    write_fasta(file1, "seq1", seq1)
    write_fasta(file2, "seq2", seq2)

    print(f"Generated {file1} ({len1} bases) and {file2} ({len2} bases).")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate two FASTA files with random DNA sequences.")
    parser.add_argument("-seq1", type=int, required=True, help="Length of the first sequence (seq1)")
    parser.add_argument("-seq2", type=int, required=True, help="Length of the second sequence (seq2)")

    args = parser.parse_args()
    generate_two_fasta_files(args.seq1, args.seq2)
