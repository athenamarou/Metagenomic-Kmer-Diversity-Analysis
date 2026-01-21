"""PROJECT 4: METAGENOMIC DATA ANALYSIS"""

import gzip as gz
import numpy as np
import itertools
from collections import Counter
import matplotlib.pyplot as plt
import argparse

# Setup Argument Parser to accept arguments from the terminal
parser = argparse.ArgumentParser(description="Analyze k-mer diversity in a FASTQ file.")
parser.add_argument("input_file", type=str, help="The input FASTQ file (.fq.gz)")
parser.add_argument(
    "-k",
    type=int,
    choices=[4, 5, 6, 7],
    help="Optional: Specific k-mer size (4-7). If omitted, runs all.",
)

args = parser.parse_args()

# Assigning values from arguments to your variables
file_name = args.input_file

if args.k:
    kmer_size = [args.k]
else:
    kmer_size = [4, 5, 6, 7]

step_size = 1000  # Check every 1000 reads

# Path sanitized: using the variable from arguments instead of the local hardcoded path
# file_name = args.input_file


# Reading file, checking and yielding sequences
def read_fastq_seqs(file):
    with gz.open(file, "rt") as f:
        for seq_header, seq, qual_header, qual in itertools.zip_longest(*[f] * 4):
            if any(line is None for line in (seq_header, seq, qual_header, qual)):
                raise Exception(
                    "Incomplete FASTQ record found. Number of lines is not a multiple of 4."
                )
            if not (seq_header.startswith("@") and qual_header.startswith("+")):
                raise Exception("Invalid FASTQ format. Expected headers not found.")
            yield seq.rstrip("\n")


# FUNCTION FOR K-MER FREQUENCY CALCULATION


# Shannon entropy calculation using counter
def shannon_entropy(counter):
    total = sum(counter.values())
    entropy = 0.0
    for count in counter.values():
        p = count / total
        entropy -= p * np.log2(p)
    return entropy


# One pass over sequences to compute k-mer frequencies
def analyze_kmers(file, kmer_size, step_size):
    kmer_counters = {k: Counter() for k in kmer_size}
    prev_entropy = {k: None for k in kmer_size}
    final_results = {}  # Store final results for each k
    plot_data = {
        k: {"reads": [], "unique_kmers": []} for k in kmer_size
    }  # Data for the plot

    active_k = set(kmer_size)  # Track active k values that haven't converged yet

    reads_proc = 0

    print(f"Initializing k-mer analysis for k={kmer_size} with step {step_size} reads.")

    # Reading file once
    for seq in read_fastq_seqs(file):
        reads_proc += 1

        # Updating counters for each active k
        for k in active_k:
            kmers = (seq[i : i + k] for i in range(len(seq) - k + 1))
            valid_kmers = (kmer for kmer in kmers if "N" not in kmer)
            kmer_counters[k].update(valid_kmers)

        # Periodically (every 1000 reads)
        if reads_proc % step_size == 0:

            # Plot data
            for k in kmer_size:
                plot_data[k]["reads"].append(reads_proc)
                plot_data[k]["unique_kmers"].append(len(kmer_counters[k]))

            # List of k values that have converged and will be removed
            converged_k = []
            for k in active_k:
                current_entropy = shannon_entropy(kmer_counters[k])

                # Convergence check
                if prev_entropy[k] is not None:
                    diff = abs(current_entropy - prev_entropy[k])

                    if diff < 0.00001:
                        print(
                            f"Convergence reached for k={k} at {reads_proc} reads with entropy difference {diff}"
                        )

                        # Saving final results
                        total_kmers = sum(kmer_counters[k].values())
                        unique_kmers = len(kmer_counters[k])
                        final_results[k] = {
                            "Final Shannon Entropy": current_entropy,
                            "Reads Processed": reads_proc,
                            "Total kmers": total_kmers,
                            "Unique kmers": unique_kmers,
                        }
                        converged_k.append(k)
                # Updating previous entropy value
                prev_entropy[k] = current_entropy

            # Removing converged k values
            for k in converged_k:
                active_k.remove(k)

            # Stop if all k values have converged
            if not active_k:
                print("\nAll k values have converged. Stop analysis.")
                break

    return final_results, plot_data


# EXECUTION OF ANALYSIS
results, plot_history = analyze_kmers(file_name, kmer_size, step_size)

print("\nFinal Results:")

print("\n" + "=" * 65)
print(f"{'k':<5} {'Shannon Entropy':<20} {'Reads Processed':<20} {'Total k-mers':<15}")
print("-" * 65)

for k in sorted(results.keys()):
    res = results[k]
    # Formatting results table
    print(
        f"{k:<5} {res['Final Shannon Entropy']:<20.5f} {res['Reads Processed']:<20,} {res['Total kmers']:<15,}"
    )
print("=" * 65 + "\n")

# GRAPHICAL REPRESENTATION: X-axis (reads processed) and Y-axis (unique k-mers)

for k in results.keys():
    plt.figure(figsize=(10, 6))
    unique_kmers = plot_history[k]["unique_kmers"]
    reads_processed = plot_history[k]["reads"]
    plt.plot(reads_processed, unique_kmers, linewidth=2, label=f"k={k}", color="teal")
    plt.title(
        f"Unique k-mers vs Reads Processed (k={k}) reflecting metagenome richness"
    )
    plt.xlabel("Reads Processed")
    plt.ylabel("Unique kmers")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.savefig(f"unique_kmers_vs_reads_processed_k{k}.png")
    plt.show()

"""Rarefaction analysis showed rapid saturation (plateau) for k-mer sizes 4 to 7. 
Given that the sample is from the human gut microbiome (high complexity), 
rapid saturation at k=7 is expected as the potential combination space (4^7 ≈ 16,000) 
is fully covered even by a small number of bacterial genomes. 
This result confirms the high data quality (lack of noise) 
and sufficient sequencing depth to recover complete information at the short sequence level."""
