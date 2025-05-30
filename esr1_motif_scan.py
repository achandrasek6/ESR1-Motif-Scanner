#!/usr/bin/env python3
"""
esr1_motif_scan.py

Scan a DNA sequence for matches to a PWM (e.g. the ESR1‐binding motif profile).
Outputs a tab‐delimited table of log‐odds scores for each 15‐mer window.

Usage:
    python esr1_motif_scan.py PROFILE_FILE SEQUENCE_FILE [--bg A C G T]

Arguments:
    PROFILE_FILE   Tab‐ or space‐delimited 4×n PWM (rows: A, C, G, T)
    SEQUENCE_FILE  FASTA or plain text file with the DNA to scan
Options:
    --bg A C G T   Optional background frequencies for A, C, G, T (default 0.25 each)
"""

import argparse
import math
import sys

# Complement map for reverse complement
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def read_pwm(path):
    """Read a 4×n PWM file and return a dict {'A': [..], 'C': [..], 'G': [..], 'T': [..]}."""
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    if len(lines) != 4:
        sys.exit("ERROR: PWM file must have exactly 4 non‐empty lines (A,C,G,T).")
    nucs = ['A','C','G','T']
    pwm = {}
    for nuc, line in zip(nucs, lines):
        vals = line.replace(",", " ").split()
        pwm[nuc] = [float(x) for x in vals]
    # Ensure all rows same length
    lengths = {len(v) for v in pwm.values()}
    if len(lengths) != 1:
        sys.exit("ERROR: PWM rows have inconsistent lengths.")
    return pwm, lengths.pop()

def read_sequence(path):
    """Read a sequence file (FASTA or plain text), concatenate all non‐header lines."""
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq).upper()

def revcomp(seq):
    """Return the reverse complement of the DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]

def score_window(window, pwm, bg):
    """
    Compute log2‐odds score for a window:
      sum_i log2( PWM[base][i] / bg[base] )
    If base not A/C/G/T, uses a tiny pseudo‐probability.
    """
    total = 0.0
    length = len(window)
    for i, b in enumerate(window):
        # handle any unexpected character
        p = pwm.get(b, [1e-6]*length)[i]
        q = bg.get(b, 1e-6)
        total += math.log2(p / q)
    return total

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("profile", help="PWM file (4 rows A,C,G,T)")
    p.add_argument("sequence", help="Sequence file (FASTA or plain text)")
    p.add_argument("--bg", nargs=4, type=float,
                   metavar=("A","C","G","T"),
                   help="Background freqs for A C G T (default=0.25 each)")
    args = p.parse_args()

    # Load PWM and sequence
    pwm, width = read_pwm(args.profile)
    seq = read_sequence(args.sequence)

    # Set up background frequencies
    if args.bg:
        bg = dict(zip(['A','C','G','T'], args.bg))
    else:
        bg = {n: 0.25 for n in "ACGT"}

    # Print header
    print("pos\tfwd_score\trev_score\tbest_score\tstrand")

    # Slide window
    for i in range(len(seq) - width + 1):
        window = seq[i:i+width]
        fwd = score_window(window, pwm, bg)
        rev = score_window(revcomp(window), pwm, bg)
        # choose best strand
        if fwd >= rev:
            best, strand = fwd, "+"
        else:
            best, strand = rev, "-"
        # Report 1‐based position
        print(f"{i+1}\t{fwd:.3f}\t{rev:.3f}\t{best:.3f}\t{strand}")

if __name__ == "__main__":
    main()
