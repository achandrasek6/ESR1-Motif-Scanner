# ESR1-Motif-Scanner

This tool scans any input DNA sequence with a user-provided position weight matrix for the ESR1-binding motif, computing and reporting log-odds scores for each 15-mer window on both strands along with the highest-scoring strand at each position.

_Rapid identification of ESR1‐binding motif occurrences in genomic sequences is enabled by this tool, thereby facilitating investigation of estrogen receptor‐mediated transcriptional regulation with implications for breast cancer research and other hormone‐responsive clinical applications._

---

## Background Information
### Estrogen Receptor 1 (ESR1)

The ESR1 gene encodes estrogen receptor α, a ligand-activated transcription factor that mediates the effects of estrogens in a variety of tissues. Upon binding 17β-estradiol, ESR1 undergoes dimerization and recognizes specific 15-base-pair estrogen response elements (EREs) in DNA, recruiting co-activator complexes to modulate transcription of downstream target genes (Smith, Nawaz, & O’Malley, 1994). Aberrant ESR1 signaling is observed in hormone-responsive breast cancers, where ESR1-driven gene expression programs contribute to increased proliferation and survival (Ali & Coombes, 2000).

### Selected Target Genes
Two established ESR1 targets—GREB1 and ADAM6—were examined:

- GREB1 (Growth Regulation by Estrogen in Breast Cancer 1) is an early estrogen-responsive gene whose expression correlates with estrogen-driven cell proliferation in breast and ovarian tissues. Multiple high-affinity EREs within the GREB1 promoter render it a sensitive indicator of ESR1 activity (Carroll et al., 2006).

- ADAM6 belongs to the “A Disintegrin And Metalloproteinase” family. Although less extensively characterized, ADAM6 transcription has been shown to be upregulated by estrogen in certain breast cancer cell lines, and its proteolytic activity may influence remodeling of the tumor microenvironment (McCulloch & Birkedal-Hansen, 2001).

Genomic intervals spanning ±100 kb around each gene’s transcriptional unit were scanned using the esr1_motif_scan.py script. Both forward and reverse strands were interrogated for canonical and variant EREs. The output columns—`fwd_score`, `rev_score`, and `best_score`—provide log-odds enrichment metrics for each 15-mer window, thereby enabling prioritization of candidate regulatory elements for subsequent experimental validation.

---

## Repository Contents

| File                      | Description                                                                                 |
|---------------------------|---------------------------------------------------------------------------------------------|
| `esr1_motif_scan.py`      | Main Python script. Parses a 4×n PWM and a DNA sequence, then slides a 15-mer window and computes log-odds scores for both strands. |
| `esr1_profile.txt`        | Tab- or space-delimited 4×n PWM for the ESR1-binding motif (rows: A, C, G, T).              |
| `GREB1.fasta`             | Genomic sequence (±100 kb flanks) for GREB1 in FASTA format.                                 |
| `ADAM6.fasta`             | Genomic sequence (±100 kb flanks) for ADAM6 in FASTA format.                                 |
| `GREB1_output.tsv`        | Tab-delimited scan results for `GREB1.fasta`.                                               |
| `ADAM6_output.tsv`        | Tab-delimited scan results for `ADAM6.fasta`.                                               |
| `Gene_list.txt`           | List of ESR1-regulated target genes used to select GREB1 and ADAM6. [View Data Source](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ESR1)|
| `file_run_instructions.txt` | High-level instructions for making the script executable and running the scans.             |
| `README.md` | This instructional file.             |

---

## Prerequisites

- **Python 3** (no external packages required; uses only the standard library: `argparse`, `math`, `sys`).
- **Unix-style shell** (Linux, macOS, or WSL on Windows).

---

## Quick Start

1. **Clone or download** this directory into your WSL or Linux environment.

2. **Make the script executable** (once):
   ```bash
   chmod +x esr1_motif_scan.py
3. **Run a scan:**
   ```bash
   ./esr1_motif_scan.py esr1_profile.txt GREB1.fasta > GREB1_output.tsv
   ./esr1_motif_scan.py esr1_profile.txt ADAM6.fasta > ADAM6_output.tsv
4. **Optional background frequencies**
   
   By default, all four nucleotides are assumed to have a background frequency of 0.25. To specify custom frequencies:
   
       ./esr1_motif_scan.py esr1_profile.txt GREB1.fasta --bg 0.30 0.20 0.20 0.30 > GREB1_output.tsv
   where the `--bg` flag takes four floats (A C G T).
5. **Modularity**

   This application may be used to scan for potential motifs in other receptor-binding domain interactions between other proteins. 

       ./esr1_motif_scan.py receptor_motif_profile.txt target_gene.fasta > target_gene_output.tsv
   
   where `receptor_motif_profile.txt` refers to the PWM for a given receptor, `target_gene.fasta` is the selected genetic sequence corresponding to the target protein, and `target_gene_output.tsv` is the log-odds output file.
---

## File Run Instructions
See `file_run_instructions.txt` for a minimal step-by-step guide to reproduce output files in directory.

---

## Output Format
| pos | fwd_score | rev_score | best_score | strand |
|-----|-----------|-----------|------------|--------|
| 1   | 2.345     | 1.123     | 2.345      | +      |
| 2   | 1.789     | 2.001     | 2.001      | -      |
| …   | …         | …         | …          | …      |

- `pos`: 1-based start position of the window in the input sequence
- `fwd_score`: log₂-odds score on the forward strand
- `rev_score`: log₂-odds score on the reverse complement
- `best_score`: the higher of forward/reverse
- `strand`: “+” if forward is higher, “–” if reverse is higher

---

## References

1. Ali, S., & Coombes, R. C. (2000). Estrogen receptor alpha in human breast cancer: occurrence and significance. Journal of Mammary Gland Biology and Neoplasia, 5(3), 271–281.

2. Carroll, J. S., Liu, X. S., Brodsky, A. S., Li, W., Meyer, C. A., Szary, A. J., … Brown, M. (2005). Chromosome-wide mapping of estrogen receptor binding reveals long-range regulation requiring the forkhead protein FOXA1. Cell, 122(1), 33–43.

3. McCulloch, D. R., & Birkedal-Hansen, H. (2001). Structures of metalloproteinases reflect their functions. Matrix Biology, 19(2), 101–112.

4. Smith, C. L., Nawaz, Z., & O’Malley, B. W. (1994). Coactivation of steroid hormone receptors: biological and molecular aspects. Endocrine Reviews, 15(3), 203–211.








