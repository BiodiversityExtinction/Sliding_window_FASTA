# Sliding-window FASTA extraction and concatenation

A lightweight shell pipeline for generating sliding-window FASTA sequences from whole-genome assemblies and preparing datasets for windowed phylogenetic and phylogenomic analyses.

This tool extracts predefined genomic intervals from multiple genome FASTA files, standardises FASTA headers, handles missing data by padding with `N`s, and produces both per-window and concatenated sequences in a deterministic and reproducible way.

---

## Features

- Extracts genomic windows using `samtools faidx`
- Supports many samples via a simple sample table
- Writes per-window FASTAs with informative headers: `>CODENAME|chr:start-end`
- Automatically pads missing windows with `N`s of the correct length
- Produces per-sample concatenated FASTA sequences
- Deterministic window order (defined entirely by the regions file)
- Detailed logging (per sample and per run)
- Minimal dependencies, no workflow framework required

---

## Requirements

- `samtools`
- `seqtk`
- Bash (tested with `bash ≥ 4`)

Genome FASTA files will be indexed automatically if a `.fai` file is not present.

---

## Installation

Clone the repository and make the script executable:

```bash
git clone https://github.com/yourusername/your-repo-name.git
cd your-repo-name
chmod +x make_windows_all.sh
```

No further installation is required.

---

## Inputs

### Sample table (`--input`)

A **tab-separated** file with two columns:

```text
CODENAME<TAB>/path/to/genome.fasta
```

Example:

```text
#codename	genome
Polar	/path/to/Polarbear.fa
Brown	/path/to/BrownBear.fa
Sun	/path/to/SunBear.fa
```

Notes:
- Lines starting with `#` are ignored
- Codenames must not contain whitespace
- Each genome must be a valid FASTA file

---

### Regions file (`--regions`)

A text file containing **one genomic interval per line**, in **samtools faidx format**:

```text
chr:start-end
```

Example:

```text
chr1:1000001-1020000
chr1:2000001-2020000
chr2:500001-520000
```

Important:
- Coordinates must be **1-based and inclusive**
- **Explicit intervals are required** (contig-only entries such as `chr1` are not supported)
- The order of regions in this file defines the concatenation order

---

## Usage

```bash
./make_windows_all.sh \
  --input samples.tsv \
  --regions regions.txt \
  --outdir OUT
```

To view help:

```bash
./make_windows_all.sh --help
```

---

## Outputs

For each sample `CODENAME`, the following files are created:

```text
OUT/
├── run.log
├── CODENAME/
│   ├── chr1:1000001-1020000.fasta
│   ├── chr1:2000001-2020000.fasta
│   ├── CODENAME_concat.fasta
│   └── CODENAME.log
```

### Per-window FASTA files

- One FASTA file per region per sample
- Header format:
  ```
  >CODENAME|chr:start-end
  ```

### Concatenated FASTA

- One concatenated FASTA per sample
- Header format:
  ```
  >CODENAME
  ```
- Sequences are concatenated in the exact order of the regions file

---

## Missing data handling

If a region is not present in a given genome assembly (e.g. missing scaffold or truncated contig), the script writes a sequence of `N`s with length:

```
end - start + 1
```

This ensures consistent alignment lengths across samples and preserves positional correspondence between windows.

All missing regions are recorded in the per-sample log file.

---

## Logging

- `OUT/run.log`  
  Summary of the full run, including start/end times and per-sample status
- `OUT/CODENAME/CODENAME.log`  
  Detailed log for each sample, including missing windows and concatenation statistics

---

## Typical applications

- Sliding-window phylogenetic tree inference
- Genome-wide concatenation analyses
- Comparative genomics across many samples
- Analyses involving ancient, historical, and modern DNA datasets

---

## Notes

- This tool assumes all regions are provided explicitly as `chr:start-end`
- No automatic BED conversion is performed
- Window generation (e.g. sliding windows) should be performed upstream

---

## Citation

If you use this workflow in a publication, please cite this repository and specify the commit hash used.
