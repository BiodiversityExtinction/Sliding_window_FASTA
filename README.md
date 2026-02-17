# PhyloSlide

PhyloSlide is a sliding-window phylogenomics pipeline for extracting genomic windows from multiple individuals and (optionally) inferring per-window phylogenetic trees in parallel across windows. It includes filters for missing data and low-information windows, optional aDNA-oriented masking (transversions-only), concordance factor analysis, and topology-based window selection for downstream dating analyses.

---

## What PhyloSlide Does

### Extraction Stage (always runs)

- Reads a sample table (`CODENAME<TAB>genome.fasta`)
- Uses a regions list in samtools faidx coordinate format (`chr:start-end`, 1-based inclusive)
- For each sample and region:
  - extracts the window using `samtools faidx`
  - pads with Ns if the region is missing
  - writes window FASTAs with header format: `>CODENAME|chr:start-end`
- Writes a per-sample concatenated sequence across all regions:
  - `OUT/<CODENAME>/<CODENAME>_concat.fasta`

### Optional Tree Pipeline (`--runtrees`)

- Filters windows by missing data (`--maxN`)
- Builds per-window multi-FASTA alignments
- Optional conservative aDNA mode (`--transversions`)
  - masks columns containing A/G or C/T variation to N
  - ensures all downstream alignments are built consistently from masked data
- Filters windows by minimum parsimony-informative sites (`--minpi`)
- Runs one IQ-TREE job per window in parallel (`--jobs`)
- Builds a reference tree:
  - concatenated IQ-TREE (`--ref concat`)
  - or ASTRAL species tree (`--ref astral`)
- Computes gCF and sCF using IQ-TREE2

### Optional Topology Filtering (`--topofilter`)

Builds a dating supermatrix using only windows that:
- passed missingness and minPI filters
- match the reference topology (exact or compatible)
- have minimum internal bootstrap ≥ `--minbs`

---

## Installation

### Recommended: Conda / Mamba Environment

Create a file named `environment.yml`:

```yaml
name: phyloslide
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.10
  - biopython>=1.80
  - samtools
  - seqtk
  - bedtools
  - iqtree
  - openjdk
  - astral-tree
```

Create and activate the environment:

```bash
mamba env create -f environment.yml
mamba activate phyloslide
```

If using conda instead of mamba:

```bash
conda env create -f environment.yml
conda activate phyloslide
```

Notes:
- `bedtools` is required only for `--makewindows`
- `astral-tree` and `openjdk` are required only for `--ref astral`
- `biopython` is required for `--ref astral` and `--topofilter`

---

## Input Format

### Sample Table (`--input`)

Tab-separated file:

```
American_black    /path/to/American_black.fa
Asian_black       /path/to/Asian_black.fa
Brown             /path/to/Brown.fa
Polar             /path/to/Polar.fa
```

- No spaces in codenames
- FASTA files must be indexed (`samtools faidx`) — PhyloSlide will create `.fai` if missing

### Regions File (`--regions`)

Each line must be:

```
chr:start-end
```

Coordinates must be 1-based inclusive.

Example:

```
chr1:1-20000
chr1:1000001-1020000
```

---

## Generating Windows Internally (`--makewindows`)

Instead of manually creating a regions file, PhyloSlide can generate one using `bedtools makewindows`.

Example: 20 kb windows sliding every 1 Mb:

```bash
python phyloslide.py \
  --input samples.tsv \
  --outdir OUT \
  --makewindows \
  --refgenome ref.fa \
  --window 20000 \
  --step 1000000
```

Restrict to long scaffolds:

```bash
python phyloslide.py \
  --input samples.tsv \
  --outdir OUT \
  --makewindows \
  --refgenome ref.fa \
  --window 20000 \
  --step 1000000 \
  --min_scaffold_len 14000000
```

Restrict to specific chromosomes:

```bash
python phyloslide.py \
  --input samples.tsv \
  --outdir OUT \
  --makewindows \
  --refgenome ref.fa \
  --window 20000 \
  --step 1000000 \
  --chroms chroms.txt
```

`chroms.txt` should contain one contig name per line.

---

## Running the Tree Pipeline

Example full workflow:

```bash
python phyloslide.py \
  --input samples.tsv \
  --outdir OUT \
  --makewindows --refgenome ref.fa --window 20000 --step 1000000 \
  --runtrees \
  --jobs 32 \
  --iqtree_threads_per_job 1
```

Parallel recommendations:
- Prefer `--iqtree_threads_per_job 1`
- Set `--jobs` equal to available cores (or slightly below)
- Reduce `--jobs` if running on slow network storage

---

## aDNA Conservative Mode (`--transversions`)

```bash
python phyloslide.py \
  --input samples.tsv \
  --regions regions.txt \
  --outdir OUT \
  --runtrees \
  --transversions
```

Behavior:
- Masks entire columns that contain both A and G OR both C and T
- All downstream alignments remain consistent with TV-only filtering

---

## Filters

### Missing Data Filter (`--maxN`)
- Default: 0.5
- Drops window if ANY sample has fraction of N > maxN

### Minimum Parsimony-Informative Sites (`--minpi`)
- Default:
  - 20 (full data)
  - 10 (with `--transversions`)
- Applied after TV masking (if enabled)

---

## Reference Trees

### Concatenated (default)

```
--ref concat
```

Runs IQ-TREE on:
- `All_concat.filtered.fasta`
- or `All_concat.filtered.tv.fasta`

### ASTRAL

```
--ref astral --astral /path/to/astral.jar
```

- Collapses branches with support < `--collapse`
- Requires Java
- Requires biopython

---

## Concordance Factors

When `--runtrees` is enabled:
- gCF and sCF are computed using IQ-TREE2
- Outputs in `OUT/trees/concordance/`

---

## Topology Filtering

```bash
python phyloslide.py \
  --input samples.tsv \
  --regions regions.txt \
  --outdir OUT \
  --runtrees \
  --topofilter \
  --minbs 90 \
  --topomode exact
```

Options:
- `--minbs` default 90
- `--topomode exact` (default)
- `--topomode compatible`

Outputs:
- Filtered region list
- Dating supermatrix alignment

---

## Output Structure

```
OUT/
  phyloslide.log
  filtering/
  Combined_windows/
  Combined/
  trees/
  <CODENAME>/
```

Key files:
- `regions.kept.final.txt`
- `All_concat.filtered.fasta`
- `all_window_trees.trs`
- Reference tree
- Concordance outputs

---

## Troubleshooting

If slow:
- Reduce `--jobs`
- Run on local scratch disk

If no windows left:
- Increase `--maxN`
- Decrease `--minpi`
- Increase window size
- Reduce `--minbs`

If IQ-TREE fails:
- Check `trees/window_logs/`
- Check `failed_windows.txt`

Windows filenames contain colons (`chr:start-end`).
Use Linux/macOS; Windows may fail due to colon in filenames.

---

## Recommended Student Workflow

1. Test on ~50 windows first.
2. Confirm filtering behavior.
3. Scale up gradually.
4. Record:
   - window parameters
   - filtering parameters
   - model and bootstrap settings
   - reference method used

---

## Dependencies Summary

Always required:
- samtools
- seqtk

If using:
- `--makewindows` → bedtools
- `--runtrees` → iqtree
- `--ref astral` → java + astral-tree
- `--topofilter` or `--ref astral` → biopython

---

## Citation

Please cite:
- samtools
- bedtools
- IQ-TREE2
- ASTRAL

PhyloSlide orchestrates these tools but does not replace them.
