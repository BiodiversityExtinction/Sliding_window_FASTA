#!/usr/bin/env python3
"""
PhyloSlide: sliding-window phylogenomics pipeline

What it does
------------
1) (Optional) Generate sliding windows from a reference FASTA index using bedtools makewindows
   and write regions as samtools faidx-style coordinates: chr:start-end (1-based inclusive).

2) Extract per-sample window FASTAs with samtools faidx (pads missing regions with Ns).
   Window header format:   >CODENAME|chr:start-end
   Per-sample concat:      >CODENAME  (concatenation of window sequences in regions-file order)

3) (Optional, --runtrees) Build per-window multi-FASTA alignments and:
   - filter windows by missing data (drop if ANY sample has N fraction > --maxN)
   - optionally mask transition-containing columns (TV-only) with --transversions
   - filter windows by minimum parsimony-informative sites (PI) with --minpi
   - infer window trees with IQ-TREE2 in parallel across windows (--jobs)
   - build reference tree (concat or ASTRAL) and compute gCF/sCF
   - optionally build a “dating” supermatrix with windows matching the reference topology
     and with minimum bootstrap support threshold (--topofilter)

Dependencies
------------
External tools (required depending on flags):
  Always:
    - samtools
    - seqtk
  If --makewindows:
    - bedtools
  If --runtrees:
    - iqtree2
  If --ref astral:
    - java
    - ASTRAL jar file (provided via --astral)
  If --topofilter or --ref astral:
    - Python package: biopython

Python:
  - Python >= 3.9 recommended
  - biopython required only for topology matching and ASTRAL collapsing steps

Author: (you)
License: (choose; e.g., MIT)
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# -----------------------------
# Logging / shell execution
# -----------------------------
def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, logfile: Path) -> None:
    line = f"[{now()}] {msg}"
    print(line)
    logfile.parent.mkdir(parents=True, exist_ok=True)
    with logfile.open("a") as f:
        f.write(line + "\n")


def which(exe: str) -> Optional[str]:
    return shutil.which(exe)


def which_or_die(exe: str, problems: List[str]) -> None:
    if which(exe) is None:
        problems.append(f"Missing external dependency in PATH: {exe}")


def run_cmd(cmd: List[str], logfile: Path, cwd: Optional[Path] = None) -> None:
    logfile.parent.mkdir(parents=True, exist_ok=True)
    with logfile.open("a") as f:
        f.write(f"[{now()}] CMD: {' '.join(cmd)}\n")
        f.flush()
        p = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=f,
            stderr=f,
            text=True,
        )
        if p.returncode != 0:
            raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}")


# -----------------------------
# FASTA helpers
# -----------------------------
def fasta_read_one(path: Path) -> Tuple[str, str]:
    header = None
    seq_parts: List[str] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()[0]
            else:
                seq_parts.append(line)
    if header is None:
        raise ValueError(f"No FASTA header found in {path}")
    return header, "".join(seq_parts)


def fasta_write_one(path: Path, header: str, seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(f">{header}\n{seq}\n")


def nfrac_seq(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 1.0
    n = sum(1 for c in seq if c == "N")
    return n / len(seq)


def sanitize_region(region: str) -> str:
    return region.replace(":", "_").replace("/", "_")


# -----------------------------
# Regions / windows
# -----------------------------
REGION_RE = re.compile(r"^([^:]+):(\d+)-(\d+)$")


def parse_region(region: str) -> Tuple[str, int, int]:
    m = REGION_RE.match(region.strip())
    if not m:
        raise ValueError(f"Invalid region format (must be chr:start-end): {region}")
    chrom, s, e = m.group(1), int(m.group(2)), int(m.group(3))
    if e < s:
        raise ValueError(f"Region end < start: {region}")
    return chrom, s, e


def region_len(region: str) -> int:
    _, s, e = parse_region(region)
    return e - s + 1


def read_regions(path: Path) -> List[str]:
    regs: List[str] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parse_region(line)  # validate
            regs.append(line)
    if not regs:
        raise ValueError(f"No regions found in file: {path}")
    return regs


def read_samples_table(path: Path) -> List[Tuple[str, Path]]:
    out: List[Tuple[str, Path]] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"Bad samples line (need 2 columns): {line}")
            code = parts[0]
            fasta = Path(parts[1])
            out.append((code, fasta))
    if not out:
        raise ValueError(f"No samples found in file: {path}")
    return out


def generate_regions_with_bedtools(
    outdir: Path,
    runlog: Path,
    refgenome: Optional[Path],
    fai: Optional[Path],
    window: int,
    step: int,
    min_len: int,
    chroms_file: Optional[Path],
    regions_out: Optional[Path],
) -> Path:
    """
    Generate regions (chr:start-end, 1-based inclusive) using bedtools makewindows.

    Source of lengths:
      - If --fai is provided: use it directly
      - Else, require --refgenome and create <refgenome>.fai via samtools faidx if missing

    Filters:
      - min_len: only include contigs/scaffolds with length >= min_len
      - chroms_file: optional whitelist of contig names (one per line)
    """
    if fai is None:
        if refgenome is None:
            raise SystemExit("ERROR: --makewindows requires --refgenome or --fai")
        if not refgenome.exists():
            raise SystemExit(f"ERROR: --refgenome not found: {refgenome}")
        fai = Path(str(refgenome) + ".fai")
        if not fai.exists():
            log("Indexing --refgenome for window generation (samtools faidx)", runlog)
            run_cmd(["samtools", "faidx", str(refgenome)], runlog)

    if not fai.exists():
        raise SystemExit(f"ERROR: FASTA index (.fai) not found: {fai}")

    allowed: Optional[set[str]] = None
    if chroms_file is not None:
        if not chroms_file.exists():
            raise SystemExit(f"ERROR: --chroms file not found: {chroms_file}")
        allowed = {ln.strip() for ln in chroms_file.read_text().splitlines() if ln.strip()}

    lengths_path = outdir / f"lengths.min{min_len}.txt"
    kept_contigs = 0
    with lengths_path.open("w") as out:
        for ln in fai.read_text().splitlines():
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            name = parts[0]
            length = int(parts[1])
            if length < min_len:
                continue
            if allowed is not None and name not in allowed:
                continue
            out.write(f"{name}\t{length}\n")
            kept_contigs += 1

    if kept_contigs == 0:
        raise SystemExit(
            "ERROR: No contigs left after applying --min_scaffold_len / --chroms filters."
        )

    if regions_out is None:
        regions_out = outdir / f"regions.w{window}.s{step}.min{min_len}.txt"

    log(
        f"Generating windows with bedtools makewindows (w={window}, step={step}), contigs={kept_contigs}",
        runlog,
    )

    # bedtools outputs BED (0-based, half-open). Convert to 1-based inclusive.
    p = subprocess.run(
        ["bedtools", "makewindows", "-g", str(lengths_path), "-w", str(window), "-s", str(step)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if p.returncode != 0:
        raise SystemExit("ERROR: bedtools makewindows failed:\n" + p.stderr)

    out_lines: List[str] = []
    for ln in p.stdout.splitlines():
        if not ln.strip():
            continue
        chrom, start0, end0 = ln.split("\t")[:3]
        s1 = int(start0) + 1
        e1 = int(end0)
        out_lines.append(f"{chrom}:{s1}-{e1}")

    if not out_lines:
        raise SystemExit("ERROR: bedtools produced zero windows (check window/step/filters).")

    regions_out.write_text("\n".join(out_lines) + "\n")
    log(f"Wrote generated regions: {regions_out} (windows={len(out_lines)})", runlog)
    return regions_out


# -----------------------------
# Alignment stats / masking
# -----------------------------
def count_parsimony_informative(seqs: List[str]) -> int:
    """
    Parsimony-informative site:
      - ignore N and gaps
      - at least 2 states, and at least 2 states each appear in >=2 taxa
    """
    if not seqs:
        return 0
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise ValueError("Sequences are not all the same length")
    pi = 0
    for i in range(L):
        counts: Dict[str, int] = {}
        for s in seqs:
            b = s[i]
            if b in "ACGT":
                counts[b] = counts.get(b, 0) + 1
        if len(counts) < 2:
            continue
        if sum(1 for v in counts.values() if v >= 2) >= 2:
            pi += 1
    return pi


def read_msa_fasta(msa_path: Path) -> Tuple[List[str], List[str]]:
    headers: List[str] = []
    seqs: List[str] = []
    with msa_path.open() as f:
        h = None
        buf: List[str] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    headers.append(h)
                    seqs.append("".join(buf).upper())
                h = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if h is not None:
            headers.append(h)
            seqs.append("".join(buf).upper())
    return headers, seqs


def write_msa_fasta(msa_path: Path, headers: List[str], seqs: List[str]) -> None:
    with msa_path.open("w") as out:
        for h, s in zip(headers, seqs):
            out.write(f">{h}\n{s}\n")


def mask_transitions_inplace(msa_path: Path) -> None:
    """
    Conservative transversions-only masking:
    For each column, if it contains both A and G OR both C and T (ignoring N/gaps),
    mask the entire column to N across all taxa.
    """
    headers, seqs = read_msa_fasta(msa_path)
    if not seqs:
        return
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise ValueError(f"Unequal sequence lengths in {msa_path}")
    seq_lists = [list(s) for s in seqs]
    for i in range(L):
        col = {seq_lists[j][i] for j in range(len(seq_lists)) if seq_lists[j][i] in {"A", "C", "G", "T"}}
        if ("A" in col and "G" in col) or ("C" in col and "T" in col):
            for j in range(len(seq_lists)):
                seq_lists[j][i] = "N"
    write_msa_fasta(msa_path, headers, ["".join(x) for x in seq_lists])


# -----------------------------
# Parallel IQ-TREE
# -----------------------------
@dataclass(frozen=True)
class WindowJob:
    region: str
    msa_path: Path
    out_prefix: Path
    iqtree_threads: int
    model: str
    bootstrap: int
    log_path: Path


def iqtree_window(job: WindowJob) -> Tuple[str, bool, str]:
    """
    Resumable: if .treefile exists and non-empty, skip.
    """
    treefile = job.out_prefix.with_suffix(".treefile")
    if treefile.exists() and treefile.stat().st_size > 0:
        return (job.region, True, "SKIP(existing treefile)")

    cmd = [
        "iqtree2",
        "-s", str(job.msa_path),
        "-m", job.model,
        "-bb", str(job.bootstrap),
        "-nt", str(job.iqtree_threads),
        "-pre", str(job.out_prefix),
    ]
    try:
        run_cmd(cmd, job.log_path)
        if treefile.exists() and treefile.stat().st_size > 0:
            return (job.region, True, "OK")
        return (job.region, False, "FAILED(no treefile produced)")
    except Exception as e:
        return (job.region, False, f"FAILED({e})")


# -----------------------------
# Preflight checks
# -----------------------------
def preflight_checks(args: argparse.Namespace) -> None:
    problems: List[str] = []

    # Required inputs
    if not args.input.exists():
        problems.append(f"--input not found: {args.input}")
    if args.regions is not None and not args.regions.exists():
        problems.append(f"--regions not found: {args.regions}")

    if args.makewindows:
        if args.fai is None and args.refgenome is None:
            problems.append("--makewindows requires --refgenome or --fai")
        if args.refgenome is not None and not args.refgenome.exists():
            problems.append(f"--refgenome not found: {args.refgenome}")
        if args.fai is not None and not args.fai.exists():
            problems.append(f"--fai not found: {args.fai}")
        if args.chroms is not None and not args.chroms.exists():
            problems.append(f"--chroms not found: {args.chroms}")

    # Regions required unless makewindows
    if not args.makewindows and args.regions is None:
        problems.append("You must provide --regions unless using --makewindows")

    # Always-required external tools
    which_or_die("samtools", problems)
    which_or_die("seqtk", problems)

    # Additional tools
    if args.makewindows:
        which_or_die("bedtools", problems)

    if args.runtrees:
        which_or_die("iqtree2", problems)

    if args.ref == "astral":
        which_or_die("java", problems)
        if args.astral is None:
            problems.append("--ref astral requires --astral /path/to/astral.jar")
        elif not args.astral.exists():
            problems.append(f"--astral jar not found: {args.astral}")

    # Python deps (only if needed)
    if args.runtrees and (args.topofilter or args.ref == "astral"):
        try:
            import Bio  # noqa: F401
        except Exception:
            problems.append("Python dependency missing: biopython (pip install biopython)")

    # Numeric sanity
    if args.window <= 0:
        problems.append("--window must be > 0")
    if args.step <= 0:
        problems.append("--step must be > 0")
    if not (0.0 <= args.maxN <= 1.0):
        problems.append("--maxN must be between 0 and 1")
    if args.minpi < 0:
        problems.append("--minpi must be >= 0")
    if args.bootstrap <= 0:
        problems.append("--bootstrap must be > 0")
    if args.iqtree_threads_per_job <= 0:
        problems.append("--iqtree_threads_per_job must be > 0")
    if args.jobs is not None and args.jobs <= 0:
        problems.append("--jobs must be > 0")
    if args.minbs < 0:
        problems.append("--minbs must be >= 0")
    if args.collapse < 0:
        problems.append("--collapse must be >= 0")

    if problems:
        msg = "Preflight check failed:\n" + "\n".join(f"  - {p}" for p in problems)
        raise SystemExit(msg)


# -----------------------------
# Main
# -----------------------------
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="phyloslide",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "PhyloSlide: sliding-window phylogenomics pipeline.\n\n"
            "Typical workflows:\n"
            "  (A) Use your own regions file (chr:start-end):\n"
            "      phyloslide --input samples.tsv --regions regions.txt --outdir OUT\n\n"
            "  (B) Generate regions internally (20 kb windows, 1 Mb step):\n"
            "      phyloslide --input samples.tsv --outdir OUT --makewindows --refgenome ref.fa \\\n"
            "                --window 20000 --step 1000000 --min_scaffold_len 14000000\n\n"
            "  Add tree inference (parallel across windows):\n"
            "      phyloslide ... --runtrees --jobs 32 --iqtree_threads_per_job 1\n\n"
            "  aDNA-style transversions-only masking (TV-only, conservative):\n"
            "      phyloslide ... --runtrees --transversions\n\n"
            "  Create a dating supermatrix from windows matching the reference topology:\n"
            "      phyloslide ... --runtrees --topofilter --minbs 90 --topomode exact\n"
        ),
    )

    # Core inputs
    ap.add_argument(
        "--input",
        required=True,
        type=Path,
        help=(
            "Sample table (TSV) with two columns:\n"
            "  CODENAME<TAB>/path/to/genome.fasta\n"
            "Lines starting with # are ignored."
        ),
    )
    ap.add_argument(
        "--regions",
        type=Path,
        default=None,
        help=(
            "Regions file (one per line) in samtools faidx coordinate format:\n"
            "  chr:start-end\n"
            "Coordinates must be 1-based inclusive.\n"
            "If omitted, you must use --makewindows to generate regions internally."
        ),
    )
    ap.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory (will be created if missing).",
    )

    # Window generation
    ap.add_argument(
        "--makewindows",
        action="store_true",
        help=(
            "Generate a regions file internally using bedtools makewindows.\n"
            "Requires --refgenome or --fai, plus bedtools in PATH.\n"
            "Output regions are written as chr:start-end (1-based inclusive)."
        ),
    )
    ap.add_argument(
        "--refgenome",
        type=Path,
        default=None,
        help=(
            "Reference genome FASTA used ONLY for window generation (bedtools makewindows).\n"
            "PhyloSlide will create <refgenome>.fai using samtools faidx if missing.\n"
            "Required for --makewindows unless --fai is provided."
        ),
    )
    ap.add_argument(
        "--fai",
        type=Path,
        default=None,
        help=(
            "FASTA index file (.fai) used for window generation (alternative to --refgenome).\n"
            "Format: contig<TAB>length<...>\n"
        ),
    )
    ap.add_argument(
        "--window",
        type=int,
        default=20000,
        help="Window size for --makewindows (bp). Default: 20000.",
    )
    ap.add_argument(
        "--step",
        type=int,
        default=1000000,
        help="Step size (slide) for --makewindows (bp). Default: 1000000.",
    )
    ap.add_argument(
        "--min_scaffold_len",
        type=int,
        default=0,
        help="Only include contigs/scaffolds with length >= this value when generating windows. Default: 0.",
    )
    ap.add_argument(
        "--chroms",
        type=Path,
        default=None,
        help="Optional whitelist of contig names to include in window generation (one name per line).",
    )
    ap.add_argument(
        "--regions_out",
        type=Path,
        default=None,
        help=(
            "Optional output path for generated regions.\n"
            "If not provided, defaults to:\n"
            "  <outdir>/regions.w<window>.s<step>.min<min_scaffold_len>.txt"
        ),
    )

    # Tree pipeline
    ap.add_argument(
        "--runtrees",
        action="store_true",
        help=(
            "Run tree pipeline after extraction:\n"
            "  - missingness filter (--maxN)\n"
            "  - optional TV-only masking (--transversions)\n"
            "  - min parsimony-informative sites filter (--minpi)\n"
            "  - IQ-TREE2 per window (parallel across windows)\n"
            "  - reference tree (concat or ASTRAL)\n"
            "  - gCF/sCF (IQ-TREE2)"
        ),
    )
    ap.add_argument(
        "--transversions",
        action="store_true",
        help=(
            "Conservative transversions-only masking (aDNA-friendly).\n"
            "For each alignment column, if both A and G OR both C and T are present (ignoring N/gaps),\n"
            "mask the entire column to N across all taxa.\n"
            "When enabled, ALL downstream alignments (concat reference / sCF / dating supermatrix)\n"
            "are built from the TV-masked window alignments (consistent TV-only analysis)."
        ),
    )
    ap.add_argument(
        "--maxN",
        type=float,
        default=0.5,
        help=(
            "Missing data filter: drop a window if ANY sample has fraction of 'N' > maxN.\n"
            "Computed on per-sample extracted window sequences (before TV masking).\n"
            "Default: 0.5"
        ),
    )
    ap.add_argument(
        "--minpi",
        type=int,
        default=None,
        help=(
            "Minimum number of parsimony-informative (PI) sites required per window.\n"
            "Computed on the per-window multi-FASTA alignment AFTER TV masking (if enabled).\n"
            "If not provided, defaults to:\n"
            "  - 20 (full data)\n"
            "  - 10 (with --transversions)"
        ),
    )

    # Reference tree / ASTRAL
    ap.add_argument(
        "--ref",
        choices=["concat", "astral"],
        default="concat",
        help=(
            "Reference tree method:\n"
            "  concat : IQ-TREE2 on concatenated alignment (default)\n"
            "  astral : ASTRAL species tree from (optionally collapsed) window trees"
        ),
    )
    ap.add_argument(
        "--astral",
        type=Path,
        default=None,
        help="Path to ASTRAL jar (required if --ref astral).",
    )
    ap.add_argument(
        "--collapse",
        type=int,
        default=30,
        help=(
            "For --ref astral: collapse gene-tree branches with support < collapse before running ASTRAL.\n"
            "Default: 30"
        ),
    )

    # IQ-TREE parameters
    ap.add_argument(
        "--model",
        default="GTR+R6",
        help="IQ-TREE2 substitution model used for window trees and concat reference. Default: GTR+R6",
    )
    ap.add_argument(
        "--bootstrap",
        type=int,
        default=1000,
        help="IQ-TREE2 ultrafast bootstrap replicates (-bb). Default: 1000",
    )
    ap.add_argument(
        "--scf",
        type=int,
        default=100,
        help="IQ-TREE2 site concordance factor replicates (--scf). Default: 100",
    )

    # Parallelism
    ap.add_argument(
        "--iqtree_threads_per_job",
        type=int,
        default=1,
        help=(
            "Threads per IQ-TREE2 window job (-nt). Default: 1\n"
            "For window-level parallelism, keep this small (1–2) and increase --jobs."
        ),
    )
    ap.add_argument(
        "--jobs",
        type=int,
        default=None,
        help=(
            "Number of windows to run concurrently (window-level parallelism).\n"
            "If omitted, auto-selected as min(32, floor(cores / iqtree_threads_per_job))."
        ),
    )

    # Topology filtering / dating supermatrix
    ap.add_argument(
        "--topofilter",
        action="store_true",
        help=(
            "Create a filtered concatenated alignment (dating supermatrix) using only windows that:\n"
            "  (1) passed missingness and minPI filters\n"
            "  (2) match the reference topology (unrooted; see --topomode)\n"
            "  (3) have minimum internal bootstrap >= --minbs\n"
            "Requires biopython."
        ),
    )
    ap.add_argument(
        "--topomode",
        choices=["exact", "compatible"],
        default="exact",
        help=(
            "Topology matching mode:\n"
            "  exact       : gene-tree splits must equal reference splits (default)\n"
            "  compatible  : gene-tree splits must be a subset of reference splits (allows unresolved trees)"
        ),
    )
    ap.add_argument(
        "--minbs",
        type=int,
        default=90,
        help="Minimum internal bootstrap support required for topofilter windows. Default: 90",
    )

    return ap


def main() -> None:
    ap = build_argparser()
    args = ap.parse_args()

    # Default minPI if not provided
    if args.minpi is None:
        args.minpi = 10 if args.transversions else 20

    # Create outdir + runlog early (but do not do any real work before preflight)
    args.outdir.mkdir(parents=True, exist_ok=True)
    runlog = args.outdir / "phyloslide.log"

    # Preflight checks must happen BEFORE running anything substantive.
    preflight_checks(args)

    log("Starting PhyloSlide", runlog)
    log(f"Command line: {' '.join(sys.argv)}", runlog)

    # Parallel design
    C = os.cpu_count() or 1
    T = max(1, int(args.iqtree_threads_per_job))
    auto_jobs = max(1, C // T)
    W = args.jobs if args.jobs is not None else min(32, auto_jobs)
    log(f"Parallel design: cores={C}, threads/job={T}, concurrent windows={W}", runlog)

    # If requested, generate regions file now (this is “real work”, but still before extraction).
    if args.makewindows:
        args.regions = generate_regions_with_bedtools(
            outdir=args.outdir,
            runlog=runlog,
            refgenome=args.refgenome,
            fai=args.fai,
            window=args.window,
            step=args.step,
            min_len=args.min_scaffold_len,
            chroms_file=args.chroms,
            regions_out=args.regions_out,
        )

    if args.regions is None:
        raise SystemExit("ERROR: No regions file available (provide --regions or use --makewindows).")

    # Load inputs
    samples = read_samples_table(args.input)
    regions = read_regions(args.regions)
    codenames = [c for c, _ in samples]

    log(f"Samples: {len(samples)}", runlog)
    log(f"Regions: {len(regions)}", runlog)

    # -----------------------------
    # Step 1: Extract per-sample windows + per-sample concat
    # -----------------------------
    log("Step 1: Extract per-sample windows (samtools faidx) + per-sample concat", runlog)

    for code, genome in samples:
        if not genome.exists():
            raise SystemExit(f"ERROR: genome not found for {code}: {genome}")

        sample_dir = args.outdir / code
        sample_dir.mkdir(parents=True, exist_ok=True)
        slog = sample_dir / f"{code}.log"

        # Ensure .fai
        fai = Path(str(genome) + ".fai")
        if not fai.exists():
            log(f"{code}: indexing genome with samtools faidx", runlog)
            run_cmd(["samtools", "faidx", str(genome)], slog)

        tmp = sample_dir / ".tmp.faidx"

        # Extract each region (skip if already exists)
        missing_padded = 0
        extracted = 0
        for r in regions:
            out_fa = sample_dir / f"{r}.fasta"
            if out_fa.exists() and out_fa.stat().st_size > 0:
                continue

            L = region_len(r)
            # samtools faidx region -> tmp file
            with tmp.open("w") as ftmp:
                p = subprocess.run(
                    ["samtools", "faidx", str(genome), r],
                    stdout=ftmp,
                    stderr=subprocess.PIPE,
                    text=True,
                )

            if p.returncode == 0 and tmp.exists() and tmp.stat().st_size > 0:
                # seqtk seq -l0 to flatten
                seq = subprocess.check_output(["seqtk", "seq", "-l0", str(tmp)], text=True)
                lines = [x for x in seq.splitlines() if x.strip()]
                if not lines or not lines[0].startswith(">"):
                    raise RuntimeError(f"{code}: seqtk produced invalid output for {r}")
                # Replace header
                header = f"{code}|{r}"
                seq_str = "".join(lines[1:]).strip()
                fasta_write_one(out_fa, header, seq_str)
                extracted += 1
            else:
                # Missing region -> pad Ns
                header = f"{code}|{r}"
                fasta_write_one(out_fa, header, "N" * L)
                missing_padded += 1
                with slog.open("a") as f:
                    f.write(f"[{now()}] MISSING padded: {r} (len={L})\n")

            if tmp.exists():
                tmp.unlink()

        # Per-sample concat across ALL regions (not filtered)
        concat = sample_dir / f"{code}_concat.fasta"
        if not concat.exists() or concat.stat().st_size == 0:
            seq_parts: List[str] = []
            for r in regions:
                _, s = fasta_read_one(sample_dir / f"{r}.fasta")
                seq_parts.append(s)
            fasta_write_one(concat, code, "".join(seq_parts))

        log(f"{code}: extracted={extracted} padded_missing={missing_padded}", runlog)

    log("Step 1 complete.", runlog)

    if not args.runtrees:
        log("Done (no --runtrees).", runlog)
        return

    # -----------------------------
    # Step 2: maxN filter (pre-TV masking)
    # -----------------------------
    filtering_dir = args.outdir / "filtering"
    filtering_dir.mkdir(parents=True, exist_ok=True)
    kept_maxn = filtering_dir / "regions.kept.maxN.txt"
    dropped_maxn = filtering_dir / "regions.dropped.maxN.txt"

    log(f"Step 2: Missingness filter (drop if ANY sample Nfrac > {args.maxN})", runlog)

    kept: List[str] = []
    with kept_maxn.open("w") as ok, dropped_maxn.open("w") as bad:
        for r in regions:
            drop_reason = None
            for code in codenames:
                f = args.outdir / code / f"{r}.fasta"
                if not f.exists():
                    drop_reason = (code, "missing_file")
                    break
                _, seq = fasta_read_one(f)
                frac = nfrac_seq(seq)
                if frac > args.maxN:
                    drop_reason = (code, f"nfrac={frac:.6f}")
                    break
            if drop_reason is None:
                ok.write(r + "\n")
                kept.append(r)
            else:
                bad.write(f"{r}\t{drop_reason[0]}\t{drop_reason[1]}\n")

    if not kept:
        raise SystemExit("ERROR: No windows left after maxN filtering.")
    log(f"maxN kept={len(kept)} dropped={len(regions)-len(kept)}", runlog)

    # -----------------------------
    # Step 3: Build per-window MSAs (kept after maxN)
    # -----------------------------
    comb_win_dir = args.outdir / "Combined_windows"
    comb_win_dir.mkdir(parents=True, exist_ok=True)

    log("Step 3: Build per-window multi-FASTA alignments (kept.maxN)", runlog)

    for r in kept:
        msa = comb_win_dir / f"{r}.fa"
        if msa.exists() and msa.stat().st_size > 0:
            continue
        with msa.open("w") as out:
            for code in codenames:
                out.write((args.outdir / code / f"{r}.fasta").read_text())

    # -----------------------------
    # Step 3b: TV-only masking (in place on MSAs)
    # -----------------------------
    if args.transversions:
        log("Step 3b: TV-only masking (conservative) on per-window MSAs (in place)", runlog)
        for r in kept:
            mask_transitions_inplace(comb_win_dir / f"{r}.fa")

    # -----------------------------
    # Step 3c: minPI filter (computed after TV masking if enabled)
    # -----------------------------
    kept_final = filtering_dir / "regions.kept.final.txt"
    dropped_pi = filtering_dir / "regions.dropped.minpi.txt"

    log(f"Step 3c: minPI filter on MSAs (keep if PI >= {args.minpi})", runlog)

    kept2: List[str] = []
    with kept_final.open("w") as ok, dropped_pi.open("w") as bad:
        for r in kept:
            msa = comb_win_dir / f"{r}.fa"
            headers, seqs = read_msa_fasta(msa)
            if not seqs:
                bad.write(f"{r}\tempty_msa\n")
                continue
            try:
                pi = count_parsimony_informative(seqs)
            except Exception as e:
                bad.write(f"{r}\tpi_error\t{e}\n")
                continue
            if pi >= args.minpi:
                ok.write(r + "\n")
                kept2.append(r)
            else:
                bad.write(f"{r}\tpi<{args.minpi}\tpi={pi}\n")

    if not kept2:
        raise SystemExit("ERROR: No windows left after minPI filtering.")
    log(f"minPI kept={len(kept2)} dropped={len(kept)-len(kept2)}", runlog)

    # -----------------------------
    # Step 4: Build concatenated alignment across FINAL kept windows
    #         If TV masking enabled, build from masked MSAs (consistent TV-only analysis)
    # -----------------------------
    comb_dir = args.outdir / "Combined"
    comb_dir.mkdir(parents=True, exist_ok=True)
    all_concat = comb_dir / ("All_concat.filtered.tv.fasta" if args.transversions else "All_concat.filtered.fasta")

    log(f"Step 4: Build concatenated alignment from FINAL windows: {all_concat.name}", runlog)

    with all_concat.open("w") as out:
        for code in codenames:
            out.write(f">{code}\n")
            for r in kept2:
                if args.transversions:
                    msa = comb_win_dir / f"{r}.fa"
                    seq = None
                    with msa.open() as f:
                        take = False
                        buf: List[str] = []
                        for line in f:
                            line = line.strip()
                            if not line:
                                continue
                            if line.startswith(">"):
                                if take:
                                    seq = "".join(buf)
                                    break
                                take = (line[1:].split()[0] == code)
                                buf = []
                            else:
                                if take:
                                    buf.append(line)
                        if seq is None and take:
                            seq = "".join(buf)
                    if seq is None:
                        raise RuntimeError(f"Could not find taxon {code} in {msa}")
                    out.write(seq + "\n")
                else:
                    _, seq = fasta_read_one(args.outdir / code / f"{r}.fasta")
                    out.write(seq + "\n")
            out.write("\n")

    # -----------------------------
    # Step 5: Window trees with IQ-TREE2 (parallel across windows)
    # -----------------------------
    tree_dir = args.outdir / "trees"
    win_tree_dir = tree_dir / "window_trees"
    win_log_dir = tree_dir / "window_logs"
    ref_dir = tree_dir / "reference"
    cf_dir = tree_dir / "concordance"
    for d in [tree_dir, win_tree_dir, win_log_dir, ref_dir, cf_dir]:
        d.mkdir(parents=True, exist_ok=True)

    log("Step 5: Infer window trees with IQ-TREE2 (parallel across windows)", runlog)

    jobs: List[WindowJob] = []
    for r in kept2:
        msa = comb_win_dir / f"{r}.fa"
        safe = sanitize_region(r)
        out_prefix = win_tree_dir / safe
        iqlog = win_log_dir / f"{safe}.iqtree.log"
        jobs.append(
            WindowJob(
                region=r,
                msa_path=msa,
                out_prefix=out_prefix,
                iqtree_threads=T,
                model=args.model,
                bootstrap=args.bootstrap,
                log_path=iqlog,
            )
        )

    ok_count = 0
    fail_count = 0
    failed: List[Tuple[str, str]] = []

    with cf.ThreadPoolExecutor(max_workers=W) as ex:
        futs = [ex.submit(iqtree_window, j) for j in jobs]
        for fut in cf.as_completed(futs):
            region, success, msg = fut.result()
            if success:
                ok_count += 1
            else:
                fail_count += 1
                failed.append((region, msg))

    log(f"Window trees complete: ok={ok_count}, failed={fail_count}", runlog)
    if fail_count > 0:
        fail_path = tree_dir / "failed_windows.txt"
        with fail_path.open("w") as f:
            for r, m in failed:
                f.write(f"{r}\t{m}\n")
        raise SystemExit(f"ERROR: {fail_count} window trees failed. See: {fail_path}")

    # Collect gene trees (deterministic order = kept2)
    all_trs = tree_dir / "all_window_trees.trs"
    log(f"Writing gene-tree set: {all_trs}", runlog)
    with all_trs.open("w") as out:
        for r in kept2:
            safe = sanitize_region(r)
            treefile = (win_tree_dir / safe).with_suffix(".treefile")
            out.write(treefile.read_text().strip() + "\n")

    # -----------------------------
    # Step 6: Reference tree (concat or ASTRAL)
    # -----------------------------
    ref_tree: Path

    if args.ref == "concat":
        log("Step 6: Build concatenated reference tree (IQ-TREE2)", runlog)
        pref = ref_dir / "Reference"
        ref_log = ref_dir / "Reference.iqtree.log"
        # Use modest threads for this single run (avoid oversubscribing)
        ref_threads = max(1, min(C, 10))
        run_cmd(
            [
                "iqtree2",
                "-s", str(all_concat),
                "-m", args.model,
                "-bb", str(args.bootstrap),
                "-nt", str(ref_threads),
                "-pre", str(pref),
            ],
            ref_log,
        )
        ref_tree = pref.with_suffix(".treefile")
    else:
        log("Step 6: Build ASTRAL reference tree (collapse then ASTRAL)", runlog)

        # Biopython is guaranteed by preflight when needed
        from io import StringIO
        from Bio import Phylo  # type: ignore

        collapsed = tree_dir / "all_window_trees.collapsed.trs"
        thr = float(args.collapse)

        def support_from_label(x: Optional[str]) -> Optional[float]:
            if x is None:
                return None
            s = str(x).strip()
            m = re.match(r"^([0-9]+(?:\.[0-9]+)?)", s)
            if not m:
                return None
            try:
                return float(m.group(1))
            except Exception:
                return None

        def get_support(clade) -> Optional[float]:
            sup = support_from_label(getattr(clade, "confidence", None))
            if sup is not None:
                return sup
            return support_from_label(getattr(clade, "name", None))

        def collapse_low_support(tree, thr: float):
            changed = True
            while changed:
                changed = False
                for clade in list(tree.find_clades(order="preorder")):
                    if clade == tree.root or clade.is_terminal():
                        continue
                    sup = get_support(clade)
                    if sup is None:
                        continue
                    if sup < thr:
                        tree.collapse(clade)
                        changed = True
                        break
            return tree

        with all_trs.open() as inp, collapsed.open("w") as out:
            for line in inp:
                line = line.strip()
                if not line:
                    continue
                t = Phylo.read(StringIO(line), "newick")
                t = collapse_low_support(t, thr)
                buf = StringIO()
                Phylo.write(t, buf, "newick")
                out.write(buf.getvalue().strip() + "\n")

        if args.astral is None:
            raise SystemExit("ERROR: --astral is required for --ref astral")
        astral_out = ref_dir / "astral_species_tree.nwk"
        astral_log = ref_dir / "astral.log"
        run_cmd(["java", "-jar", str(args.astral), "-i", str(collapsed), "-o", str(astral_out)], astral_log)
        ref_tree = astral_out

    log(f"Reference tree: {ref_tree}", runlog)

    # -----------------------------
    # Step 7: Concordance factors (gCF/sCF)
    # -----------------------------
    log("Step 7: Compute gCF/sCF (IQ-TREE2)", runlog)
    cf_pref = cf_dir / "Finaltrees"
    cf_log = cf_dir / "Finaltrees.cf.log"
    cf_threads = max(1, min(C, 10))
    run_cmd(
        [
            "iqtree2",
            "-t", str(ref_tree),
            "--gcf", str(all_trs),
            "-s", str(all_concat),
            "--scf", str(args.scf),
            "-nt", str(cf_threads),
            "-pre", str(cf_pref),
        ],
        cf_log,
    )

    # -----------------------------
    # Step 8: Optional topology filter -> dating supermatrix
    # -----------------------------
    if args.topofilter:
        log("Step 8: Topology filter + min bootstrap + dating supermatrix", runlog)

        from io import StringIO
        from Bio import Phylo  # type: ignore

        def read_tree(path: Path):
            s = path.read_text().strip()
            return Phylo.read(StringIO(s), "newick")

        def tips(tree):
            return {t.name for t in tree.get_terminals()}

        def clade_tips(clade):
            return {t.name for t in clade.get_terminals()}

        def unrooted_splits(tree):
            Tset = tips(tree)
            n = len(Tset)
            splits = set()
            for clade in tree.get_nonterminals(order="preorder"):
                s = clade_tips(clade)
                k = len(s)
                if k <= 1 or k >= n - 1:
                    continue
                if k > n - k:
                    s = Tset - s
                splits.add(frozenset(s))
            return splits, Tset

        def support_from_label(x: Optional[str]) -> Optional[float]:
            if x is None:
                return None
            s = str(x).strip()
            m = re.match(r"^([0-9]+(?:\.[0-9]+)?)", s)
            if not m:
                return None
            try:
                return float(m.group(1))
            except Exception:
                return None

        def get_support(clade) -> Optional[float]:
            sup = support_from_label(getattr(clade, "confidence", None))
            if sup is not None:
                return sup
            return support_from_label(getattr(clade, "name", None))

        def min_internal_support(tree) -> Optional[float]:
            Tset = tips(tree)
            n = len(Tset)
            vals = []
            for clade in tree.get_nonterminals(order="preorder"):
                if clade == tree.root:
                    continue
                s = clade_tips(clade)
                k = len(s)
                if k <= 1 or k >= n - 1:
                    continue
                sup = get_support(clade)
                if sup is None:
                    return None
                vals.append(sup)
            return min(vals) if vals else None

        ref = read_tree(ref_tree)
        ref_splits, ref_tips = unrooted_splits(ref)

        suffix = "tv" if args.transversions else "full"
        topomatch = filtering_dir / f"regions.topomatch.minbs{args.minbs}.{suffix}.txt"
        topofail = filtering_dir / f"regions.topofail.minbs{args.minbs}.{suffix}.txt"

        kept_topo: List[str] = []
        with topomatch.open("w") as ok, topofail.open("w") as bad:
            for r in kept2:
                safe = sanitize_region(r)
                tf = (win_tree_dir / safe).with_suffix(".treefile")
                if not tf.exists():
                    bad.write(f"{r}\tmissing_treefile\n")
                    continue
                gt = read_tree(tf)
                gt_splits, gt_tips = unrooted_splits(gt)
                if gt_tips != ref_tips:
                    bad.write(f"{r}\ttaxa_mismatch\n")
                    continue
                mbs = min_internal_support(gt)
                if mbs is None:
                    bad.write(f"{r}\tmissing_support\n")
                    continue
                if mbs < args.minbs:
                    bad.write(f"{r}\tmin_bootstrap<{args.minbs}\tmin={mbs}\n")
                    continue

                if args.topomode == "exact":
                    ok_topo = (gt_splits == ref_splits)
                else:
                    ok_topo = gt_splits.issubset(ref_splits)

                if ok_topo:
                    ok.write(r + "\n")
                    kept_topo.append(r)
                else:
                    bad.write(f"{r}\ttopology_mismatch\n")

        if not kept_topo:
            raise SystemExit("ERROR: No windows passed topology+minBS filter.")

        dating = comb_dir / f"All_concat.topomatch.minbs{args.minbs}.{suffix}.fasta"
        log(f"Writing dating supermatrix: {dating}", runlog)

        with dating.open("w") as out:
            for code in codenames:
                out.write(f">{code}\n")
                for r in kept_topo:
                    if args.transversions:
                        msa = comb_win_dir / f"{r}.fa"
                        seq = None
                        with msa.open() as f:
                            take = False
                            buf: List[str] = []
                            for line in f:
                                line = line.strip()
                                if not line:
                                    continue
                                if line.startswith(">"):
                                    if take:
                                        seq = "".join(buf)
                                        break
                                    take = (line[1:].split()[0] == code)
                                    buf = []
                                else:
                                    if take:
                                        buf.append(line)
                            if seq is None and take:
                                seq = "".join(buf)
                        if seq is None:
                            raise RuntimeError(f"Could not find taxon {code} in {msa}")
                        out.write(seq + "\n")
                    else:
                        _, seq = fasta_read_one(args.outdir / code / f"{r}.fasta")
                        out.write(seq + "\n")
                out.write("\n")

    log("Done.", runlog)
    log("Key outputs:", runlog)
    log(f"  regions used: {args.regions}", runlog)
    log(f"  kept.final:   {filtering_dir / 'regions.kept.final.txt'}", runlog)
    log(f"  concat:       {all_concat}", runlog)
    log(f"  gene trees:   {all_trs}", runlog)
    log(f"  ref tree:     {ref_tree}", runlog)
    log(f"  CF prefix:    {cf_dir / 'Finaltrees'}", runlog)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        raise SystemExit("Interrupted by user.")
