#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Sliding-window FASTA extraction + concatenation (multi-sample)
#
# REQUIRED ARGUMENTS
#   --input   samples.tsv   (CODENAME<TAB>/path/to/genome.fa)
#   --regions regions.txt   (chr:start-end, 1-based, inclusive)
#   --outdir  OUTDIR
#
# OUTPUTS (per sample CODENAME)
#   OUTDIR/CODENAME/<region>.fasta        header: >CODENAME|region
#   OUTDIR/CODENAME/CODENAME_concat.fasta header: >CODENAME
#   OUTDIR/CODENAME/CODENAME.log
#   OUTDIR/run.log
#
# REQUIREMENTS
#   samtools, seqtk
# ============================================================

usage() {
  cat <<EOF
Usage:
  $0 --input samples.tsv --regions regions.txt --outdir OUTDIR

Required arguments:
  --input     Sample table (CODENAME<TAB>/path/to/genome.fa)
  --regions   Regions file (chr:start-end, one per line)
  --outdir    Output directory

Notes:
  - Regions must be explicit faidx intervals (chr:start-end)
  - Missing regions are padded with Ns of the correct length
EOF
}

# ----------------------------
# Parse arguments
# ----------------------------
INPUT=""
REGIONS=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      INPUT="$2"
      shift 2
      ;;
    --regions)
      REGIONS="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# ----------------------------
# Validate arguments
# ----------------------------
if [[ -z "$INPUT" || -z "$REGIONS" || -z "$OUTDIR" ]]; then
  echo "ERROR: --input, --regions, and --outdir are required" >&2
  usage
  exit 1
fi

[[ -f "$INPUT" ]]   || { echo "ERROR: input file not found: $INPUT" >&2; exit 1; }
[[ -f "$REGIONS" ]] || { echo "ERROR: regions file not found: $REGIONS" >&2; exit 1; }

mkdir -p "$OUTDIR"
RUNLOG="$OUTDIR/run.log"

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

log_run() {
  echo "[$(timestamp)] $*" | tee -a "$RUNLOG"
}

log_sample() {
  echo "[$(timestamp)] $*" | tee -a "$SAMPLE_LOG" >/dev/null
}

# Compute expected length from chr:start-end
region_len() {
  local region="$1"
  if [[ "$region" =~ ^[^:]+:([0-9]+)-([0-9]+)$ ]]; then
    local start="${BASH_REMATCH[1]}"
    local end="${BASH_REMATCH[2]}"
    (( end >= start )) && echo $(( end - start + 1 )) || echo 0
  else
    echo 0
  fi
}

emit_ns() {
  local L="$1"
  awk -v n="$L" 'BEGIN{ s=""; for(i=1;i<=n;i++) s=s "N"; print s }'
}

# ----------------------------
# Run
# ----------------------------
log_run "Starting run"
log_run "Input samples : $INPUT"
log_run "Regions file : $REGIONS"
log_run "Output dir   : $OUTDIR"
log_run "----------------------------------------"

REGION_COUNT="$(grep -cve '^\s*$' "$REGIONS" || true)"
log_run "Regions (non-empty): $REGION_COUNT"

# Process samples table (TAB-separated; lines starting with # ignored)
while IFS=$'\t' read -r CODENAME GENOME; do
  [[ -z "${CODENAME:-}" ]] && continue
  [[ "${CODENAME:0:1}" == "#" ]] && continue

  if [[ -z "${GENOME:-}" ]]; then
    log_run "ERROR: Missing genome path for codename '$CODENAME'"
    exit 1
  fi
  if [[ ! -f "$GENOME" ]]; then
    log_run "ERROR: Genome FASTA not found for $CODENAME: $GENOME"
    exit 1
  fi

  SAMPLE_DIR="$OUTDIR/$CODENAME"
  mkdir -p "$SAMPLE_DIR"
  SAMPLE_LOG="$SAMPLE_DIR/${CODENAME}.log"

  log_run "==> Processing $CODENAME"
  log_sample "Genome: $GENOME"
  log_sample "Regions file: $REGIONS"

  # Index genome if needed
  if [[ ! -f "${GENOME}.fai" ]]; then
    log_sample "Indexing genome with samtools faidx"
    samtools faidx "$GENOME"
  fi

  extracted=0
  padded=0

  while IFS= read -r region; do
    [[ -z "$region" ]] && continue

    # Enforce explicit interval format
    if [[ ! "$region" =~ ^[^:]+:[0-9]+-[0-9]+$ ]]; then
      log_sample "ERROR: Invalid region format (must be chr:start-end): $region"
      exit 1
    fi

    out="$SAMPLE_DIR/${region}.fasta"
    L="$(region_len "$region")"

    if samtools faidx "$GENOME" "$region" > "$SAMPLE_DIR/.tmp.faidx" 2>>"$SAMPLE_LOG"; then
      seqtk seq -l0 "$SAMPLE_DIR/.tmp.faidx" \
        | sed "1s/^>.*/>${CODENAME}|${region}/" \
        > "$out"
      extracted=$((extracted + 1))
    else
      padded=$((padded + 1))
      log_sample "MISSING: $region (padded with N length=$L)"
      {
        echo ">${CODENAME}|${region}"
        emit_ns "$L"
      } > "$out"
    fi
  done < "$REGIONS"

  rm -f "$SAMPLE_DIR/.tmp.faidx"

  log_sample "Extracted windows: $extracted"
  log_sample "Missing windows padded: $padded"

  CONCAT="$SAMPLE_DIR/${CODENAME}_concat.fasta"
  {
    echo ">${CODENAME}"
    while IFS= read -r region; do
      [[ -z "$region" ]] && continue
      grep -v '^>' "$SAMPLE_DIR/${region}.fasta"
    done < "$REGIONS"
    echo
  } > "$CONCAT"

  concat_len="$(grep -v '^>' "$CONCAT" | awk '{n+=length($0)} END{print n+0}')"
  log_sample "Wrote concat FASTA: $CONCAT"
  log_sample "Concat length (bp): $concat_len"

  log_run "<== Done $CODENAME (padded windows: $padded)"

done < "$INPUT"

log_run "Finished run"
