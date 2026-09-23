#!/usr/bin/env bash
#
# run_grid.sh
#
# Setup: ONE reference index (.mmi) and THREE ONT fasta query files
# (e.g. complete genome, locus assembly, ground truth) that all get
# aligned to that same reference.
#
#   PART A: cross-file comparison
#     Same reference, same baseline params, different query fasta.
#     Answers: "does minimap2 behave differently depending on which
#     fasta file I align?"
#
#   PART B: per-file parameter sweep
#     For EACH query fasta separately, run several mapping-time parameter
#     variants and compare them against each other.
#     Answers: "which parameters are good for this particular file?"
#
# NOTE: -k/-w/-H are frozen into the .mmi at index time (whoever built it).
# This script only varies mapping-time flags (-N, --secondary, -A/-B/-O/-E,
# etc.) since you're aligning against an existing .mmi. If you need to sweep
# -k/-w themselves, you'd have to rebuild the index from its source FASTA
# per value -- not covered here since you only have the .mmi in hand.
#
# Usage:
#   ./run_grid.sh reference.mmi outdir \
#       complete_genome=/path/complete_genome.fasta \
#       locus_assembly=/path/locus_assembly.fasta \
#       ground_truth=/path/ground_truth.fasta
#
set -euo pipefail

REF="$1"; shift
OUTDIR="$1"; shift
QUERIES=("$@")   # remaining args: name=path.fasta

THREADS="${THREADS:-8}"
mkdir -p "$OUTDIR"

# Mapping-time parameter variants to test within each query file.
# name|extra_minimap2_args
VARIANTS=(
  "baseline|"
  "no_secondary|--secondary=no"
  "high_N|-N 20"
  "relaxed_scoring|-A2 -B4 -O4,24 -E2,1"
)

run_one() {
  local query_name="$1" query_path="$2" variant_name="$3" variant_args="$4"
  local bam="$OUTDIR/${query_name}__${variant_name}.sorted.bam"
  if [[ -f "$bam" ]]; then
    echo "  [skip, exists] $bam"
    return
  fi
  echo "  aligning: query=$query_name variant=$variant_name args='$variant_args'"
  # shellcheck disable=SC2086
  minimap2 -t "$THREADS" -ax map-ont $variant_args "$REF" "$query_path" \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"
}

echo "=============================================="
echo "Aligning each query fasta x each variant, to reference: $REF"
echo "=============================================="
# Parallel arrays instead of an associative array (declare -A needs bash 4+,
# which macOS doesn't ship by default -- this works on bash 3.2 too).
QUERY_NAMES=()
QUERY_PATHS=()
for entry in "${QUERIES[@]}"; do
  name="${entry%%=*}"
  path="${entry#*=}"
  QUERY_NAMES+=("$name")
  QUERY_PATHS+=("$path")
  echo "-- query: $name ($path) --"
  for v in "${VARIANTS[@]}"; do
    IFS="|" read -r vname vargs <<< "$v"
    run_one "$name" "$path" "$vname" "$vargs"
  done
done

echo
echo "=============================================="
echo "PART A: cross-file comparison (baseline params only)"
echo "=============================================="
baseline_bams=()
baseline_labels=()
for name in "${QUERY_NAMES[@]}"; do
  baseline_bams+=("$OUTDIR/${name}__baseline.sorted.bam")
  baseline_labels+=("$name")
done
python3 validate_minimap2.py "${baseline_bams[@]}" \
  --preset map-ont \
  --labels "${baseline_labels[@]}" \
  --csv "$OUTDIR/cross_file_comparison.csv" \
  | tee "$OUTDIR/cross_file_comparison.txt"

echo
echo "=============================================="
echo "PART B: per-file parameter sweep"
echo "=============================================="
for name in "${QUERY_NAMES[@]}"; do
  echo
  echo "---- query: $name ----"
  variant_bams=()
  variant_labels=()
  for v in "${VARIANTS[@]}"; do
    IFS="|" read -r vname vargs <<< "$v"
    variant_bams+=("$OUTDIR/${name}__${vname}.sorted.bam")
    variant_labels+=("$vname")
  done
  python3 validate_minimap2.py "${variant_bams[@]}" \
    --preset map-ont \
    --labels "${variant_labels[@]}" \
    --csv "$OUTDIR/${name}__param_sweep.csv" \
    | tee "$OUTDIR/${name}__param_sweep.txt"
done

echo
echo "Done. Outputs in $OUTDIR/:"
echo "  cross_file_comparison.{txt,csv}    -- does the query fasta matter?"
echo "  <query_name>__param_sweep.{txt,csv} -- best params per query file"