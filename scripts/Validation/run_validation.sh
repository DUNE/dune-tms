#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_validation.sh <input file|glob|directory> [output name|output.root] [num_events] [max_slices]
  ./run_validation.sh <input file|glob|directory> [num_events] [max_slices]
  ./run_validation.sh <input file|glob|directory> [num_events] [max_slices] [output name|output.root]

If the second argument is non-numeric, it is used as the output name/stem.
Directory inputs are scanned recursively by Tracking_Validation.
EOF
}

is_integer() {
  [[ "$1" =~ ^-?[0-9]+$ ]]
}

validation_dir="/exp/dune/data/users/${USER}/dune-tms/Validation/Tracking_Validation"

if [[ $# -lt 1 ]]; then
  usage >&2
  exit 1
fi

infile="$1"
shift

output_name=""
if [[ $# -gt 0 ]] && ! is_integer "$1"; then
  output_name="$1"
  shift
fi

num_events="${1:--1}"
if [[ $# -gt 0 ]]; then
  shift
fi

max_slices="${1:--1}"
if [[ $# -gt 0 ]]; then
  shift
fi

if [[ $# -gt 0 ]] && [[ -z "$output_name" ]]; then
  output_name="$1"
  shift
fi

if [[ $# -gt 0 ]]; then
  echo "Unexpected extra arguments: $*" >&2
  usage >&2
  exit 1
fi

if [[ -z "$output_name" ]]; then
  base_filename=$(basename "$infile")
  output_name="${base_filename%.root}"
fi

if [[ "$output_name" == */* ]]; then
  outfile="$output_name"
else
  outfile="${validation_dir}/${output_name}"
fi

if [[ "$outfile" != *.root ]]; then
  outfile="${outfile}.root"
fi

mkdir -p "$(dirname "$outfile")"

make
./Tracking_Validation "$infile" "$num_events" "$max_slices" "$outfile"
python simply_draw_everything.py "$outfile"

echo "Output should be in"
echo "${outfile/.root/_images}"
