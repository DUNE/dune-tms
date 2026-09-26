#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_validation.sh [--cluster3d] <input file|filelist.txt|glob|directory> [output name|output.root] [num_events] [max_slices]
  ./run_validation.sh <input file|filelist.txt|glob|directory> [num_events] [max_slices]
  ./run_validation.sh <input file|filelist.txt|glob|directory> [num_events] [max_slices] [output name|output.root]

If the second argument is non-numeric, it is used as the output name/stem.
Directory inputs are scanned recursively by Tracking_Validation.
Text-file inputs are treated as one ROOT file path or URL per line.

--cluster3d validates the Cluster3D reconstruction instead of the legacy one:
it reads Reco_Tree_C3D / Truth_Info_C3D (written by ConvertToTMSTree with
[Recon.Cluster3D] Enabled). Equivalent to setting TMS_VALIDATION_RECO_TREE and
TMS_VALIDATION_TRUTH_TREE.
EOF
}

is_integer() {
  [[ "$1" =~ ^-?[0-9]+$ ]]
}

validation_dir="/exp/dune/data/users/${USER}/dune-tms/Validation/Tracking_Validation"

if [[ $# -gt 0 ]] && [[ "$1" == "--cluster3d" ]]; then
  export TMS_VALIDATION_RECO_TREE=Reco_Tree_C3D
  export TMS_VALIDATION_TRUTH_TREE=Truth_Info_C3D
  shift
fi

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
  output_name="${output_name%.txt}"
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
