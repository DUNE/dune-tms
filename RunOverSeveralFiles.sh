#!/bin/bash

# Directory containing the files
mydir="/pnfs/dune/persistent/users/abooth/nd-production/MicroProdN1p2/output/run-spill-build/MicroProdN1p2_NDLAr_1E18_RHC.spill.nu/EDEPSIM_SPILLS/0000000/0000000"
# Directory to move FooBar* output
outdir="/exp/dune/data/users/${USER}/dune-tms/2024-09-19_test_liam_lmao"
# Directory to store logs
logdir="${outdir}/logs"

NFILES=20

# Create log directory if it doesn't exist
mkdir -p "$logdir"

# Initialize counters
nonzero_exit_count=0
total_time=0
num_files=0

# Loop over files in the directory found by 'find'
for file in $(find "$mydir" -type f | sort | head -n ${NFILES}); do
  # Get the base name of the file (without path)
  basefile=$(basename "$file")
  logfile="$logdir/$basefile.log"

  # Run the command and capture the start time
  start_time=$(date +%s)
  echo "Running $file"
  timeout 5m stdbuf -oL -eL ConvertToTMSTree.exe "$file" >"$logfile" 2>&1
  exit_code=$?
  end_time=$(date +%s)

  # Calculate the time taken for this command
  elapsed_time=$((end_time - start_time))
  total_time=$((total_time + elapsed_time))
  num_files=$((num_files + 1))

  # Check the exit code
  if [ $exit_code -ne 0 ]; then
    echo "Nonzero exit code ($exit_code) for file: $file"
    nonzero_exit_count=$((nonzero_exit_count + 1))
    rm -f MicroProdN1p2*
  else
    echo "Got exit code ($exit_code) for file: $file"
    # Move FooBar* output files to mydir2 if they exist
    mv MicroProdN1p2* "$outdir" 2>/dev/null
  fi
done

# Print the summary
echo "Total files processed: $num_files"
echo "Nonzero exit codes: $nonzero_exit_count"
if [ $num_files -gt 0 ]; then
  average_time=$((total_time / num_files))
  echo "Average time per command: $average_time seconds"
fi
