#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Script Name: run_make_heatmap.sh
# Author: Sakuntha Devaka Gunarathna
# Date: 2025-07-31
# Description: 
#   Batch processing script to run "make_heatmap" on all .singleFrag.bed files 
#   in the input directory. Results and logs are generated per sample. 
#
# Usage:
#   ./run_make_heatmap.sh <input_directory> <reference_file> <output_prefix>
#
# Example:
#   ./run_make_heatmap.sh ./cfDNA_samples reference_regions.txt ATAC_10kbp
# ----------------------------------------------------------------------

# ------------------ ARGUMENTS ------------------
INPUT_DIR="${1:-./}"                  # Default: current directory
REFERENCE_FILE="${2:-reference_regions.txt}"   # Default reference file
OUTPUT_PREFIX="${3:-Results}"         # Default prefix for outputs
MAX_CONCURRENT_JOBS=5                 # Maximum number of concurrent jobs
THREADS=10                            # Number of threads per job
TOOL_PATH="./make_heatmap"            # Path to make_heatmap binary (in repo)
WINDOW_SIZE=5000                      # Upstream/downstream window size
STEP_SIZE=20                          # Step size (bp)
BIN_SIZE=500                          # Bin size (bp)
# ------------------------------------------------

# Validate reference file
if [[ ! -f "$REFERENCE_FILE" ]]; then
  echo "ERROR: Reference file '$REFERENCE_FILE' not found."
  exit 1
fi

# Create or clear logs
> successful_samples.log
> failed_samples.log

job_count=0

# Loop through all ".singleFrag.bed" files
for file in "${INPUT_DIR}"/*.singleFrag.bed; do
  [[ -e "$file" ]] || { echo "No .singleFrag.bed files found in $INPUT_DIR"; exit 1; }

  sample_name=$(basename "${file}" .singleFrag.bed)
  output_file="${sample_name}_${OUTPUT_PREFIX}.txt"
  log_file="${sample_name}_${OUTPUT_PREFIX}.log"

  echo "Starting ${sample_name}..."

  (
    "${TOOL_PATH}"       -t "${THREADS}" -b c -h b -l c -a u -v t -d p -s b --       "${file}" "${REFERENCE_FILE}" "${output_file}"       -${WINDOW_SIZE} "${STEP_SIZE}" "${BIN_SIZE}" > "${log_file}" 2>&1

    if [[ $? -eq 0 ]]; then
      echo "${sample_name}" >> successful_samples.log
    else
      echo "${sample_name}" >> failed_samples.log
      echo "FAILED: ${sample_name}. See ${log_file} for details."
    fi
  ) &

  ((job_count++))

  if [[ "$job_count" -ge "$MAX_CONCURRENT_JOBS" ]]; then
    wait
    job_count=0
  fi
done

# Wait for any remaining jobs
wait
echo "All jobs completed. See 'successful_samples.log' and 'failed_samples.log'."
