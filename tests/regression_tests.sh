#!/usr/bin/env bash

work_dir=tests

SUDO_OPT=
while getopts ":s" opt; do
  case ${opt} in
    s ) SUDO_OPT="sudo "
      ;;
    \? ) echo "Usage: $0 [-s to run with sudo]"
         exit 1
      ;;
  esac
done

echo "Starting regression tests..."
echo ""

# Run pipeline on test input data
${SUDO_OPT}nextflow -log nextflow_test.log run palidis.nf --manifest ${work_dir}/test_manifest.txt --batch_name ${work_dir}/regression_test_data/output -profile github_ci

cat nextflow_test.log
echo ""

function file_diff {
    local test=${1}
    local reference=${2}

    diff ${test} ${reference} > /dev/null 2>&1
    status=$?

    # If files different then warning
    if [[ ${status} -gt 0 ]]; then
        if [[ ${status} -eq 1 ]]; then
            echo ""
            echo "The contents of ${test} is not expected."
        elif [[ ${status} -eq 2 ]]; then
            echo ""
            echo "Unable to perform differences check."
            if [[ ! -f $test ]]; then
                echo "Expected output ${test} is missing."
            fi
        fi
        return 1
    fi
}

error_status=0

file_diff "${work_dir}/regression_test_data/output/test_insertion_sequences_info.txt" "${work_dir}/regression_test_data/output/reference_insertion_sequences_info.txt"
out=$?
error_status=$(($error_status | $out))

file_diff "${work_dir}/regression_test_data/output/test_insertion_sequences.fasta" "${work_dir}/regression_test_data/output/reference_insertion_sequences.fasta"
out=$?
error_status=$(($error_status | $out))

# Error if any output files missing or not expected
if [[ ${error_status} -eq 1 ]]; then
    echo ""
    echo "Test failed. Outputs listed may be missing or their contents not expected."
    exit 1
else
    echo "Test passed. All outputs expected."
fi
