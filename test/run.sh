#!/bin/bash

ROOT_DIR="$(dirname "${0}")/.."
pushd "${ROOT_DIR}"
ROOT_DIR=$(pwd)
popd

TARGET="${ROOT_DIR}/pstools"

set -eux

"${TARGET}" 2>&1 | grep -q Usage

pushd test/data/qv
# The command returns non-zero even when the command is success.
"${TARGET}" qv ./true_hap2.fasta ./hap2.1.fastq ./hap2.2.fastq > qv_score.txt || :
grep QV_RAW -A 1 qv_score.txt
grep -E '^19\.628\s+-1\.000' qv_score.txt
popd
