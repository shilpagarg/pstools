#!/bin/bash

ROOT_DIR="$(dirname "${0}")/.."
pushd "${ROOT_DIR}"
ROOT_DIR=$(pwd)
popd

TARGET="${ROOT_DIR}/pstools"
# TEST_COMMANDS="qv resolve_haplotypes"
TEST_COMMANDS="qv"
# TEST_COMMANDS="resolve_haplotypes"
# Debug option.
TEST_DEBUG=${TEST_DEBUG-}
GDB_ARGS=
if [ -n "${TEST_DEBUG}" ]; then
    GDB_ARGS='gdb -q -ex "set breakpoint pending on" --args'
fi

set -eux

test_qv() {
    # The command returns non-zero even when the command is success.
    "${TARGET}" qv ./true_hap2.fasta ./hap2.1.fastq ./hap2.2.fastq > qv_score.txt || :
    grep QV_RAW -A 1 qv_score.txt
    grep -E '^19\.628\s+-1\.000' qv_score.txt
}

test_resolve_haplotypes() {
    rm -rf out
    ${GDB_ARGS} "${TARGET}" resolve_haplotypes -t 8 ./hic_name_connection.output ./icHarAxyr1.r_utg.gfa out
}

"${TARGET}" 2>&1 | grep -q Usage

for cmd in ${TEST_COMMANDS}; do
    pushd test/data/${cmd}
    test_${cmd}
    popd
done
