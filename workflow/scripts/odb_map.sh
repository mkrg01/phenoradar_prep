#!/usr/bin/env bash
# One local ODB job, invoked inside a resource allocation by run_odb_chunk.py.
set -Eeuo pipefail
ulimit -s unlimited
command=$1
manifest=$2
label=$3
node=$4
version=$5
data=$6
work=$7
jobs=$8
batch=$9
export ODBMAPPER_WORK="$work/odbmapper"
project="$ODBMAPPER_WORK/$version/pipeline"
mkdir -p "$work/tmp"
export TMPDIR="$work/tmp"

if [[ ! -e $ODBMAPPER_WORK/$version/data && ! -L $ODBMAPPER_WORK/$version/data ]]; then
    "$command" SETUP
fi
if [[ -L $ODBMAPPER_WORK/$version/data ]]; then
    [[ $(readlink -f "$ODBMAPPER_WORK/$version/data") == "$data" ]]
else
    rmdir "$ODBMAPPER_WORK/$version/data/tarfiles"
    rmdir "$ODBMAPPER_WORK/$version/data"
    ln -s "$data" "$ODBMAPPER_WORK/$version/data"
fi
"$command" SETUP
[[ $("$command" CONFIG project) == "$project" ]]
config="$project/orthologer_conf.sh"
set_scalar() {
    grep -q "^${1}=" "$config"
    sed -i "s|^${1}=.*|${1}=${2}|" "$config"
}
set_scalar SCHEDULER_LABEL NONE
set_scalar OP_NJOBMAX_BATCH "$batch"
set_scalar OP_NJOBMAX_LOCAL "$jobs"
set_scalar OP_SAVE_JOBLOG 0
set_scalar SKIP_REMAKE_CHECK 0
set_scalar STEP_SLEEP 0
set_scalar TMP_DIR_BASE "\"$work/tmp\""
for step in PREPROC MASKER SELECT STATS FORMATDB ALIGNMENT MAKEBRH MAKEINPAR MAKEINPARSEL; do
    grep -Fq "OP_STEP_NPARALLEL[$step]=" "$config"
    sed -i "s|^OP_STEP_NPARALLEL\[$step\]=.*|OP_STEP_NPARALLEL[$step]=$jobs|" "$config"
done
"$command" MAP "$label" "$manifest" "$node"
"$command" REPORT "$label" > "$work/report.txt"
"$command" CONFIG > "$work/odbmapper_config.txt"
