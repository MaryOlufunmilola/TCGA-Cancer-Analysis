#!/usr/bin/env bash
# Run the unit tests, the toy example, and both analysis scripts in order,
# logging each step to logs/ and recording its run time and peak memory.
#
# Usage (from the repo root, inside the tcga-cancer-analysis environment):
#   bash run_all.sh               # tests, toy example, analysis.R, mutation_analysis.R
#   bash run_all.sh --report      # also render report.Rmd to report.html
#   bash run_all.sh --skip-tests  # skip the unit tests and toy example
#
# Run time and peak memory for each step go to logs/timing.tsv and are printed
# at the end. analysis.R's time includes the GDC download on the first run.

set -euo pipefail

RUN_TESTS=true
RUN_REPORT=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --report) RUN_REPORT=true; shift ;;
    --skip-tests) RUN_TESTS=false; shift ;;
    -h|--help) sed -n '2,11p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

cd "$(dirname "$0")"
mkdir -p logs

TIMING="logs/timing.tsv"
printf "step\tseconds\tpeak_memory_mb\n" > "$TIMING"
TOTAL_START=$SECONDS

# Peak memory comes from /usr/bin/time: -l on macOS (bytes), -v on Linux (kilobytes).
if [[ "$(uname)" == "Darwin" ]]; then TIME_FLAG="-l"; else TIME_FLAG="-v"; fi
if [[ ! -x /usr/bin/time ]]; then TIME_FLAG=""; fi

run_step() {
  local name="$1"; shift
  local log="logs/${name}.log"
  echo ""
  echo "==> ${name}"
  local start=$SECONDS
  if [[ -n "$TIME_FLAG" ]]; then
    /usr/bin/time "$TIME_FLAG" "$@" 2>&1 | tee "$log"
  else
    "$@" 2>&1 | tee "$log"
  fi
  local secs=$(( SECONDS - start ))
  local mem="NA"
  if [[ "$TIME_FLAG" == "-l" ]]; then
    mem=$(awk '/maximum resident set size/ {printf "%.0f", $1/1048576}' "$log")
  elif [[ "$TIME_FLAG" == "-v" ]]; then
    mem=$(awk -F: '/Maximum resident set size/ {printf "%.0f", $2/1024}' "$log")
  fi
  printf "%s\t%s\t%s\n" "$name" "$secs" "${mem:-NA}" >> "$TIMING"
  echo "    ${name}: ${secs} s, peak memory ${mem:-NA} MB"
}

if $RUN_TESTS; then
  run_step tests Rscript tests/run_tests.R
  run_step toy_example Rscript example/toy_example.R
fi
run_step analysis Rscript analysis.R
run_step mutation_analysis Rscript mutation_analysis.R
if $RUN_REPORT; then
  run_step report Rscript -e 'rmarkdown::render("report.Rmd", params = list(project = "TCGA-UCEC"))'
fi

TOTAL_SECS=$(( SECONDS - TOTAL_START ))
printf "total\t%s\tNA\n" "$TOTAL_SECS" >> "$TIMING"

echo ""
echo "Timing (also in ${TIMING}):"
column -t -s $'\t' "$TIMING" 2>/dev/null || cat "$TIMING"
echo "Total: $(( TOTAL_SECS / 60 )) min $(( TOTAL_SECS % 60 )) s"
echo ""
echo "Done. Outputs are in results/; logs are in logs/."
