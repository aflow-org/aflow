#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 <build> [-i] [-j <jobs>]"
  echo " <build>  where to find compile_commands.json"
  echo "      -i  apply changes in-place (default: dry-run)"
  echo "      -j  number of parallel jobs to run (defualt: 0) *note* need GNU parallel to use this"
  exit 1
}

if [[ $# -lt 1 ]]; then
  usage
fi
bld=$1
shift 1;

DO_INPLACE=false
njobs=0

while [[ $# -gt 0 ]]; do
  opt="$1"
  case "$opt" in
    "-i" )
      DO_INPLACE=true
      shift 1;;
    "-j" )
      njobs=$2
      shift 2;;
    "-h" | "--help" )
      usage;
      exit 0;;
    *)
      echo "invalid option $opt"; exit 1;;
  esac
done

if $DO_INPLACE; then
  arg="--fix"
else
  arg=""
  echo "-i not given, performing dry run"
fi

# macOS-specific arguments for clang-tidy
DARWIN_ARGS=""
if [[ "$(uname)" == "Darwin" ]]; then
  SDK=$(xcrun --sdk macosx --show-sdk-path)
  RES=$(clang++ -print-resource-dir)
  DARWIN_ARGS="-extra-arg=-isysroot -extra-arg=$SDK -extra-arg=-resource-dir=$RES"
fi

# Produce one clang-tidy command line per affected file (but don't run them yet).
# Each printed line looks like:
# clang-tidy --use-color=true --config-file .clang-tidy-lint -p build/gcc-debug --format-style=file --line-filter=[{"name":"src/aflow.cpp","lines":[[226,270]]}] "src/aflow.cpp"
tmp_cmds=$(mktemp)
trap 'rm -f "$tmp_cmds"' EXIT
git diff -U3 --no-color staging -- '*.cpp' '*.h' '*.tpp' \
  | awk -v arg="$arg" -v bld="$bld" -v darwin_args="$DARWIN_ARGS" '
    /^\+\+\+ b\// { file=substr($0,7); next }
    /^@@/ && file != "" {
      if (match($0, /\+[0-9]+/)) {
        start_str = substr($0, RSTART+1, RLENGTH-1)
        start = start_str + 0
        start_pos = RSTART + RLENGTH
        rest = substr($0, start_pos)
        if (match(rest, /,[0-9]+/)) {
          count_str = substr(rest, RSTART+1, RLENGTH-1)
          count = count_str + 0
        } else {
          count = 1
        }
        if (count == 0) next
        end = start + count
        ranges[file] = ranges[file]  ",[" start "," end "]"
        files[file] = 1
      }
    }
    END {
      for (f in files) {
        printf("clang-tidy " arg " --use-color=true --config-file .clang-tidy-lint -p " bld " --format-style=file " darwin_args " --line-filter='\''[{\"name\":\"%s\",\"lines\":[%s]}]'\'' \"%s\"\n", f, ranges[f], f)
      }
    }' \
  | sed -e 's/:\[,/:\[/g' > "$tmp_cmds"

CMDS=()
while IFS= read -r line; do
  CMDS+=("$line")
done < "$tmp_cmds"
rm -f "$tmp_cmds"

# No changed lines -> nothing to do.
if [[ ${#CMDS[@]} -eq 0 ]]; then
  echo "No changed lines to tidy."
  exit 0
fi

ANY_ISSUES=0

if [[ "$njobs" -gt 0 ]]; then
  command -v parallel >/dev/null 2>&1 || { echo "need GNU parallel to run in parallel"; exit 1; }
  if printf '%s\n' "${CMDS[@]}" | parallel --jobs "$njobs" --group --will-cite 'bash -lc {}'
  then
    :
  else
    ANY_ISSUES=1;
  fi
else
  for cmd in "${CMDS[@]}"; do
    if eval "$cmd"; then
      # success for this file (no immediate output)
      :
    else
      # unlikely: clang-format -i usually returns 0; capture if it didn't
      ANY_ISSUES=1;
    fi
  done
fi

if $DO_INPLACE; then
  if [[ $ANY_ISSUES -eq 0 ]]; then
    echo "Fixes applied."
    exit 0
  else
    echo "Some fixes applied. Some issues detected."
    echo "Should check input/output and run again."
    exit 2
  fi
else
  # Dry-run mode: report results
  if [[ $ANY_ISSUES -eq 0 ]]; then
    echo "Looks good! No tidy fixes required for changed lines."
    exit 0
  else
    echo "Tidy issues detected for changed lines"
    echo "Run with -i to apply fixes: $0 $@ -i"
    exit 2
  fi
fi

