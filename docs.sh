#!/usr/bin/env bash
# Build the reference manual without configuring or compiling the application.
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
output="$project_dir/docs/generated"
open_docs=false
while (($#)); do
    case "$1" in
        --open) open_docs=true; shift ;;
        --output)
            if (($# < 2)); then echo '--output requires a directory' >&2; exit 2; fi
            output="$2"; shift 2 ;;
        --help|-h)
            echo 'Usage: docs.sh [--open] [--output DIRECTORY]'
            echo 'Requires Doxygen 1.9.8+ and Python 3. Set DOXYGEN to select its executable.'
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
doxygen_command="${DOXYGEN:-doxygen}"
if ! command -v "$doxygen_command" >/dev/null; then
    echo 'Doxygen is missing. On Ubuntu: sudo apt install doxygen' >&2
    exit 1
fi
mkdir -p -- "$output"
export OCTOTIGERII_DOC_OUTPUT="$(cd -- "$output" && pwd)"
cd -- "$project_dir"
python3 docs/check_references.py
"$doxygen_command" Doxyfile
printf 'Documentation: %s/html/index.html\n' "$OCTOTIGERII_DOC_OUTPUT"
if $open_docs; then
    if command -v xdg-open >/dev/null; then
        xdg-open "$OCTOTIGERII_DOC_OUTPUT/html/index.html"
    elif command -v open >/dev/null; then
        open "$OCTOTIGERII_DOC_OUTPUT/html/index.html"
    else
        echo 'Open that index.html in your browser.'
    fi
fi
