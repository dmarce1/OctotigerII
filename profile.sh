#!/usr/bin/env bash
set -euo pipefail

if (($# == 0)) || [[ "$1" == --help || "$1" == -h ]]; then
    cat <<'EOF'
Usage: ./profile.sh EXECUTABLE [ARGUMENTS...]

Run one locality with APEX screen, CSV, and TAU-format profiles. Every invocation
gets its own report directory under OCTOTIGERII_PROFILE_DIR (default: ./profiles).
Use this wrapper for each process in a distributed launch. The working
directory and application arguments are preserved.

APEX_PAPI_METRICS selects optional hardware events. Existing APEX_SCREEN_OUTPUT,
APEX_CSV_OUTPUT, and APEX_PROFILE_OUTPUT settings override this wrapper's defaults.
EOF
    exit 0
fi

profile_root="${OCTOTIGERII_PROFILE_DIR:-$PWD/profiles}"
mkdir -p -- "$profile_root"
profile_root="$(cd -- "$profile_root" && pwd)"
report_dir="$(mktemp -d "$profile_root/$(hostname).XXXXXXXX")"
export APEX_OUTPUT_FILE_PATH="$report_dir"
export APEX_SCREEN_OUTPUT="${APEX_SCREEN_OUTPUT:-1}"
export APEX_CSV_OUTPUT="${APEX_CSV_OUTPUT:-1}"
export APEX_PROFILE_OUTPUT="${APEX_PROFILE_OUTPUT:-1}"
printf 'APEX reports: %s\n' "$report_dir"
exec "$@"
