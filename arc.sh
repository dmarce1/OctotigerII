#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_name="$(basename -- "$project_dir")"
project_parent="$(dirname -- "$project_dir")"
archive="${1:-$project_parent/${project_name}-source.tar.gz}"

if (($# > 1)); then
    printf 'Usage: %s [archive.tar.gz]\n' "${0##*/}" >&2
    exit 2
fi
if [[ ! -d "$project_dir/.git" ]]; then
    printf 'No .git directory in %s; cannot include Git history\n' "$project_dir" >&2
    exit 1
fi
if [[ "$archive" != /* ]]; then
    archive="$PWD/$archive"
fi
archive_parent="$(dirname -- "$archive")"
mkdir -p -- "$archive_parent"
temporary="$(mktemp "$archive_parent/.${project_name}-source.XXXXXXXX.tar.gz")"
trap 'rm -f -- "$temporary"' EXIT

tar -C "$project_parent" -czf "$temporary" \
    --exclude="$project_name/packages" \
    --exclude="$project_name/release" \
    --exclude="$project_name/debug" \
    --exclude="$project_name/relwithdebinfo" \
    --exclude="$project_name/build-*" \
    --exclude="$project_name/output*" \
    "$project_name"
mv -f -- "$temporary" "$archive"
printf 'Created %s\n' "$archive"
