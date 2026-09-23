#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./build.sh [release|debug|relwithdebinfo] --problem NAME --ndim N [-j JOBS]

Build HPX 1.11.0 in packages/TYPE/hpx and OctotigerII in TYPE/PROBLEM/NDIMd.
Defaults: --problem sod --ndim 1. See bin/*/*/CMakeLists.txt for supported dimensions.
The compiler, CMake, Git, Boost, and hwloc must be installed already.
HPX fetches its Asio headers during configuration.
Install Silo with HDF5; set CC, CXX, CMAKE_PREFIX_PATH, and Silo_ROOT as needed.
EOF
}

problem=sod
ndim=1
build_type=Release
build_dir_name=release
jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '2')"
while (($#)); do
    case "$1" in
        release|Release) build_type=Release; build_dir_name=release ;;
        debug|Debug) build_type=Debug; build_dir_name=debug ;;
        relwithdebinfo|RelWithDebInfo) build_type=RelWithDebInfo; build_dir_name=relwithdebinfo ;;
        --problem|--ndim)
            (($# >= 2)) || { usage >&2; exit 2; }
            if [[ "$1" == --problem ]]; then problem="$2"; else ndim="$2"; fi
            shift
            ;;
        --problem=*) problem="${1#*=}" ;;
        --ndim=*) ndim="${1#*=}" ;;
        -j|--jobs)
            (($# >= 2)) || { usage >&2; exit 2; }
            jobs="$2"
            shift
            ;;
        -h|--help) usage; exit 0 ;;
        *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { printf 'Invalid job count: %s\n' "$jobs" >&2; exit 2; }

[[ "$problem" =~ ^[a-z][a-z0-9_-]*$ ]] || { echo 'Invalid problem name' >&2; exit 2; }
[[ "$ndim" =~ ^[123]$ ]] || { echo 'ndim must be 1, 2, or 3' >&2; exit 2; }

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
hpx_dir="$project_dir/packages/$build_dir_name/hpx"
hpx_src="$hpx_dir/src"
hpx_build="$hpx_dir/build"
hpx_install="$hpx_dir/install"
octo_build="$project_dir/$build_dir_name/$problem/${ndim}d"

for program in cmake git; do
    command -v "$program" >/dev/null || { printf 'Missing required program: %s\n' "$program" >&2; exit 1; }
done

# Reject unsupported choices before fetching or building dependencies.
cmake "-DPROJECT_SOURCE_DIR=$project_dir" \
    "-DOCTOTIGERII_PROBLEM=$problem" "-DOCTOTIGERII_NDIM=$ndim" \
    -P "$project_dir/cmake/Problems.cmake"

# Use one compiler pair for both CMake projects, including repeated builds.
compiler_args=()
[[ -z "${CC:-}" ]] || compiler_args+=("-DCMAKE_C_COMPILER=$CC")
[[ -z "${CXX:-}" ]] || compiler_args+=("-DCMAKE_CXX_COMPILER=$CXX")

if [[ ! -d "$hpx_src/.git" ]]; then
    [[ ! -e "$hpx_src" ]] || { printf 'Existing HPX source is not a Git checkout: %s\n' "$hpx_src" >&2; exit 1; }
    mkdir -p "$hpx_dir"
    git clone --branch v1.11.0 --depth 1 https://github.com/STEllAR-GROUP/hpx.git "$hpx_src"
fi
[[ "$(git -C "$hpx_src" describe --tags --exact-match HEAD 2>/dev/null)" == v1.11.0 ]] || {
    printf 'Expected HPX v1.11.0 in %s\n' "$hpx_src" >&2
    exit 1
}

printf 'Building HPX %s in %s\n' "$build_type" "$hpx_dir"
cmake -S "$hpx_src" -B "$hpx_build" \
    "-DCMAKE_BUILD_TYPE=$build_type" \
    "-DCMAKE_INSTALL_PREFIX=$hpx_install" \
    -DHPX_WITH_CXX_STANDARD=20 \
    -DHPX_WITH_DISTRIBUTED_RUNTIME=ON \
    -DHPX_WITH_NETWORKING=ON \
    -DHPX_WITH_PARCELPORT_TCP=ON \
    -DHPX_WITH_MALLOC=system \
    -DHPX_WITH_TESTS=OFF \
    -DHPX_WITH_EXAMPLES=OFF \
    -DHPX_WITH_FETCH_BOOST=OFF \
    -DHPX_WITH_FETCH_ASIO=ON \
    -DHPX_WITH_FETCH_HWLOC=OFF \
    "${compiler_args[@]}"
cmake --build "$hpx_build" --parallel "$jobs"
cmake --install "$hpx_build"
hpx_config="$(find "$hpx_install" -name HPXConfig.cmake -print -quit)"
[[ -n "$hpx_config" ]] || { printf 'HPXConfig.cmake missing from %s\n' "$hpx_install" >&2; exit 1; }

# Resolve the HPX just installed, rather than a different HPX on the host.
printf 'Building OctotigerII %s in %s\n' "$build_type" "$octo_build"
cmake -S "$project_dir" -B "$octo_build" \
    "-DCMAKE_BUILD_TYPE=$build_type" \
    "-DHPX_DIR=$(dirname "$hpx_config")" \
    "-DCMAKE_PREFIX_PATH=$hpx_install${CMAKE_PREFIX_PATH:+;$CMAKE_PREFIX_PATH}" \
    -DOCTOTIGERII_WITH_HPX=ON \
    "-DOCTOTIGERII_PROBLEM=$problem" \
    "-DOCTOTIGERII_NDIM=$ndim" \
    "${compiler_args[@]}"
cmake --build "$octo_build" --parallel "$jobs"
printf 'Executable: %s/octotigerII-%s-%sd\n' "$octo_build" "$problem" "$ndim"
