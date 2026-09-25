#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./build.sh [release|debug|relwithdebinfo] [-j JOBS] [-DOCTOII_WITH_HYDRO=OFF ...]

Build octoII-1d, octoII-2d, octoII-3d and the octoII link in TYPE/.
HYDRO, RADIATION and GRAVITY are all ON by default.
Requires a C++20 compiler, CMake, and Git. Uses installed Boost and hwloc,
tries environment modules when available, then lets HPX fetch missing ones.
HPX also fetches Asio and APEX. HPX APEX support and Octo-II profiling
are enabled by default. No system packages are installed or changed.
Install Silo with HDF5; set CC, CXX, CMAKE_PREFIX_PATH, and Silo_ROOT as needed.
On clusters, load your preferred compiler module first. To select particular
dependency modules, set OCTOTIGERII_BOOST_MODULE and OCTOTIGERII_HWLOC_MODULE.
EOF
}

physics_args=()
build_type=Release
build_dir_name=release
jobs="${SLURM_CPUS_PER_TASK:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '2')}"
while (($#)); do
    case "$1" in
        release|Release) build_type=Release; build_dir_name=release ;;
        debug|Debug) build_type=Debug; build_dir_name=debug ;;
        relwithdebinfo|RelWithDebInfo) build_type=RelWithDebInfo; build_dir_name=relwithdebinfo ;;
        -DOCTOII_WITH_HYDRO=*|-DOCTOII_WITH_RADIATION=*|-DOCTOII_WITH_GRAVITY=*) physics_args+=("$1") ;;
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

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
hpx_dir="$project_dir/packages/$build_dir_name/hpx"
hpx_src="$hpx_dir/src"
hpx_build="$hpx_dir/build"
hpx_install="$hpx_dir/install"
octo_build="$project_dir/$build_dir_name"

for program in cmake git; do
    command -v "$program" >/dev/null || { printf 'Missing required program: %s\n' "$program" >&2; exit 1; }
done

# CMake arguments use semicolons; environment modules may update this path.
cmake_prefix_path=
refresh_prefix_path() {
    cmake_prefix_path="${CMAKE_PREFIX_PATH:-}"
    cmake_prefix_path="${cmake_prefix_path//:/;}"
}

# Use one compiler pair for both CMake projects, including repeated builds.
compiler_args=()
[[ -z "${CC:-}" ]] || compiler_args+=("-DCMAKE_C_COMPILER=$CC")
[[ -z "${CXX:-}" ]] || compiler_args+=("-DCMAKE_CXX_COMPILER=$CXX")

# A tiny separate CMake project checks the *same* discovery paths that the
# real builds will use. Module systems need not be installed on a workstation.
probe_dir="$(mktemp -d)"
trap 'rm -rf -- "$probe_dir"' EXIT
cat > "$probe_dir/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.18)
project(octotiger_dependency_probe C CXX)
if(CHECK_BOOST)
    find_package(Boost 1.71 QUIET)
    if(NOT Boost_FOUND)
        message(FATAL_ERROR "Boost 1.71+ not found")
    endif()
    set(found_units FALSE)
    foreach(dir IN LISTS Boost_INCLUDE_DIRS)
        if(EXISTS "${dir}/boost/units/quantity.hpp")
            set(found_units TRUE)
        endif()
    endforeach()
    if(NOT found_units)
        message(FATAL_ERROR "Boost.Units headers not found")
    endif()
else()
    find_path(HWLOC_INCLUDE_DIR hwloc.h
        HINTS "${Hwloc_ROOT}" PATH_SUFFIXES include)
    find_library(HWLOC_LIBRARY hwloc
        HINTS "${Hwloc_ROOT}" PATH_SUFFIXES lib lib64)
    if(NOT HWLOC_INCLUDE_DIR OR NOT HWLOC_LIBRARY)
        message(FATAL_ERROR "hwloc headers or library not found")
    endif()
endif()
EOF

# LONI uses Environment Modules. Loading a module changes this script's shell
# environment only; it requires no root privileges or system package manager.
if ! type module >/dev/null 2>&1 && [[ -r /etc/profile.d/modules.sh ]]; then
    source /etc/profile.d/modules.sh
fi

dependency_args=()
refresh_dependency_roots() {
    local boost_root="${Boost_ROOT:-${BOOST_ROOT:-${EBROOTBOOST:-}}}"
    local hwloc_root="${Hwloc_ROOT:-${HWLOC_ROOT:-${EBROOTHWLOC:-}}}"
    dependency_args=()
    [[ -z "$boost_root" ]] || dependency_args+=("-DBoost_ROOT=$boost_root")
    [[ -z "$hwloc_root" ]] || dependency_args+=("-DHwloc_ROOT=$hwloc_root")
}

probe_dependency() {
    local requested="$1"
    rm -rf -- "$probe_dir/build"
    refresh_prefix_path
    refresh_dependency_roots
    cmake -S "$probe_dir" -B "$probe_dir/build" \
        "-DCHECK_BOOST=$requested" \
        "-DCMAKE_PREFIX_PATH=$cmake_prefix_path" \
        "${compiler_args[@]}" "${dependency_args[@]}" \
        > "$probe_dir/log" 2>&1
}

fetch_boost=OFF
if ! probe_dependency ON; then
    if type module >/dev/null 2>&1; then
        module load "${OCTOTIGERII_BOOST_MODULE:-boost}" || true
    fi
    if ! probe_dependency ON; then
        fetch_boost=ON
        printf 'Boost not found; HPX will fetch it.\n'
    fi
fi

fetch_hwloc=OFF
if ! probe_dependency OFF; then
    if type module >/dev/null 2>&1; then
        module load "${OCTOTIGERII_HWLOC_MODULE:-hwloc}" || true
    fi
    if ! probe_dependency OFF; then
        fetch_hwloc=ON
        printf 'hwloc not found; HPX will fetch it.\n'
    fi
fi
refresh_prefix_path
refresh_dependency_roots

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
    -DHPX_WITH_APEX=ON \
    -DHPX_WITH_FETCH_APEX=ON \
    "-DHPX_WITH_FETCH_BOOST=$fetch_boost" \
    -DHPX_WITH_FETCH_ASIO=ON \
    "-DHPX_WITH_FETCH_HWLOC=$fetch_hwloc" \
    "-DCMAKE_PREFIX_PATH=$cmake_prefix_path" \
    "${compiler_args[@]}" "${dependency_args[@]}"
# The APEX revision fetched by HPX 1.11.0 uses uint64_t in gzstream.hpp
# without including <cstdint>. Apply the missing include after FetchContent
# has populated the source, including on fresh builds.
apex_gzstream="$hpx_build/_deps/apex-src/src/apex/gzstream.hpp"
if [[ -f "$apex_gzstream" ]] &&
    ! grep -Eq '^[[:space:]]*#[[:space:]]*include[[:space:]]*[<"](cstdint|stdint\.h)[>"]' "$apex_gzstream"; then
    printf 'Adding missing <cstdint> include to APEX gzstream.hpp\n'
    sed -i '1i#include <cstdint>' "$apex_gzstream"
fi

cmake --build "$hpx_build" --parallel "$jobs"
cmake --install "$hpx_build"
hpx_config="$(find "$hpx_install" -name HPXConfig.cmake -print -quit)"
[[ -n "$hpx_config" ]] || { printf 'HPXConfig.cmake missing from %s\n' "$hpx_install" >&2; exit 1; }

# OctotigerII also uses Boost.Units. When HPX had to fetch Boost, expose its
# installation to the second CMake project as well.
if [[ "$fetch_boost" == ON ]]; then
    fetched_boost="$hpx_install"
    units_header="$(find "$hpx_install" "$hpx_build/_deps" -type f \
        -path '*/boost/units/quantity.hpp' -print -quit 2>/dev/null || true)"
    if [[ -z "$units_header" ]]; then
        printf 'HPX fetched Boost, but Boost.Units headers are missing.\n' >&2
        printf 'Set Boost_ROOT to a complete Boost 1.71+ installation and retry.\n' >&2
        exit 1
    fi
    dependency_args+=("-DBoost_INCLUDEDIR=${units_header%/boost/units/quantity.hpp}")
    dependency_args+=("-DBoost_ROOT=$fetched_boost")
    cmake_prefix_path="$fetched_boost${cmake_prefix_path:+;$cmake_prefix_path}"
fi

# Resolve the HPX just installed, rather than a different HPX on the host.
printf 'Building OctotigerII %s in %s\n' "$build_type" "$octo_build"
cmake -S "$project_dir" -B "$octo_build" \
    "-DCMAKE_BUILD_TYPE=$build_type" \
    "-DHPX_DIR=$(dirname "$hpx_config")" \
    "-DCMAKE_PREFIX_PATH=$hpx_install${cmake_prefix_path:+;$cmake_prefix_path}" \
    -DOCTOTIGERII_WITH_HPX=ON \
    -DOCTOTIGERII_WITH_PROFILING=ON \
    -DOCTOTIGERII_BUILD_TESTS=ON \
    "${physics_args[@]}" \
    "${compiler_args[@]}" "${dependency_args[@]}"
cmake --build "$octo_build" --parallel "$jobs"
printf 'Executables: %s/octoII-{1d,2d,3d}; octoII -> octoII-3d\n' "$octo_build"
