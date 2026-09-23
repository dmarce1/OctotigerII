#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./build.sh [release|debug|relwithdebinfo] --problem NAME --ndim N [-j JOBS]

Build HPX 1.11.0 in packages/TYPE/hpx and OctotigerII in TYPE/PROBLEM/NDIMd.
Defaults: --problem sod --ndim 1. See bin/*/*/CMakeLists.txt for supported dimensions.
Requires a C++20 compiler, CMake, and Git. Uses installed Boost and hwloc,
tries environment modules when available, then lets HPX fetch missing ones.
HPX fetches Asio and its matching APEX version. APEX and PAPI are enabled.
PAPI is reused from installed packages/modules or built locally if missing.
No system packages are installed or changed (no sudo).
Install Silo with HDF5; set CC, CXX, CMAKE_PREFIX_PATH, and Silo_ROOT as needed.
On clusters, load your preferred compiler module first. To select particular
dependency modules, set OCTOTIGERII_BOOST_MODULE, OCTOTIGERII_HWLOC_MODULE,
and OCTOTIGERII_PAPI_MODULE. Set Papi_ROOT (or PAPI_ROOT) for a custom PAPI.
OCTOTIGERII_PAPI=auto|system|build selects reuse/fallback, reuse-only, or a
local PAPI 7.2.0 build. Default: auto. See docs/profiling.md for running profiles.
EOF
}

problem=sod
ndim=1
build_type=Release
build_dir_name=release
jobs="${SLURM_CPUS_PER_TASK:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '2')}"
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
papi_dir="$project_dir/packages/$build_dir_name/papi"
papi_mode="${OCTOTIGERII_PAPI:-auto}"
case "$papi_mode" in
    auto|system|build) ;;
    *) printf 'OCTOTIGERII_PAPI must be auto, system, or build\n' >&2; exit 2 ;;
esac
papi_root="${Papi_ROOT:-${PAPI_ROOT:-${EBROOTPAPI:-}}}"

for program in cmake git; do
    command -v "$program" >/dev/null || { printf 'Missing required program: %s\n' "$program" >&2; exit 1; }
done

# CMake arguments use semicolons; environment modules may update this path.
cmake_prefix_path=
refresh_prefix_path() {
    cmake_prefix_path="${CMAKE_PREFIX_PATH:-}"
    cmake_prefix_path="${cmake_prefix_path//:/;}"
}

# Reject unsupported choices before fetching or building dependencies.
cmake "-DPROJECT_SOURCE_DIR=$project_dir" \
    "-DOCTOTIGERII_PROBLEM=$problem" "-DOCTOTIGERII_NDIM=$ndim" \
    -P "$project_dir/cmake/Problems.cmake"

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
if(CHECK_DEPENDENCY STREQUAL "BOOST")
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
elseif(CHECK_DEPENDENCY STREQUAL "HWLOC")
    find_path(HWLOC_INCLUDE_DIR hwloc.h
        HINTS "${Hwloc_ROOT}" PATH_SUFFIXES include)
    find_library(HWLOC_LIBRARY hwloc
        HINTS "${Hwloc_ROOT}" PATH_SUFFIXES lib lib64)
    if(NOT HWLOC_INCLUDE_DIR OR NOT HWLOC_LIBRARY)
        message(FATAL_ERROR "hwloc headers or library not found")
    endif()
elseif(CHECK_DEPENDENCY STREQUAL "PAPI")
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(PC_PAPI QUIET papi)
    endif()
    find_path(Papi_INCLUDE_DIR papi.h
        HINTS "${Papi_ROOT}" ${PC_PAPI_INCLUDE_DIRS} PATH_SUFFIXES include)
    find_library(Papi_LIBRARY NAMES papi
        HINTS "${Papi_ROOT}" ${PC_PAPI_LIBRARY_DIRS} PATH_SUFFIXES lib lib64)
    if(NOT Papi_INCLUDE_DIR OR NOT Papi_LIBRARY)
        message(FATAL_ERROR "PAPI headers or library not found")
    endif()
    # Compile/link only: counter permissions need not be granted on login nodes.
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_INCLUDES "${Papi_INCLUDE_DIR}")
    set(CMAKE_REQUIRED_LIBRARIES "${Papi_LIBRARY}")
    check_cxx_source_compiles("#include <papi.h>
        int main() { return PAPI_library_init(PAPI_VER_CURRENT) < 0; }" PAPI_LINKS)
    if(NOT PAPI_LINKS)
        message(FATAL_ERROR "PAPI cannot be linked with the selected compiler")
    endif()
    file(WRITE "${CMAKE_BINARY_DIR}/papi-include" "${Papi_INCLUDE_DIR}\n")
    file(WRITE "${CMAKE_BINARY_DIR}/papi-library" "${Papi_LIBRARY}\n")
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
    papi_root="${papi_root:-${Papi_ROOT:-${PAPI_ROOT:-${EBROOTPAPI:-}}}}"
    dependency_args=()
    [[ -z "$boost_root" ]] || dependency_args+=("-DBoost_ROOT=$boost_root")
    [[ -z "$hwloc_root" ]] || dependency_args+=("-DHwloc_ROOT=$hwloc_root")
    [[ -z "$papi_root" ]] || dependency_args+=("-DPapi_ROOT=$papi_root" "-DPAPI_ROOT=$papi_root")
}

probe_dependency() {
    local requested="$1"
    rm -rf -- "$probe_dir/build"
    refresh_prefix_path
    refresh_dependency_roots
    cmake -S "$probe_dir" -B "$probe_dir/build" \
        "-DCHECK_DEPENDENCY=$requested" \
        "-DCMAKE_PREFIX_PATH=$cmake_prefix_path" \
        "${compiler_args[@]}" "${dependency_args[@]}" \
        > "$probe_dir/log" 2>&1
}

fetch_boost=OFF
if ! probe_dependency BOOST; then
    if type module >/dev/null 2>&1; then
        module load "${OCTOTIGERII_BOOST_MODULE:-boost}" || true
    fi
    if ! probe_dependency BOOST; then
        fetch_boost=ON
        printf 'Boost not found; HPX will fetch it.\n'
    fi
fi

fetch_hwloc=OFF
if ! probe_dependency HWLOC; then
    if type module >/dev/null 2>&1; then
        module load "${OCTOTIGERII_HWLOC_MODULE:-hwloc}" || true
    fi
    if ! probe_dependency HWLOC; then
        fetch_hwloc=ON
        printf 'hwloc not found; HPX will fetch it.\n'
    fi
fi
# PAPI is a separate development library, not a standard Linux installation.
# Let site modules/system packages take precedence; build without root if absent.
papi_found=OFF
if [[ "$papi_mode" != build ]]; then
    if probe_dependency PAPI; then
        papi_found=ON
    else
        if type module >/dev/null 2>&1; then
            module load "${OCTOTIGERII_PAPI_MODULE:-papi}" || true
        fi
        if probe_dependency PAPI; then papi_found=ON; fi
    fi
fi
if [[ "$papi_found" == OFF && "$papi_mode" == system ]]; then
    cat "$probe_dir/log" >&2
    printf 'PAPI development files missing; load a module or set Papi_ROOT.\n' >&2
    exit 1
fi
if [[ "$papi_found" == OFF ]]; then
    papi_root="$papi_dir/install"
    if [[ "$papi_mode" != build ]] && probe_dependency PAPI; then
        papi_found=ON
    else
        for program in make cc; do
            # CC can name a non-default compiler; the configure step uses it.
            [[ "$program" != cc || -z "${CC:-}" ]] || continue
            command -v "$program" >/dev/null || { printf 'Missing required program: %s\n' "$program" >&2; exit 1; }
        done
        papi_src="$papi_dir/src"
        papi_tag=papi-7-2-0-t
        if [[ ! -d "$papi_src/.git" ]]; then
            [[ ! -e "$papi_src" ]] || { printf 'Existing PAPI source is not a Git checkout: %s\n' "$papi_src" >&2; exit 1; }
            mkdir -p "$papi_dir"
            git clone --depth 1 --branch "$papi_tag" https://github.com/icl-utk-edu/papi.git "$papi_src"
        fi
        [[ "$(git -C "$papi_src" describe --tags --exact-match HEAD 2>/dev/null)" == "$papi_tag" ]] || {
            printf 'Expected PAPI %s in %s\n' "$papi_tag" "$papi_src" >&2; exit 1;
        }
        printf 'Building PAPI 7.2.0 in %s\n' "$papi_dir"
        (
            cd -- "$papi_src/src"
            ./configure --prefix="$papi_root" --with-shared-lib=yes --with-static-lib=yes --with-tests=no
            make -j "$jobs"
            make install
        )
        if ! probe_dependency PAPI; then
            cat "$probe_dir/log" >&2
            printf 'Local PAPI installation failed the compile/link check.\n' >&2
            exit 1
        fi
    fi
fi
# Pin both find modules to the same checked library, including lib64/multiarch.
IFS= read -r papi_include < "$probe_dir/build/papi-include"
IFS= read -r papi_library < "$probe_dir/build/papi-library"
papi_args=("-DPapi_INCLUDE_DIR=$papi_include" "-DPapi_LIBRARY=$papi_library"
    "-DPAPI_INCLUDE_DIR=$papi_include" "-DPAPI_LIBRARY=$papi_library")
printf 'Using PAPI: %s\n' "$papi_library"
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
    -DHPX_WITH_APEX=ON \
    -DHPX_WITH_FETCH_APEX=ON \
    -DHPX_WITH_PAPI=ON \
    -DAPEX_WITH_PAPI=ON \
    -DHPX_WITH_THREAD_DESCRIPTION=ON \
    -DHPX_WITH_THREAD_IDLE_RATES=ON \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
    -DHPX_WITH_TESTS=OFF \
    -DHPX_WITH_EXAMPLES=OFF \
    "-DHPX_WITH_FETCH_BOOST=$fetch_boost" \
    -DHPX_WITH_FETCH_ASIO=ON \
    "-DHPX_WITH_FETCH_HWLOC=$fetch_hwloc" \
    "-DCMAKE_PREFIX_PATH=$cmake_prefix_path" \
    "${compiler_args[@]}" "${dependency_args[@]}" "${papi_args[@]}"
cmake "-DHPX_BINARY_DIR=$hpx_build" -P "$project_dir/cmake/PatchApex.cmake"
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
    "-DOCTOTIGERII_PROBLEM=$problem" \
    "-DOCTOTIGERII_NDIM=$ndim" \
    "${compiler_args[@]}" "${dependency_args[@]}" "${papi_args[@]}"
cmake --build "$octo_build" --parallel "$jobs"
printf 'Executable: %s/octotigerII-%s-%sd\n' "$octo_build" "$problem" "$ndim"
