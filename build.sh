#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./build.sh [release|debug|relwithdebinfo] [--hpx-debug|--no-hpx-debug] [-j JOBS] [-DOCTOII_WITH_HYDRO=OFF ...]

Build octoII-1d, octoII-2d, octoII-3d and the octoII link in TYPE/.
HPX diagnostics default OFF for release, ON for debug and relwithdebinfo.
--hpx-debug enables HPX lock verification, lock backtraces, thread
debug information, spinlock detection, and scheduling counters. Diagnostic
builds use TYPE-hpxdebug/ and packages/TYPE-hpxdebug/, keeping normal builds.
--no-hpx-debug disables those options and uses TYPE/.
HYDRO, RADIATION and GRAVITY are all ON by default.
All C and C++ targets are optimized for the build host with -march=native.
Requires a C++20 compiler, CMake, and Git. Uses installed Boost and hwloc,
tries environment modules when available, then lets HPX fetch missing ones.
HPX also fetches Asio and APEX. HPX APEX support and Octo-II profiling
are enabled by default. No system packages are installed or changed.
Uses installed HDF5 and Silo when available, or builds them locally.
Set OCTOTIGERII_BUILD_HDF5=ON to force a local HDF5 build, for example when
the cluster's HDF5 module uses a different compiler. Set CC, CXX,
CMAKE_PREFIX_PATH, HDF5_ROOT, and Silo_ROOT to select installations.
On clusters, load your preferred compiler module first. To select particular
dependency modules, set OCTOTIGERII_BOOST_MODULE and OCTOTIGERII_HWLOC_MODULE.
EOF
}

physics_args=()
build_type=Release
build_dir_name=release
hpx_debug=
jobs="${SLURM_CPUS_PER_TASK:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '2')}"
while (($#)); do
    case "$1" in
        release|Release) build_type=Release; build_dir_name=release ;;
        debug|Debug) build_type=Debug; build_dir_name=debug ;;
        relwithdebinfo|RelWithDebInfo) build_type=RelWithDebInfo; build_dir_name=relwithdebinfo ;;
        --hpx-debug) hpx_debug=ON ;;
        --no-hpx-debug) hpx_debug=OFF ;;
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
if [[ -z "$hpx_debug" ]]; then
    if [[ "$build_type" == Release ]]; then hpx_debug=OFF; else hpx_debug=ON; fi
fi
if [[ "$hpx_debug" == ON ]]; then build_dir_name+="-hpxdebug"; fi

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
local_prefix_path=
refresh_prefix_path() {
    cmake_prefix_path="${CMAKE_PREFIX_PATH:-}"
    cmake_prefix_path="${cmake_prefix_path//:/;}"
    cmake_prefix_path="${local_prefix_path:+$local_prefix_path;}$cmake_prefix_path"
}

# Use one compiler pair for both CMake projects, including repeated builds.
compiler_args=()
[[ -z "${CC:-}" ]] || compiler_args+=("-DCMAKE_C_COMPILER=$CC")
[[ -z "${CXX:-}" ]] || compiler_args+=("-DCMAKE_CXX_COMPILER=$CXX")

# OctotigerII is an HPC application built for the machine where it will run.
# Apply the native ISA to HPX/APEX and the application so the compiler can use
# the host's vector width and fused operations consistently.
native_arch_args=(
    "-DCMAKE_C_FLAGS=${CFLAGS:+$CFLAGS }-march=native"
    "-DCMAKE_CXX_FLAGS=${CXXFLAGS:+$CXXFLAGS }-march=native"
)

# A tiny separate CMake project checks the *same* discovery paths that the
# real builds will use. Module systems need not be installed on a workstation.
probe_dir="$(mktemp -d)"
trap 'rm -rf -- "$probe_dir"' EXIT
cat > "$probe_dir/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.18)
project(octotiger_dependency_probe C CXX)
if(CHECK_HDF5)
    find_package(HDF5 QUIET COMPONENTS C)
    if(NOT HDF5_FOUND)
        message(FATAL_ERROR "HDF5 C headers or library not found")
    endif()
elseif(CHECK_SILO)
    find_path(Silo_INCLUDE_DIR NAMES silo.h HINTS "${Silo_ROOT}" PATH_SUFFIXES include)
    find_library(Silo_LIBRARY NAMES siloh5 HINTS "${Silo_ROOT}" PATH_SUFFIXES lib lib64)
    if(NOT Silo_INCLUDE_DIR OR NOT Silo_LIBRARY)
        message(FATAL_ERROR "Silo headers or HDF5-enabled library not found")
    endif()
elseif(CHECK_BOOST)
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
silo_args=()
hdf5_args=()
silo_hdf5_args=()
local_hdf5_install=
refresh_dependency_roots() {
    local boost_root="${Boost_ROOT:-${BOOST_ROOT:-${EBROOTBOOST:-}}}"
    local hwloc_root="${Hwloc_ROOT:-${HWLOC_ROOT:-${EBROOTHWLOC:-}}}"
    local silo_root="${Silo_ROOT:-${SILO_ROOT:-${EBROOTSILO:-}}}"
    local hdf5_root="${local_hdf5_install:-${HDF5_ROOT:-${EBROOTHDF5:-}}}"
    dependency_args=()
    silo_args=()
    hdf5_args=()
    silo_hdf5_args=()
    [[ -z "$boost_root" ]] || dependency_args+=("-DBoost_ROOT=$boost_root")
    [[ -z "$hwloc_root" ]] || dependency_args+=("-DHwloc_ROOT=$hwloc_root")
    [[ -z "$silo_root" ]] || silo_args+=("-DSilo_ROOT=$silo_root")
    [[ -z "$hdf5_root" ]] || hdf5_args+=("-DHDF5_ROOT=$hdf5_root")
    [[ -z "$hdf5_root" ]] || silo_hdf5_args+=("-DSILO_HDF5_DIR=$hdf5_root")
}

probe_dependency() {
    local requested="$1"
    local check_boost="$requested"
    local check_silo=OFF
    local check_hdf5=OFF
    if [[ "$requested" == SILO || "$requested" == HDF5 ]]; then
        check_boost=OFF
        [[ "$requested" != SILO ]] || check_silo=ON
        [[ "$requested" != HDF5 ]] || check_hdf5=ON
    fi
    rm -rf -- "$probe_dir/build"
    refresh_prefix_path
    refresh_dependency_roots
    cmake -S "$probe_dir" -B "$probe_dir/build" \
        "-DCHECK_BOOST=$check_boost" \
        "-DCHECK_SILO=$check_silo" \
        "-DCHECK_HDF5=$check_hdf5" \
        "-DCMAKE_PREFIX_PATH=$cmake_prefix_path" \
        "${compiler_args[@]}" "${dependency_args[@]}" "${silo_args[@]}" "${hdf5_args[@]}" \
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

if [[ "${OCTOTIGERII_BUILD_HDF5:-OFF}" == ON ]] || ! probe_dependency HDF5; then
    hdf5_dir="$project_dir/packages/$build_dir_name/hdf5"
    hdf5_src="$hdf5_dir/src"
    hdf5_build="$hdf5_dir/build"
    hdf5_install="$hdf5_dir/install"
    if [[ ! -d "$hdf5_src/.git" ]]; then
        [[ ! -e "$hdf5_src" ]] || { printf 'Existing HDF5 source is not a Git checkout: %s\n' "$hdf5_src" >&2; exit 1; }
        mkdir -p "$hdf5_dir"
        git clone --branch hdf5_1.14.6 --depth 1 https://github.com/HDFGroup/hdf5.git "$hdf5_src"
    fi
    git -C "$hdf5_src" tag --points-at HEAD | grep -Fxq hdf5_1.14.6 || {
        printf 'Expected HDF5 1.14.6 in %s\n' "$hdf5_src" >&2
        exit 1
    }
    printf 'Building HDF5 %s in %s\n' "$build_type" "$hdf5_dir"
    cmake -S "$hdf5_src" -B "$hdf5_build" \
        "-DCMAKE_BUILD_TYPE=$build_type" \
        "-DCMAKE_INSTALL_PREFIX=$hdf5_install" \
        -DHDF5_INSTALL_LIB_DIR=lib \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_STATIC_LIBS=OFF \
        -DHDF5_ENABLE_PARALLEL=OFF \
        -DHDF5_ENABLE_SZIP_SUPPORT=OFF \
        -DHDF5_ENABLE_Z_LIB_SUPPORT=ON \
        -DHDF5_BUILD_HL_LIB=OFF \
        -DHDF5_BUILD_CPP_LIB=OFF \
        -DHDF5_BUILD_FORTRAN=OFF \
        -DHDF5_BUILD_TOOLS=OFF \
        -DHDF5_BUILD_EXAMPLES=OFF \
        -DBUILD_TESTING=OFF \
        "${native_arch_args[@]}" \
        "${compiler_args[@]}"
    cmake --build "$hdf5_build" --parallel "$jobs"
    cmake --install "$hdf5_build"
    [[ -f "$hdf5_install/include/hdf5.h" && -f "$hdf5_install/lib/libhdf5.so" ]] || {
        printf 'HDF5 header or library missing from %s\n' "$hdf5_install" >&2
        exit 1
    }
    local_hdf5_install="$hdf5_install"
    local_prefix_path="$hdf5_install${local_prefix_path:+;$local_prefix_path}"
    refresh_prefix_path
    refresh_dependency_roots
fi

if [[ -z "$local_hdf5_install" ]] && ! probe_dependency SILO; then
    if type module >/dev/null 2>&1; then
        module load "${OCTOTIGERII_SILO_MODULE:-silo}" || true
    fi
fi
if [[ -n "$local_hdf5_install" ]] || ! probe_dependency SILO; then
        silo_dir="$project_dir/packages/$build_dir_name/silo"
        silo_src="$silo_dir/src"
        silo_build="$silo_dir/build"
        silo_install="$silo_dir/install"
        if [[ ! -d "$silo_src/.git" ]]; then
            [[ ! -e "$silo_src" ]] || { printf 'Existing Silo source is not a Git checkout: %s\n' "$silo_src" >&2; exit 1; }
            mkdir -p "$silo_dir"
            git clone --branch 4.12.1 --depth 1 https://github.com/LLNL/Silo.git "$silo_src"
        fi
        [[ "$(git -C "$silo_src" describe --tags --exact-match HEAD 2>/dev/null)" == 4.12.1 ]] || {
            printf 'Expected Silo 4.12.1 in %s\n' "$silo_src" >&2
            exit 1
        }
        printf 'Building Silo %s in %s\n' "$build_type" "$silo_dir"
        cmake -S "$silo_src" -B "$silo_build" \
            "-DCMAKE_BUILD_TYPE=$build_type" \
            "-DCMAKE_INSTALL_PREFIX=$silo_install" \
            -DSILO_ENABLE_HDF5=ON \
            -DSILO_ENABLE_FORTRAN=OFF \
            -DSILO_ENABLE_BROWSER=OFF \
            -DSILO_ENABLE_SILOCK=OFF \
            -DSILO_ENABLE_ZFP=OFF \
            -DBUILD_TESTING=OFF \
            "-DCMAKE_PREFIX_PATH=$cmake_prefix_path" \
            "${native_arch_args[@]}" \
            "${compiler_args[@]}" "${hdf5_args[@]}" "${silo_hdf5_args[@]}"
        cmake --build "$silo_build" --parallel "$jobs"
        cmake --install "$silo_build"
        [[ -f "$silo_install/include/silo.h" ]] || { printf 'Silo header missing from %s\n' "$silo_install" >&2; exit 1; }
        silo_args=("-DSilo_ROOT=$silo_install")
        local_prefix_path="$silo_install${local_prefix_path:+;$local_prefix_path}"
        refresh_prefix_path
fi

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
printf 'HPX diagnostics: %s\n' "$hpx_debug"
# QueenBee4 reported ENOMEM while mapping HPX coroutine stacks.
cmake -S "$hpx_src" -B "$hpx_build" \
    "-DCMAKE_BUILD_TYPE=$build_type" \
    "-DCMAKE_INSTALL_PREFIX=$hpx_install" \
    -DHPX_WITH_CXX_STANDARD=20 \
    -DHPX_WITH_DISTRIBUTED_RUNTIME=ON \
    -DHPX_WITH_NETWORKING=ON \
    -DHPX_WITH_PARCELPORT_TCP=ON \
    -DHPX_WITH_THREAD_STACK_MMAP=OFF \
    "-DHPX_WITH_VERIFY_LOCKS=$hpx_debug" \
    "-DHPX_WITH_VERIFY_LOCKS_BACKTRACE=$hpx_debug" \
    -DHPX_WITH_THREAD_BACKTRACE_ON_SUSPENSION=OFF \
    "-DHPX_WITH_THREAD_DEBUG_INFO=$hpx_debug" \
    "-DHPX_WITH_SPINLOCK_DEADLOCK_DETECTION=$hpx_debug" \
    -DHPX_WITH_THREAD_DESCRIPTION_FULL=OFF \
    "-DHPX_WITH_THREAD_QUEUE_WAITTIME=$hpx_debug" \
    "-DHPX_WITH_THREAD_IDLE_RATES=$hpx_debug" \
    "-DHPX_WITH_THREAD_CREATION_AND_CLEANUP_RATES=$hpx_debug" \
    "-DHPX_WITH_THREAD_STEALING_COUNTS=$hpx_debug" \
    "-DHPX_WITH_COROUTINE_COUNTERS=$hpx_debug" \
    "-DHPX_WITH_PARCELPORT_COUNTERS=$hpx_debug" \
    "-DHPX_WITH_PARCELPORT_ACTION_COUNTERS=$hpx_debug" \
    -DHPX_WITH_STACKTRACES=ON \
    -DHPX_WITH_MALLOC=system \
    -DHPX_WITH_TESTS=OFF \
    -DHPX_WITH_EXAMPLES=OFF \
    -DHPX_WITH_APEX=ON \
    -DHPX_WITH_FETCH_APEX=ON \
    "-DHPX_WITH_FETCH_BOOST=$fetch_boost" \
    -DHPX_WITH_FETCH_ASIO=ON \
    "-DHPX_WITH_FETCH_HWLOC=$fetch_hwloc" \
    "-DCMAKE_PREFIX_PATH=$cmake_prefix_path" \
    "${native_arch_args[@]}" \
    "${compiler_args[@]}" "${dependency_args[@]}"
for option in HPX_WITH_VERIFY_LOCKS HPX_WITH_VERIFY_LOCKS_BACKTRACE \
              HPX_WITH_THREAD_DEBUG_INFO \
              HPX_WITH_SPINLOCK_DEADLOCK_DETECTION \
              HPX_WITH_THREAD_QUEUE_WAITTIME HPX_WITH_THREAD_IDLE_RATES \
              HPX_WITH_THREAD_CREATION_AND_CLEANUP_RATES HPX_WITH_THREAD_STEALING_COUNTS \
              HPX_WITH_COROUTINE_COUNTERS HPX_WITH_PARCELPORT_COUNTERS \
              HPX_WITH_PARCELPORT_ACTION_COUNTERS; do
    grep -Fx "$option:BOOL=$hpx_debug" "$hpx_build/CMakeCache.txt" >/dev/null || {
        printf 'HPX diagnostic option did not take effect: %s=%s\n' "$option" "$hpx_debug" >&2
        exit 1
    }
done
grep -Fx 'HPX_WITH_THREAD_BACKTRACE_ON_SUSPENSION:BOOL=OFF' "$hpx_build/CMakeCache.txt" >/dev/null
grep -Fx 'HPX_WITH_THREAD_DESCRIPTION_FULL:BOOL=OFF' "$hpx_build/CMakeCache.txt" >/dev/null
grep -Fx 'HPX_WITH_STACKTRACES:BOOL=ON' "$hpx_build/CMakeCache.txt" >/dev/null
grep -Fx 'HPX_WITH_THREAD_STACK_MMAP:BOOL=OFF' "$hpx_build/CMakeCache.txt" >/dev/null
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
    "${native_arch_args[@]}" \
    "${compiler_args[@]}" "${dependency_args[@]}" "${silo_args[@]}" "${hdf5_args[@]}"
cmake --build "$octo_build" --parallel "$jobs"
printf 'Executables: %s/octoII-{1d,2d,3d}; octoII -> octoII-3d\n' "$octo_build"
