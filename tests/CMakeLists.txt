# Use an installed GoogleTest first (including module/package-manager prefixes).
# The pinned fallback requires no root access. Offline users can disable fetching
# or supply FETCHCONTENT_SOURCE_DIR_GOOGLETEST to an unpacked source tree.
option(OCTOTIGERII_FETCH_GOOGLETEST "Fetch GoogleTest if no installation is found" ON)
find_package(GTest 1.12 QUIET)
if(NOT GTest_FOUND)
  if(NOT OCTOTIGERII_FETCH_GOOGLETEST)
    message(FATAL_ERROR "GoogleTest 1.12+ is required for tests. Set GTest_ROOT/CMAKE_PREFIX_PATH, enable OCTOTIGERII_FETCH_GOOGLETEST, or disable OCTOTIGERII_BUILD_TESTS.")
  endif()
  include(FetchContent)
  if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
  endif()
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
    URL_HASH SHA256=8ad598c73ad796e0d8280b082cebd82a630d73e73cd3c70057938a6501bba5d7)
  FetchContent_MakeAvailable(googletest)
endif()
include(GoogleTest)

set(testRuntimeArgs)
if(OCTOTIGERII_WITH_HPX)
  set(testRuntimeArgs --hpx:threads=2 --hpx:bind=none)
endif()

function(octo_test target library labels)
  add_executable(${target} ${target}.cpp googleTestMain.cpp)
  target_link_libraries(${target} PRIVATE ${library} GTest::gtest)
  # Runtime flags and compile definitions come from Application only for tests
  # that need distributed storage. List-tests always runs without HPX.
  gtest_discover_tests(${target}
    TEST_PREFIX "${OCTOTIGERII_PROBLEM}.${OCTOTIGERII_NDIM}d."
    EXTRA_ARGS ${testRuntimeArgs}
    DISCOVERY_MODE PRE_TEST
    PROPERTIES TIMEOUT 300 LABELS "${labels}")
endfunction()

octo_test(unitsChecks OctotigerII::Subgrid "unit")
octo_test(meshChecks OctotigerII::Mesh "unit")
octo_test(mathChecks OctotigerII::Mesh "unit")
octo_test(optionsChecks OctotigerII::Subgrid "unit")
octo_test(boundaryChecks OctotigerII::Application "unit;integration")
octo_test(limiterChecks OctotigerII::Physics "unit")
if(PROBLEM_HYDRO)
  octo_test(hydroChecks OctotigerII::Hydro "unit")
endif()
if(PROBLEM_RADIATION)
  octo_test(radiationChecks OctotigerII::Radiation "unit")
endif()
if(PROBLEM_HYDRO OR PROBLEM_RADIATION)
  octo_test(finiteVolumeChecks OctotigerII::Application "integration")
endif()
if(PROBLEM_GRAVITY)
  octo_test(gravityChecks OctotigerII::Gravity "unit")
  octo_test(ewaldChecks OctotigerII::Gravity "unit;gravity")
  octo_test(directGravityChecks OctotigerII::Application "integration")
  # Has a custom main to preserve the existing 'benchmark [level] [order]' CLI.
  add_executable(gravityParallelChecks gravityParallelChecks.cpp)
  target_link_libraries(gravityParallelChecks PRIVATE OctotigerII::Application GTest::gtest)
  gtest_discover_tests(gravityParallelChecks
    TEST_PREFIX "${OCTOTIGERII_PROBLEM}.${OCTOTIGERII_NDIM}d."
    EXTRA_ARGS ${testRuntimeArgs} DISCOVERY_MODE PRE_TEST
    PROPERTIES TIMEOUT 300 LABELS "integration")
endif()
octo_test(numericalChecks OctotigerII::Application "integration")
if(OCTOTIGERII_PROBLEM STREQUAL "rayleigh-taylor")
  octo_test(rayleighTaylorChecks OctotigerII::Application "integration;rayleigh-taylor")
endif()
if(OCTOTIGERII_PROBLEM STREQUAL "streaming" AND OCTOTIGERII_NDIM STREQUAL "1")
  target_compile_definitions(numericalChecks PRIVATE OCTOTIGERII_TEST_STREAMING_PROFILE=1)
endif()
octo_test(storageChecks OctotigerII::Application "integration")
octo_test(amrChecks OctotigerII::Application "unit;integration;amr")
octo_test(analyticChecks OctotigerII::Application "integration")
if(OCTOTIGERII_PROBLEM MATCHES "^(sod|gravity-sphere|gravity-gaussian|streaming)$")
  target_compile_definitions(analyticChecks PRIVATE OCTOTIGERII_TEST_CONVERGENCE=1)
endif()
octo_test(siloChecks OctotigerII::Application "integration")
target_link_libraries(siloChecks PRIVATE ${Silo_LIBRARY})
target_include_directories(siloChecks PRIVATE ${Silo_INCLUDE_DIR})
if(OCTOTIGERII_WITH_HPX)
  octo_test(serializationChecks OctotigerII::Application "unit")
endif()

# Whole-process CLI, report, and selected-problem smoke checks complement the
# GoogleTest suites and remain Python/CTest tests.
add_test(NAME example.${OCTOTIGERII_PROBLEM}.${OCTOTIGERII_NDIM}d COMMAND octotigerII
  --config=${PROBLEM_DIRECTORY}/inputs --mesh.cells=4 --mesh.level=1 --output.enabled=off ${testRuntimeArgs})
set_tests_properties(example.${OCTOTIGERII_PROBLEM}.${OCTOTIGERII_NDIM}d PROPERTIES TIMEOUT 120 LABELS "application")
if(OCTOTIGERII_PROBLEM STREQUAL "rayleigh-taylor")
  set_tests_properties(example.${OCTOTIGERII_PROBLEM}.${OCTOTIGERII_NDIM}d PROPERTIES TIMEOUT 180)
endif()
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME application COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/application.py
    $<TARGET_FILE:octotigerII> ${PROBLEM_DIRECTORY}/inputs ${CMAKE_CURRENT_BINARY_DIR}/application-output
    ${OCTOTIGERII_PROBLEM} ${OCTOTIGERII_NDIM} ${testRuntimeArgs})
  set_tests_properties(application PROPERTIES TIMEOUT 180 LABELS "application")
  if(OCTOTIGERII_WITH_HPX)
    foreach(count 2 3)
      foreach(target storageChecks numericalChecks boundaryChecks amrChecks)
        add_test(NAME distributed.${target}.${count} COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/distributed.py $<TARGET_FILE:${target}>)
        set_tests_properties(distributed.${target}.${count} PROPERTIES
          ENVIRONMENT "OCTOTIGERII_TEST_LOCALITIES=${count}" TIMEOUT 360 LABELS "distributed" RUN_SERIAL TRUE)
      endforeach()
      if(PROBLEM_GRAVITY)
        add_test(NAME gravity.distributed.${count} COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/distributed.py $<TARGET_FILE:gravityParallelChecks>)
        set_tests_properties(gravity.distributed.${count} PROPERTIES
          ENVIRONMENT "OCTOTIGERII_TEST_LOCALITIES=${count}" TIMEOUT 360 LABELS "distributed;gravity" RUN_SERIAL TRUE)
      endif()
      if(OCTOTIGERII_PROBLEM STREQUAL "rayleigh-taylor")
        add_test(NAME rayleigh-taylor.distributed.${count} COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/distributed.py $<TARGET_FILE:rayleighTaylorChecks>
          --gtest_filter=-RayleighTaylor.PerturbationGrowsUnderGravity)
        set_tests_properties(rayleigh-taylor.distributed.${count} PROPERTIES
          ENVIRONMENT "OCTOTIGERII_TEST_LOCALITIES=${count}" TIMEOUT 360 LABELS "distributed;rayleigh-taylor" RUN_SERIAL TRUE)
      endif()
    endforeach()
  endif()
endif()
