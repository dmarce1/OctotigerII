cmake_minimum_required(VERSION 3.22)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(OctotigerIIHelmholtz VERSION 0.1 LANGUAGES CXX)
  include(CTest)
endif()
add_library(octotigerII_helmholtz STATIC helmholtz.cpp table.cpp direct.cpp)
add_library(OctotigerII::helmholtz ALIAS octotigerII_helmholtz)
target_compile_features(octotigerII_helmholtz PUBLIC cxx_std_17)
set_target_properties(octotigerII_helmholtz PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(octotigerII_helmholtz PRIVATE -Wall -Wextra -Wpedantic)
endif()
target_include_directories(octotigerII_helmholtz PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
# Pin the supplied data; never substitute a locally generated table.
set(timmes_table "${CMAKE_CURRENT_SOURCE_DIR}/data/helm_table.dat")
file(SHA256 "${timmes_table}" timmes_table_sha256)
if(NOT timmes_table_sha256 STREQUAL "c9a57c26c6fd2b2b378b9d5295ca1214022f6fec6289d038b47bf8c8938881a1")
  message(FATAL_ERROR "Timmes supplied Helmholtz table checksum mismatch")
endif()
target_compile_definitions(octotigerII_helmholtz PRIVATE
  OCTOTIGERII_HELMHOLTZ_TABLE="${timmes_table}")
# Research tool only; omitted from normal builds.
add_executable(helmholtz_build_table EXCLUDE_FROM_ALL build_table.cpp)
target_link_libraries(helmholtz_build_table PRIVATE octotigerII_helmholtz)
# Explicit target: full table generation is never hidden in normal builds.
add_custom_target(helmholtz_table
  COMMAND helmholtz_build_table ${CMAKE_CURRENT_BINARY_DIR}/helm_table.experimental.dat
  DEPENDS helmholtz_build_table
  COMMENT "Generating Helmholtz table from Timmes Fermi-Dirac quadrature"
  VERBATIM)
if(BUILD_TESTING OR OCTOTIGERII_BUILD_TESTS)
  add_executable(helmholtz_checks tests/checks.cpp)
  target_link_libraries(helmholtz_checks PRIVATE octotigerII_helmholtz)
  target_include_directories(helmholtz_checks PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME helmholtz_checks COMMAND helmholtz_checks)
  add_executable(helmholtz_probe tests/probe.cpp)
  target_link_libraries(helmholtz_probe PRIVATE octotigerII_helmholtz)
  target_include_directories(helmholtz_probe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()
