cmake_minimum_required(VERSION 3.22)
project(OctotigerII VERSION 0.1.0 LANGUAGES C CXX)

option(OCTOTIGERII_WITH_HPX "Use distributed HPX field storage and locality schedulers" ON)
option(OCTOTIGERII_WITH_PROFILING "Name HPX tasks and enable APEX regions when available" ON)
option(OCTOTIGERII_BUILD_TESTS "Build GoogleTest unit/integration suites and application checks" ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
include(GNUInstallDirs)
foreach(module HYDRO RADIATION GRAVITY MASS_FRACTIONS)
  option(OCTOII_WITH_${module} "Compile the ${module} physics module" ON)
endforeach()
if(DEFINED OCTOTIGERII_PROBLEM OR DEFINED OCTOTIGERII_NDIM)
  message(FATAL_ERROR "Remove OCTOTIGERII_PROBLEM/OCTOTIGERII_NDIM from your cache. All three dimensions are built; select --problem.name at runtime.")
endif()
include(cmake/Problems.cmake)
find_package(Threads REQUIRED)
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()
find_package(Boost 1.71 REQUIRED)


if(OCTOTIGERII_BUILD_TESTS)
  enable_testing()
endif()
foreach(OCTOTIGERII_NDIM RANGE 1 3)
  add_subdirectory(cmake/dimension "${OCTOTIGERII_NDIM}d")
endforeach()
add_custom_target(octoII ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink octoII-3d "${PROJECT_BINARY_DIR}/octoII"
  DEPENDS octoII-3d VERBATIM)
install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink octoII-3d \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/octoII\")")
install(DIRECTORY bin/ DESTINATION ${CMAKE_INSTALL_DATADIR}/OctotigerII
  FILES_MATCHING PATTERN inputs)

# Documentation stays optional for the application build. The standalone
# docs.sh path needs neither configured physics dependencies nor HPX.
add_custom_target(docs
  COMMAND ${PROJECT_SOURCE_DIR}/docs.sh --output ${CMAKE_CURRENT_BINARY_DIR}/docs
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMENT "Generating the OctotigerII developer reference"
  VERBATIM)
