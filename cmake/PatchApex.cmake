# HPX 1.11.0's pinned APEX gzstream.hpp uses uint64_t without its header.
# Repair the fetched copy after configuration, before compilation. Respect a
# cached FetchContent source override and leave already-fixed headers alone.
cmake_minimum_required(VERSION 3.18)
load_cache("${HPX_BINARY_DIR}" READ_WITH_PREFIX profiling_
  FETCHCONTENT_SOURCE_DIR_APEX FETCHCONTENT_BASE_DIR)
if(profiling_FETCHCONTENT_SOURCE_DIR_APEX)
  set(apex_source "${profiling_FETCHCONTENT_SOURCE_DIR_APEX}")
else()
  set(apex_source "${profiling_FETCHCONTENT_BASE_DIR}/apex-src")
endif()
set(header "${apex_source}/src/apex/gzstream.hpp")
if(NOT EXISTS "${header}")
  message(FATAL_ERROR "Cannot locate fetched APEX header: ${header}")
endif()
file(READ "${header}" contents)
if(NOT contents MATCHES "#include [<\"](cstdint|stdint.h)[>\"]")
  if(NOT contents MATCHES "#include <streambuf>")
    message(FATAL_ERROR "APEX gzstream.hpp changed; review the integer-header fix")
  endif()
  string(REPLACE "#include <streambuf>" "#include <cstdint>\n#include <streambuf>" contents "${contents}")
  file(WRITE "${header}" "${contents}")
  message(STATUS "Added missing <cstdint> include to HPX's fetched APEX")
endif()
