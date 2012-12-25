cmake_minimum_required (VERSION 2.6)
project (rogsolv)

set (CMAKE_VERBOSE_MAKEFILE ON)
set (CMAKE_CXX_COMPILER mpic++.openmpi)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -pedantic -Wextra -Winline -Wdisabled-optimization -mfpmath=sse -msse2 --param inline-unit-growth=100")
# -Wunreachable-code -Wcast-qual -Wcast-align -Wconversion

aux_source_directory (src library_sources)

add_executable (rogsolv ${library_sources})
target_link_libraries (rogsolv rt boost_filesystem boost_system)

