cmake_minimum_required (VERSION 2.6)
project (rogsolv)

set (CMAKE_VERBOSE_MAKEFILE ON)
set (CMAKE_CXX_COMPILER mpic++.openmpi)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -mfpmath=sse -msse2 -pedantic -Wall -Wextra -Winline -Wdisabled-optimization --param inline-unit-growth=100 -Wfloat-equal -Wunreachable-code -Wcast-qual -Wcast-align -Wconversion")
# -Wfloat-equal -Wunreachable-code -Wcast-qual -Wcast-align -Wconversion

aux_source_directory (src library_sources)
aux_source_directory (src/base library_sources)
aux_source_directory (src/ci library_sources)
aux_source_directory (src/ci_griders library_sources)
aux_source_directory (src/containers library_sources)
aux_source_directory (src/schemes library_sources)
aux_source_directory (src/workers library_sources)
aux_source_directory (src/writers library_sources)

add_executable (rogsolv ${library_sources})
target_link_libraries (rogsolv boost_filesystem-mt boost_system-mt boost_chrono-mt)

