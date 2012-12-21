cmake_minimum_required (VERSION 2.6)
project (rogsolv)

set (CMAKE_VERBOSE_MAKEFILE ON)
set (CMAKE_CXX_COMPILER /usr/bin/mpic++.openmpi)
set (CMAKE_CXX_FLAGS "-std=c++0x -Wall -pedantic -Wextra -Winline -Wdisabled-optimization -mfpmath=sse -msse2 --param inline-unit-growth=100")
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -g0 -DNDEBUG") #-DNDEBUG
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")
set (CMAKE_EXE_LINKER_FLAGS "-lrt -lboost_system -lboost_filesystem")
# -Wunreachable-code -Wcast-qual -Wcast-align -Wconversion
# -lboost_filesystem-mt

set (MPI_LINK_FLAGS "-lrt")

aux_source_directory (src library_sources)

add_executable (rogsolv ${library_sources})

