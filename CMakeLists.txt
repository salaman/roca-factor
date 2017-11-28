cmake_minimum_required(VERSION 3.5)
project(roca_factor)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3 -ffast-math -DNDEBUG")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_library(GMP_LIBRARY gmp REQUIRED)
find_path(GMP_INCLUDE_DIR gmp.h)

find_library(FPLLL_LIBRARY fplll REQUIRED)
find_path(FPLLL_INCLUDE_DIR fplll.h)

find_library(FLINT_LIBRARY flint REQUIRED)
find_path(FLINT_INCLUDE_DIR flint.h PATH_SUFFIXES flint)

find_library(MPFR_LIBRARY mpfr REQUIRED)

add_executable(roca_factor main.cpp DiscreteLog.cpp DiscreteLog.h)

target_link_libraries(roca_factor ${GMP_LIBRARY} ${GMPXX_LIBRARY} ${FPLLL_LIBRARY} ${FLINT_LIBRARY} ${MPFR_LIBRARY})
target_include_directories(roca_factor PUBLIC ${GMP_INCLUDE_DIR} ${GMPXX_INCLUDE_DIR} ${FPLLL_INCLUDE_DIR} ${FLINT_INCLUDE_DIR})
