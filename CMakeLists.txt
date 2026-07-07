cmake_minimum_required(VERSION 3.15)
project(fesim LANGUAGES C CXX)

# C++11 required by all targets
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- BLAS detection ---
# Intel MKL (HPC): module load intel sets MKLROOT
# macOS: Apple Accelerate framework
# Other: user specifies BLA_VENDOR or BLAS_LIBRARIES
if(DEFINED ENV{MKLROOT})
    set(BLA_VENDOR "Intel10_64lp")
elseif(APPLE)
    set(BLA_VENDOR "Apple")
endif()
find_package(BLAS REQUIRED)

# Determine which BLAS header to use based on the library CMake found,
# NOT on compiler identity (GCC+MKL should use <mkl.h>, not <cblas.h>)
if(BLA_VENDOR MATCHES "^Intel" OR BLAS_LIBRARIES MATCHES "mkl")
    set(BLAS_HEADER_DEFINE "BLAS_MKL")
elseif(BLA_VENDOR STREQUAL "Apple" OR BLAS_LIBRARIES MATCHES "Accelerate")
    set(BLAS_HEADER_DEFINE "BLAS_ACCELERATE")
else()
    set(BLAS_HEADER_DEFINE "BLAS_GENERIC")
endif()
message(STATUS "BLAS header: ${BLAS_HEADER_DEFINE}")

# --- MKL include directory ---
# FindBLAS often finds MKL libraries but not its include dir.
# Try MKLROOT first (set by 'module load intel'), then infer from library path.
if(BLAS_HEADER_DEFINE STREQUAL "BLAS_MKL")
    if(NOT BLAS_INCLUDE_DIRS)
        # 1) $MKLROOT/include (if module loaded)
        if(DEFINED ENV{MKLROOT} AND IS_DIRECTORY "$ENV{MKLROOT}/include")
            set(BLAS_INCLUDE_DIRS "$ENV{MKLROOT}/include")
        else()
            # 2) Infer from library path: .../mkl/lib/intel64_lin/... -> .../mkl/include
            list(GET BLAS_LIBRARIES 0 _first_blas_lib)
            string(REGEX REPLACE "/lib/.*" "/include" _blas_include_guess "${_first_blas_lib}")
            if(IS_DIRECTORY "${_blas_include_guess}")
                set(BLAS_INCLUDE_DIRS "${_blas_include_guess}")
            else()
                # 3) Search common locations
                find_path(BLAS_INCLUDE_DIRS mkl.h
                    PATHS
                        "$ENV{MKLROOT}/include"
                        /opt/intel/mkl/include
                        /opt/intel/oneapi/mkl/latest/include
                )
            endif()
        endif()
    endif()
    if(BLAS_INCLUDE_DIRS)
        message(STATUS "MKL include: ${BLAS_INCLUDE_DIRS}")
    else()
        message(WARNING "BLAS=MKL but could not find mkl.h include directory. "
                        "Set MKLROOT or pass -DBLAS_INCLUDE_DIRS=/path/to/mkl/include")
    endif()
endif()

# --- inih (INI config parser) ---
add_library(inih STATIC
    external/inih/ini.c
    external/inih/cpp/INIReader.cpp
)
target_include_directories(inih PUBLIC
    external/inih
    external/inih/cpp
)

# --- subdirectories ---
add_subdirectory(fesim)
add_subdirectory(qmat)
