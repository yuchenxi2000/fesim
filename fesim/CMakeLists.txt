# Apple Clang doesn't ship OpenMP; Homebrew libomp is at a known path
if(APPLE AND NOT OpenMP_ROOT)
    set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
endif()
find_package(OpenMP REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX ${BLAS_LIBRARIES} inih)
target_include_directories(main PRIVATE ${BLAS_INCLUDE_DIRS})
target_compile_definitions(main PRIVATE ${BLAS_HEADER_DEFINE})

# CMake's FindOpenMP may omit the include dir with Apple Clang
if(APPLE AND OpenMP_CXX_FOUND AND NOT OpenMP_CXX_INCLUDE_DIRS)
    target_include_directories(main PRIVATE "${OpenMP_ROOT}/include")
endif()
