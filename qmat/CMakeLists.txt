find_package(MPI REQUIRED)

add_executable(qmat main.cpp)
target_link_libraries(qmat PRIVATE MPI::MPI_CXX ${BLAS_LIBRARIES} inih)
target_include_directories(qmat PRIVATE ${BLAS_INCLUDE_DIRS})
target_compile_definitions(qmat PRIVATE ${BLAS_HEADER_DEFINE})
