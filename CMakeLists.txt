cmake_minimum_required(VERSION 3.15)
project(SpBench LANGUAGES CXX CUDA)
project(SpBench C CXX)

# Exposed options. Disable something, if do not want to build cuda or cl stuff.
option(BENCH_WITH_CUBOOL      "Add cubool lib and related benchmarks" OFF)
option(BENCH_WITH_CUSP        "Add cusp lib and related benchmarks" OFF)
option(BENCH_WITH_CLSPARSE    "Add clSPARSE lib and related benchmarks" OFF)
option(BENCH_WITH_CUSPARSE    "Add cuSPARSE lib and related benchmarks" OFF)
option(BENCH_WITH_CLBOOL      "Add clbool lib and related benchmarks" OFF)
option(BENCH_WITH_SUITESPARSE "Add GraphBLAS:SuiteSparse lib and related benchmarks" ON)

add_library(sp_bench_base INTERFACE)
target_include_directories(sp_bench_base INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src)

# Append here all benchmark targets
set(TARGETS)

if (BENCH_WITH_SUITESPARSE)
    set(SUITESPARSE_TARGETS)

    add_executable(suitesparse_multiply src/suitesparse_multiply.cpp)
    add_executable(suitesparse_msbfs src/suitesparse_msbfs.cpp)
    list(APPEND SUITESPARSE_TARGETS suitesparse_multiply suitesparse_msbfs)

    foreach(SUITESPARSE_TARGET ${SUITESPARSE_TARGETS})
        target_link_libraries(${SUITESPARSE_TARGET} PUBLIC sp_bench_base)
        target_link_libraries(${SUITESPARSE_TARGET} PUBLIC lagraph)
        target_link_libraries(${SUITESPARSE_TARGET} PUBLIC graphblas)

        # The graphblas is installed here, must explicitly add path to be sure, that link will be successful
        target_link_directories(${SUITESPARSE_TARGET} PUBLIC /usr/local/lib)
        target_compile_features(${SUITESPARSE_TARGET} PUBLIC cxx_std_14)
        set_target_properties(${SUITESPARSE_TARGET} PROPERTIES CXX_STANDARD 17)
        set_target_properties(${SUITESPARSE_TARGET} PROPERTIES CXX_STANDARD_REQUIRED ON)

        list(APPEND TARGETS ${SUITESPARSE_TARGET})
    endforeach()
endif()

# Some fancy stuff here
foreach(TARGET ${TARGETS})
    message(STATUS "Build target benchmark ${TARGET}")
endforeach()

# Copy data into build directory
file(COPY ${CMAKE_CURRENT_LIST_DIR}/data/ DESTINATION ${CMAKE_BINARY_DIR}/data/)

# Copy scripts into build directory
file(GLOB SCRIPTS "*.sh")
file(COPY ${SCRIPTS} DESTINATION ${CMAKE_BINARY_DIR})
