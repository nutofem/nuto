# Creates a benchmark and links all nuto-libraries and the benchmark lib to it.
# Special compiler setup is used to avoid warnings from the benchmark lib
function(add_benchmark BenchmarkName)
    if(benchmark_FOUND)
        add_executable(${BenchmarkName} EXCLUDE_FROM_ALL ${BenchmarkName}.cpp)
        target_include_directories(${BenchmarkName} PRIVATE ${CMAKE_SOURCE_DIR}/src)
        target_link_libraries(${BenchmarkName} benchmark Mechanics Math Base Visualize)
        append_to_benchmarks(${BenchmarkName})
    endif()
endfunction()
