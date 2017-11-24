# Creates a benchmark and links all nuto-libraries and the benchmark lib to it.
# Special compiler setup is used to avoid warnings from the benchmark lib
function(add_benchmark BenchmarkName)
    add_executable(${BenchmarkName} EXCLUDE_FROM_ALL ${BenchmarkName}.cpp)
    target_include_directories(${BenchmarkName} PRIVATE ${CMAKE_SOURCE_DIR}/src)
    target_link_libraries(${BenchmarkName} benchmark Mechanics Math Base Visualize)
    # Disable missing override warnings
    set_target_properties(${BenchmarkName} PROPERTIES COMPILE_FLAGS "-Wno-suggest-override")
endfunction()
