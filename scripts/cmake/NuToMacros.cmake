include(AddBenchmark)
include(AddUnitTest)
include(AddIntegrationtest)
include(CheckForDependencies)
include(CreateSymlink)
include(CreateNutoModule)
include(GetGitRevisionDescription)
include(NutoSwigModule)
include(SetCompilerFlags)
include(Warning)

function(append_to_benchmarks Benchmark)
    set(all_benchmarks "${all_benchmarks};${Benchmark}"
        CACHE INTERNAL "The names of all the benchmarks")
endfunction()

function(append_to_examples Example)
    set(all_examples "${all_examples};${Example}"
        CACHE INTERNAL "The names of all the examples")
endfunction()

# get commit hash in base64
add_custom_target(turtles
    COMMAND python3 -c "import base64;print(base64.b64decode('G1sxQRtbMksbWzMybUxlb25hcmRvLCBSYXBoYWVsLCBEb25hdGVsbG8gYW5kIE1pY2hlbGFuZ28=').decode())" VERBATIM)
