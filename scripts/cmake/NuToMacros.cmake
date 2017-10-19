include(AddUnitTest)
include(AddIntegrationtest)
include(CheckForDependencies)
include(CreateObjectLib)
include(CreateSymlink)
include(CreateNutoModule)
include(GetGitRevisionDescription)
include(GetTargetFromSource)
include(NutoSwigModule)
include(SetCompilerFlags)
include(SourcesToObjects)
include(Warning)

function(append_to_benchmarks Benchmark)
    set(all_benchmarks "${all_benchmarks};${Benchmark}"
        CACHE INTERNAL "The names of all the benchmarks")
endfunction()

function(append_to_examples Example)
    set(all_examples "${all_examples};${Example}"
        CACHE INTERNAL "The names of all the examples")
endfunction()
