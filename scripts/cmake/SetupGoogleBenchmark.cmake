function(setup_google_benchmark)
        execute_process(COMMAND "./scripts/setup_google_benchmark.sh" WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
        add_subdirectory(external/benchmark)
endfunction()
