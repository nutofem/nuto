# save output from `git rev-parse HEAD` to variable `commit_nr`
execute_process(COMMAND
    git rev-parse HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE commit_nr)

# remove trailing newline
string(REGEX REPLACE "\n$" "" commit_nr "${commit_nr}")
