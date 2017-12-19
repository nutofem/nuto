# `add_unit_test(SomeClass FilesItNeedsToLink)` builds the unit test of the
# name `SomeClass` and all the files required (`FilesItNeedsToLink`), links
# Boost unit test framework to it, and adds it to the test suite under an
# appropriate name
function(add_unit_test ClassName)
    # find relative path, i.e. remove the `.../test/`
    string(REPLACE "${CMAKE_SOURCE_DIR}/test/" ""
        relpath ${CMAKE_CURRENT_SOURCE_DIR})

    # if there are additional arguments, transform them from their filename
    # to the appropriate object library targets
    if(ARGN)
        foreach(filename ${ARGN})
            get_target_from_source(${CMAKE_SOURCE_DIR}/src/${filename} target)
            set(AdditionalObjects
                "${AdditionalObjects};$<TARGET_OBJECTS:${target}>")
        endforeach()
    endif()

    # look for the corresponding object to the unit test
    string(REPLACE "/" "." relpath ${relpath})
    set(srcObject "${relpath}.${ClassName}")
    # if there is an object (the class is not "header-only"), link it as well
    if(TARGET "${srcObject}")
        add_executable(${ClassName} ${ClassName}.cpp
            $<TARGET_OBJECTS:${srcObject}> ${AdditionalObjects})
    else()
        add_executable(${ClassName} ${ClassName}.cpp
            ${AdditionalObjects})
    endif()
    # link the unit test framework to the unit test
    target_link_libraries(${ClassName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    target_include_directories(${ClassName} PRIVATE ${CMAKE_SOURCE_DIR}/src)
    target_include_directories(${ClassName}
        PRIVATE ${CMAKE_SOURCE_DIR}/test/tools)
    target_link_libraries(${ClassName} Fakeit)

    # generate a ctest name for the test
    string(REPLACE "." "::" testname ${relpath})
    add_test(unit::${testname}::${ClassName} ${ClassName} --log_level=message)
    set(all_unit_tests "${all_unit_tests};${ClassName}"
        CACHE INTERNAL "The names of all the unit tests")
endfunction()
