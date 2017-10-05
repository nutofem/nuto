# nuto_swig_module(<module> <interface file> <install_path> <dependencies>)
#   Function to set up a SWIG module for a NuTo library.
#   <module>        SWIG module name
#   <interface file>    SWIG .i file
#   <dependencies>      Dependencies needed by that module.
#                       Typically the NuTo library it wraps.
function(nuto_swig_module module_name interface_file)
    set(libraries ${ARGN})

    set_source_files_properties(${interface_file} PROPERTIES CPLUSPLUS ON)
    set_source_files_properties(${interface_file}
        PROPERTIES SWIG_FLAGS "${NuTo_SWIG_FLAGS}")
    swig_add_module(${module_name} python ${interface_file})
    # link library
    swig_link_libraries(${module_name} ${libraries} ${PYTHON_LIBRARIES})
    # check for unresolved symbols
    set_target_properties(${SWIG_MODULE_${module_name}_REAL_NAME}
        PROPERTIES LINK_FLAGS -Wl,-z,defs)
    # Additional build flags for module
    set_source_files_properties("${swig_generated_file_fullname}"
        PROPERTIES COMPILE_FLAGS "${PYTHON_CXX_FLAGS}")
endfunction()

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
        add_executable(${ClassName} EXCLUDE_FROM_ALL ${ClassName}.cpp
            $<TARGET_OBJECTS:${srcObject}> ${AdditionalObjects})
    else()
        add_executable(${ClassName} EXCLUDE_FROM_ALL ${ClassName}.cpp
            ${AdditionalObjects})
    endif()
    # link the unit test framework to the unit test
    target_link_libraries(${ClassName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    target_include_directories(${ClassName} PUBLIC ${CMAKE_SOURCE_DIR}/src)

    # generate a ctest name for the test
    string(REPLACE "." "::" testname ${relpath})
    add_test(unit::${testname}::${ClassName} ${ClassName} --log_level=message)
    set(all_unit_tests "${all_unit_tests};${ClassName}"
        CACHE INTERNAL "The names of all the unit tests")
endfunction()

function(create_object_lib Source Object)
    # give a source file, and create an object lib for it
    # return the link target $<TARGET_OBJECT:target> to Object
    get_target_from_source(${Source} libName)
    add_library(${libName} OBJECT ${Source})
    target_include_directories(${libName} PUBLIC ${CMAKE_SOURCE_DIR}/src)

    set(objectLib "$<TARGET_OBJECTS:${libName}>")
    set(${Object} ${objectLib} PARENT_SCOPE)
endfunction()

# in three steps go from "/home/user/nuto/src/mechanics/MyClass.cpp"
# to "mechanics.MyClass", and export this to Target
function(get_target_from_source Source Target)
    string(REPLACE "${CMAKE_SOURCE_DIR}/src/" "" libName ${Source})
    string(REPLACE ".cpp" "" libName ${libName})
    string(REPLACE "/" "." libName ${libName})
    set(${Target} ${libName} PARENT_SCOPE)
endfunction()

# give it a list of sources, and it creates an object lib for each source,
# and returns the list of link targets
function(sources_to_objects Sources Objects)
    foreach(source ${Sources})
        create_object_lib(${CMAKE_CURRENT_SOURCE_DIR}/${source} object)
        list(APPEND objectList ${object})
    endforeach()
    set(${Objects} ${objectList} PARENT_SCOPE)
endfunction()


# symlink a given file/directory to the build directory
function(create_symlink pathName)
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
        "${CMAKE_CURRENT_SOURCE_DIR}/${pathName}"
        "${CMAKE_CURRENT_BINARY_DIR}/${pathName}")
endfunction()

function(create_nuto_module ModuleName ModuleSources)
    sources_to_objects("${ModuleSources}" ModuleObjects)

    add_library(${ModuleName} ${ModuleObjects} ${ARGN})
    set_target_properties(${ModuleName} PROPERTIES
        OUTPUT_NAME NuTo${ModuleName}
        )
    target_include_directories(${ModuleName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include/nuto>
        )

    install(TARGETS ${ModuleName} EXPORT NuToTargets
        LIBRARY DESTINATION lib
        INCLUDES DESTINATION include/nuto
        )
endfunction()

function(warning)
    string(ASCII 27 Esc)
    set(ColourReset "${Esc}[m")
    set(Red         "${Esc}[31m")
    message(STATUS "${Red}${ARGV}${ColourReset}")
endfunction()

function(append_to_tests TestName)
    set(all_integration_tests "${all_integration_tests};${TestName}"
        CACHE INTERNAL "The names of all the integration tests")
endfunction()

function(append_to_benchmarks Benchmark)
    set(all_benchmarks "${all_benchmarks};${Benchmark}"
        CACHE INTERNAL "The names of all the benchmarks")
endfunction()
