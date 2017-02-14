# NuTo specific cmake macros
#
# Stefan Eckardt, Institute of Structural Mechanics, Bauhaus-University Weimar, August 2009
#
# This file defines the following macros:
#
#   NUTO_INSTALL_PYTHON_FILE(python_file destination_dir)
#   NUTO_INSTALL_SWIG_PYTHON_MODULE(module_name destination_dir)


# macro for installing single python files (including byte compiling)
# name ... is the name of the python file
# destination ... is the directory in which the python file will be installed
MACRO(NUTO_INSTALL_PYTHON_FILE name destination)
  # get some informations about the python file
  GET_FILENAME_COMPONENT(NUTO_INSTALL_PYTHON_FILE_BASE_PATH ${name} PATH)
  GET_FILENAME_COMPONENT(NUTO_INSTALL_PYTHON_FILE_NAME ${name} NAME)
  GET_FILENAME_COMPONENT(NUTO_INSTALL_PYTHON_FILE_NAME_WE ${name} NAME_WE)
  # determine absolute and relative path
  SET(NUTO_INSTALL_PYTHON_FILE_ABSOLUTE_PATH "${CMAKE_CURRENT_BINARY_DIR}/${NUTO_INSTALL_PYTHON_FILE_BASE_PATH}")
  STRING(REGEX REPLACE "^${CMAKE_BINARY_DIR}/" "" NUTO_INSTALL_PYTHON_FILE_RELATIVE_PATH "${NUTO_INSTALL_PYTHON_FILE_ABSOLUTE_PATH}")
  # remove slash from end of the path
  STRING(REGEX REPLACE "/$" "" NUTO_INSTALL_PYTHON_FILE_RELATIVE_PATH ${NUTO_INSTALL_PYTHON_FILE_RELATIVE_PATH})
  STRING(REGEX REPLACE "/$" "" NUTO_INSTALL_PYTHON_FILE_ABSOLUTE_PATH ${NUTO_INSTALL_PYTHON_FILE_ABSOLUTE_PATH})
  #MESSAGE(STATUS "NUTO_INSTALL_PYTHON_FILE_ABSOLUTE_PATH: ${NUTO_INSTALL_PYTHON_FILE_ABSOLUTE_PATH}")
  #MESSAGE(STATUS "NUTO_INSTALL_PYTHON_FILE_RELATIVE_PATH: ${NUTO_INSTALL_PYTHON_FILE_RELATIVE_PATH}")
  # generate individual target name
  STRING(REPLACE "/" "_" NUTO_INSTALL_PYTHON_FILE_TARGET_NAME "${NUTO_INSTALL_PYTHON_FILE_RELATIVE_PATH}/${NUTO_INSTALL_PYTHON_FILE_NAME_WE}")
  # generate target
  ADD_CUSTOM_TARGET(nuto_python_byte_compile_${NUTO_INSTALL_PYTHON_FILE_TARGET_NAME} ALL
                    DEPENDS ${NUTO_INSTALL_PYTHON_FILE_NAME_WE}.pyc
                    WORKING_DIRECTORY ${NUTO_INSTALL_PYTHON_FILE_ABSOLUTE_PATH}
                   )
  # add rule for byte compilation
  ADD_CUSTOM_COMMAND(COMMAND ${PYTHON_EXECUTABLE}
                     ARGS -c "import py_compile;py_compile.compile(\"${name}\");"
                     DEPENDS ${NUTO_INSTALL_PYTHON_FILE_NAME}
                     OUTPUT ${NUTO_INSTALL_PYTHON_FILE_NAME_WE}.pyc
                     WORKING_DIRECTORY ${NUTO_INSTALL_PYTHON_FILE_ABSOLUTE_PATH}
                     VERBATIM
                     COMMENT "Byte compiling Python library: ${NUTO_INSTALL_PYTHON_FILE_RELATIVE_PATH}/${NUTO_INSTALL_PYTHON_FILE_NAME}"
                     )
  # install both files
  INSTALL(FILES ${NUTO_INSTALL_PYTHON_FILE_ABSOLUTE_PATH}/${NUTO_INSTALL_PYTHON_FILE_NAME} DESTINATION ${destination})
  INSTALL(FILES ${NUTO_INSTALL_PYTHON_FILE_ABSOLUTE_PATH}/${NUTO_INSTALL_PYTHON_FILE_NAME_WE}.pyc DESTINATION ${destination})
ENDMACRO(NUTO_INSTALL_PYTHON_FILE name destination)

# macro for installing a python module which was built using swig
# module_name ... is the name of the module
# destination ... is the directory in which the
MACRO(NUTO_INSTALL_SWIG_PYTHON_MODULE module_name destination)
  # get some informations about the python file
  GET_FILENAME_COMPONENT(NUTO_INSTALL_SWIG_PYTHON_MODULE_BASE_PATH ${module_name} PATH)
  GET_FILENAME_COMPONENT(NUTO_INSTALL_SWIG_PYTHON_MODULE_NAME ${module_name} NAME)
  # determine absolute and relative path
  SET(NUTO_INSTALL_SWIG_PYTHON_MODULE_ABSOLUTE_PATH "${CMAKE_CURRENT_BINARY_DIR}/${NUTO_INSTALL_SWIG_PYTHON_MODULE_BASE_PATH}")
  STRING(REGEX REPLACE "^${CMAKE_BINARY_DIR}/" "" NUTO_INSTALL_SWIG_PYTHON_MODULE_RELATIVE_PATH "${NUTO_INSTALL_SWIG_PYTHON_MODULE_ABSOLUTE_PATH}")
  # remove slash from end of the path
  STRING(REGEX REPLACE "/$" "" NUTO_INSTALL_SWIG_PYTHON_MODULE_RELATIVE_PATH "${NUTO_INSTALL_SWIG_PYTHON_MODULE_RELATIVE_PATH}")
  STRING(REGEX REPLACE "/$" "" NUTO_INSTALL_SWIG_PYTHON_MODULE_ABSOLUTE_PATH "${NUTO_INSTALL_SWIG_PYTHON_MODULE_ABSOLUTE_PATH}")
  #MESSAGE(STATUS "NUTO_INSTALL_SWIG_PYTHON_MODULE_ABSOLUTE_PATH: ${NUTO_INSTALL_SWIG_PYTHON_MODULE_ABSOLUTE_PATH}")
  #MESSAGE(STATUS "NUTO_INSTALL_SWIG_PYTHON_MODULE_RELATIVE_PATH: ${NUTO_INSTALL_SWIG_PYTHON_MODULE_RELATIVE_PATH}")
  # generate individual target name
  STRING(REPLACE "/" "_" NUTO_INSTALL_SWIG_PYTHON_MODULE_TARGET_NAME "${NUTO_INSTALL_SWIG_PYTHON_MODULE_RELATIVE_PATH}/${NUTO_INSTALL_SWIG_PYTHON_MODULE_NAME}")
  # generate target
  ADD_CUSTOM_TARGET(nuto_python_byte_compile_${NUTO_INSTALL_SWIG_PYTHON_MODULE_TARGET_NAME} ALL
                    DEPENDS ${NUTO_INSTALL_SWIG_PYTHON_MODULE_NAME}.pyc
                    WORKING_DIRECTORY ${NUTO_INSTALL_SWIG_PYTHON_MODULE_ABSOLUTE_PATH}
                   )
  # add rule for byte compilation
  ADD_CUSTOM_COMMAND(COMMAND ${PYTHON_EXECUTABLE}
                     ARGS -c "import py_compile;py_compile.compile(\"${NUTO_INSTALL_SWIG_PYTHON_MODULE_NAME}.py\");"
                     DEPENDS _${NUTO_INSTALL_SWIG_PYTHON_MODULE_NAME} ${NUTO_INSTALL_SWIG_PYTHON_MODULE_NAME}.py
                     OUTPUT ${NUTO_INSTALL_SWIG_PYTHON_MODULE_NAME}.pyc
                     WORKING_DIRECTORY ${NUTO_INSTALL_SWIG_PYTHON_MODULE_ABSOLUTE_PATH}
                     VERBATIM
                     COMMENT "Byte compiling Python library: ${NUTO_INSTALL_SWIG_PYTHON_MODULE_RELATIVE_PATH}/${NUTO_INSTALL_SWIG_PYTHON_MODULE_NAME}.py"
                     )
  # install all files
  INSTALL(TARGETS _${module_name} LIBRARY DESTINATION ${destination})
  INSTALL(FILES ${NUTO_INSTALL_SWIG_PYTHON_MODULE_ABSOLUTE_PATH}/${NUTO_INSTALL_SWIG_PYTHON_MODULE_NAME}.py DESTINATION ${destination})
  INSTALL(FILES ${NUTO_INSTALL_SWIG_PYTHON_MODULE_ABSOLUTE_PATH}/${NUTO_INSTALL_SWIG_PYTHON_MODULE_NAME}.pyc DESTINATION ${destination})
ENDMACRO(NUTO_INSTALL_SWIG_PYTHON_MODULE)

# NUTO_SWIG_MODULE(<module> <interface file> <install_path> <dependencies>)
#	Function to set up a SWIG module for a NuTo library.
#	<module>		SWIG module name
#	<interface file>	SWIG .i file
#	<dependencies>		Dependencies needed by that module. Typically the NuTo library it wraps.
#	<install_path>		Where to install the module. Relative to NUTO_PYTHON_MODULES_INSTALL_PATH.
function(NUTO_SWIG_MODULE module_name interface_file install_path)
  set(libraries ${ARGN})

  SET_SOURCE_FILES_PROPERTIES(${interface_file} PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(${interface_file} PROPERTIES SWIG_FLAGS "${NuTo_SWIG_FLAGS}")
  SWIG_ADD_MODULE(${module_name} python ${interface_file})
  # link library
  SWIG_LINK_LIBRARIES(${module_name} ${libraries} ${PYTHON_LIBRARIES})
  # check for unresolved symbols
  SET_TARGET_PROPERTIES(${SWIG_MODULE_${module_name}_REAL_NAME} PROPERTIES LINK_FLAGS -Wl,-z,defs)
  # Additional build flags for module
  SET_SOURCE_FILES_PROPERTIES("${swig_generated_file_fullname}" PROPERTIES COMPILE_FLAGS "${PYTHON_C_FLAGS}")

  # installation
  NUTO_INSTALL_SWIG_PYTHON_MODULE(${module_name} ${NUTO_PYTHON_MODULES_INSTALL_PATH}/${install_path})
endfunction()

# `add_unit_test(SomeClass FilesItNeedsToLink)` builds the unit test of the
# name `SomeClass` and all the files required (`FilesItNeedsToLink`), links
# NuToBase and the Boost unit test framework to it, and adds it to the test
# suite under an appropriate name
function(add_unit_test ClassName)
    # find relative path, i.e. remove the `.../test/`
    string(REPLACE "${CMAKE_SOURCE_DIR}/test/" "" relpath ${CMAKE_CURRENT_SOURCE_DIR})

    # if there are additional arguments, transform them from their filename
    # to the appropriate object library targets
    if(ARGN)
        foreach(filename ${ARGN})
            get_target_from_source(${CMAKE_SOURCE_DIR}/src/${filename} target)
            set(AdditionalObjects "${AdditionalObjects};$<TARGET_OBJECTS:${target}>")
        endforeach()
    endif()

    # look for the corresponding object to the unit test
    string(REPLACE "/" "." relpath ${relpath})
    set(srcObject "${relpath}.${ClassName}")
    # if there is an object (the class is not "header-only"), link it as well
    if(TARGET "${srcObject}")
        add_executable(${ClassName} ${ClassName}.cpp $<TARGET_OBJECTS:${srcObject}> ${AdditionalObjects})
    else()
        add_executable(${ClassName} ${ClassName}.cpp ${AdditionalObjects})
    endif()
    # link the unit test framework to the unit test
    target_link_libraries(${ClassName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} NuToBase)

    # generate a ctest name for the test
    string(REPLACE "." "::" testname ${relpath})
    add_test(unit::${testname}::${ClassName} ${CMAKE_CURRENT_BUILD_DIR}/${ClassName})
endfunction()

function(create_object_lib Source Object)
    # give a source file, and create an object lib for it
    # return the link target $<TARGET_OBJECT:target> to Object
    get_target_from_source(${Source} libName)
    add_library(${libName} OBJECT ${Source})

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
