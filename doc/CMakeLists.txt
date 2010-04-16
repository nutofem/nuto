# $Id$

get_target_property(DOC_TARGET doc TYPE)
if(NOT DOC_TARGET)
    add_custom_target(doc)
endif()

MESSAGE(STATUS "Checking for doxygen ...")
find_package(Doxygen)
if(DOXYGEN_FOUND)
    # configure Doxyfile
    set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/nuto)
    set(DOXYGEN_INPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/src/nuto)
    set(DOXYGEN_IMAGE_DIRECTORY ${CMAKE_SOURCE_DIR}/doc/images)
    set(DOXYGEN_HAVE_DOT NO)
    if(DOXYGEN_DOT_FOUND)
        set(DOXYGEN_HAVE_DOT YES)
    endif()
    #message(STATUS "NuTo_DOXYGEN_PREDEFINED = ${NuTo_DOXYGEN_PREDEFINED}")
    STRING(REPLACE ";" " " DOXYGEN_PREDEFINED "${NuTo_DOXYGEN_PREDEFINED}")
    #message(STATUS "DOXYGEN_PREDEFINED = ${DOXYGEN_PREDEFINED}")
    configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # clean doc path
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_CURRENT_BINARY_DIR}/nuto")

    # create target doxgen
    add_custom_target(doxygen ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # configure Doxyfile for python interface
    set(DOXYGEN_PYTHON_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/nuto_python)
    set(DOXYGEN_PYTHON_INPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../src/nuto)
    configure_file(DoxyfilePython.in ${CMAKE_CURRENT_BINARY_DIR}/DoxyfilePython)

    # clean doc path
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_CURRENT_BINARY_DIR}/nuto_python")

    # create target doxgen
    add_custom_target(doxygen_python ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/DoxyfilePython)

else()
    # create target doxgen
    add_custom_target(doxygen echo "WARNING: Doxygen not found -- source code documentation is not possible!")
    add_custom_target(doxygen_python echo "WARNING: Doxygen not found -- source code documentation is not possible!")
endif()

# add target doxygen as depedency of target doc
add_dependencies(doc doxygen doxygen_python)
