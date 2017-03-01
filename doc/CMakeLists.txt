message(STATUS "Checking for doxygen ...")
find_package(Doxygen)
if(DOXYGEN_FOUND)

    # configure Doxyfile
    set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/nuto)
    set(DOXYGEN_INPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/examples ${CMAKE_SOURCE_DIR}/doc")
    set(DOXYGEN_IMAGES_DIRECTORY "${CMAKE_SOURCE_DIR}/doc/images")

    set(DOXYGEN_HAVE_DOT NO)
    if(DOXYGEN_DOT_FOUND)
        set(DOXYGEN_HAVE_DOT YES)
    endif()

    string(REPLACE ";" " " DOXYGEN_PREDEFINED "${NuTo_DOXYGEN_PREDEFINED}")

    configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # clean doc path
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_CURRENT_BINARY_DIR}/nuto")

    # create target doxgen
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

else()
    # create target doxgen
    add_custom_target(doc echo "WARNING: Doxygen not found -- source code documentation is not possible!")
endif()
