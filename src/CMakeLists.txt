include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# add library
add_subdirectory(nuto)
install(DIRECTORY nuto DESTINATION ${NUTO_INCLUDE_INSTALL_PATH} FILES_MATCHING PATTERN "*.h")

# add GUI
if(ENABLE_GUI)
    add_subdirectory(gui)
endif()
