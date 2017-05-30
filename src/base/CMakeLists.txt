include(Version)
configure_file(Version.cpp ${CMAKE_CURRENT_BINARY_DIR}/Version.cpp)

set(NuToBaseSources
    Logger.cpp
    CallbackInterface.cpp
    serializeStream/SerializeStreamBase.cpp
    serializeStream/SerializeStreamIn.cpp
    serializeStream/SerializeStreamOut.cpp
    Timer.cpp
    )

sources_to_objects("${NuToBaseSources}" NuToBaseObjects)

# build static library
add_library(NuToBase ${NuToBaseObjects} ${CMAKE_CURRENT_BINARY_DIR}/Version.cpp)
target_link_libraries(NuToBase ${Boost_LIBRARIES})

# install library
install(TARGETS NuToBase
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )
