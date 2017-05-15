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
add_library(NuToBase ${NuToBaseObjects})
target_link_libraries(NuToBase ${Boost_LIBRARIES})

# install library
install(TARGETS NuToBase
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )
