set(NuToBase_SOURCES
    Logger.cpp
    CallbackInterface.cpp
    serializeStream/SerializeStreamBase.cpp
    serializeStream/SerializeStreamIn.cpp
    serializeStream/SerializeStreamOut.cpp
    )

# build static library
add_library(NuToBase ${NuToBase_SOURCES})
target_link_libraries(NuToBase ${Boost_LIBRARIES})
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set_target_properties(NuToBase PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

# install library
install(TARGETS NuToBase
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )
