configure_file(Version.cpp ${CMAKE_CURRENT_BINARY_DIR}/Version.cpp)

set(BaseSources
    Logger.cpp
    CallbackInterface.cpp
    serializeStream/SerializeStreamBase.cpp
    serializeStream/SerializeStreamIn.cpp
    serializeStream/SerializeStreamOut.cpp
    Timer.cpp
    )

sources_to_objects("${BaseSources}" BaseObjects)

add_library(Base ${BaseObjects} ${CMAKE_CURRENT_BINARY_DIR}/Version.cpp)
set_target_properties(Base PROPERTIES
    OUTPUT_NAME NuToBase
    )

target_link_libraries(Base ${Boost_LIBRARIES})
target_include_directories(Base PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/nuto>
    )

install(TARGETS Base EXPORT NuToTargets
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include/nuto
    )
