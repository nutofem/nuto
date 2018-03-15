configure_file(Version.cpp ${CMAKE_CURRENT_BINARY_DIR}/Version.cpp)

add_sources(${CMAKE_CURRENT_BINARY_DIR}/Version.cpp)
add_sources(Logger.cpp)
add_sources(serializeStream/SerializeStreamBase.cpp)
add_sources(serializeStream/SerializeStreamIn.cpp)
add_sources(serializeStream/SerializeStreamOut.cpp)
add_sources(Timer.cpp)
