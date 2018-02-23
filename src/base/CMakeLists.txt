configure_file(Version.cpp ${CMAKE_CURRENT_BINARY_DIR}/Version.cpp)

set(BaseSources
    Logger.cpp
    serializeStream/SerializeStreamBase.cpp
    serializeStream/SerializeStreamIn.cpp
    serializeStream/SerializeStreamOut.cpp
    Timer.cpp
    )

create_nuto_module(Base "${BaseSources}"
    ${CMAKE_CURRENT_BINARY_DIR}/Version.cpp)
target_link_libraries(Base PRIVATE Eigen3::Eigen)
