configure_file(Version.cpp ${CMAKE_CURRENT_BINARY_DIR}/Version.cpp)

set(BaseSources
    Logger.cpp
    CallbackInterface.cpp
    serializeStream/SerializeStreamBase.cpp
    serializeStream/SerializeStreamIn.cpp
    serializeStream/SerializeStreamOut.cpp
    tools/GNUPlot.cpp
    tools/PlotData.cpp
    tools/PlotEnum.cpp
    Timer.cpp
    )

create_nuto_module(Base "${BaseSources}" ${CMAKE_CURRENT_BINARY_DIR}/Version.cpp)
target_link_libraries(Base ${Boost_LIBRARIES})
