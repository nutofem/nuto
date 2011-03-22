add_definitions(-DRESULTVIEWER_EXPORTED=NUTOGUI_EXPORT)

file(GLOB ResultViewer_art art/*.png)
bundle_artprovider_art(ResultViewer_art_sources "image/png" ${ResultViewer_art})

file(GLOB ResultViewer_sources *.h *.cpp View3D/*.h  View3D/*.cpp)
add_library(ResultViewer STATIC ${ResultViewer_sources} ${ResultViewer_art_sources})
target_link_libraries(ResultViewer vtkwx vtkWidgets vtkCharts)
if(STATIC_LIBRARY_FLAGS)
  set_target_properties(ResultViewer PROPERTIES COMPILE_FLAGS ${STATIC_LIBRARY_FLAGS})
endif()
