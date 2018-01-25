add_library(Fakeit INTERFACE)
target_include_directories(Fakeit SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(Json INTERFACE)
target_include_directories(Json SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
