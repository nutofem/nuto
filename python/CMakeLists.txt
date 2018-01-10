pybind11_add_module(metamodel NeuralNetwork.cpp ${CMAKE_SOURCE_DIR}/src/metamodel/NeuralNetwork.cpp)
target_include_directories(metamodel PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(metamodel PRIVATE Metamodel)
