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

# build python modules
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_install_python_file(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/base)
nuto_swig_module(ModulNuToBase ModulNuToBase.i nuto/base NuToBase)
