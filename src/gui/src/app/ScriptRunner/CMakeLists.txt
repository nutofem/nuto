include(FindPythonLibs)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
find_package_message(SWIG "Found SWIG: ${SWIG_EXECUTABLE}" ${SWIG_EXECUTABLE})

# Generate SWIG runtime
set(SWIG_RUNTIME_DIR "${CMAKE_CURRENT_BINARY_DIR}/SwigRuntime")
set(SWIG_RUNTIME_HEADER "${SWIG_RUNTIME_DIR}/SwigRuntime.h")
file(MAKE_DIRECTORY "${SWIG_RUNTIME_DIR}")
add_custom_command(OUTPUT ${SWIG_RUNTIME_HEADER}
		   COMMAND "${SWIG_EXECUTABLE}" -python -c++ -external-runtime "${SWIG_RUNTIME_HEADER}")
set_source_files_properties(${SWIG_RUNTIME_HEADER} PROPERTIES GENERATED TRUE)
include_directories(${SWIG_RUNTIME_DIR})

file(GLOB ScriptRunner_sources *.h *.cpp)
add_library(ScriptRunner STATIC ${ScriptRunner_sources} ${SWIG_RUNTIME_HEADER})
target_link_libraries(ScriptRunner NuToMechanics NuToMath NuToVisualize NuToBase)
target_link_libraries(ScriptRunner ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
