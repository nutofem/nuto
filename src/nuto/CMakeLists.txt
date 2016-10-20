configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_install_python_file(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto)

add_subdirectory(base)
add_subdirectory(math)
if(ENABLE_VISUALIZE)
    add_subdirectory(visualize)
endif()

if(ENABLE_MECHANICS)
    add_subdirectory(mechanics)
else()
    set(SEARCH_REGEX "from nuto.mechanics.ModulMechanics")
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/__init__.py FILE_CONTENT)
    string(REGEX REPLACE "${SEARCH_REGEX}" "#from nuto.mechanics.ModulMechanics"
        MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "${MODIFIED_FILE_CONTENT}")
endif()

if(ENABLE_GEOMETRYCONCRETE)
    add_subdirectory(geometryConcrete)
else()
    set(SEARCH_REGEX "from nuto.geometryConcrete.ModulGeometryConcrete")
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/__init__.py FILE_CONTENT)
    string(REGEX REPLACE "${SEARCH_REGEX}" "#from nuto.geometryConcrete.ModulGeometryConcrete"
        MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "${MODIFIED_FILE_CONTENT}")
endif()

if(ENABLE_METAMODEL)
    add_subdirectory(metamodel)
else()
    set(SEARCH_REGEX "from nuto.metamodel.ModulMetaModel")
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/__init__.py FILE_CONTENT)
    string(REGEX REPLACE "${SEARCH_REGEX}" "#from nuto.metamodel.ModulMetaModel"
        MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "${MODIFIED_FILE_CONTENT}")
endif()

if(ENABLE_OPTIMIZE)
    add_subdirectory(optimize)
else()
    set(SEARCH_REGEX "from nuto.optimize.ModulOptimizer")
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/__init__.py FILE_CONTENT)
    string(REGEX REPLACE "${SEARCH_REGEX}" "#from nuto.optimize.ModulOptimizer"
        MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "${MODIFIED_FILE_CONTENT}")
endif()
