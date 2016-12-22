configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_install_python_file(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/metamodel)
nuto_swig_module(ModulMetaModel ModulMetaModel.i nuto/metamodel NuToMetamodel)
