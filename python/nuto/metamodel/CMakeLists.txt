configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_swig_module(ModuleMetamodel ModuleMetamodel.i Metamodel)
