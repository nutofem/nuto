if(ENABLE_MECHANICS)
    add_test(examples::py::Brick8NCoupling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Brick8NCoupling.py)
    add_test(examples::py::DamageBar ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/DamageBar.py --no-plot)
    add_test(examples::py::SparseDirectSolverMUMPS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SparseDirectSolverMUMPS.py)
    add_test(examples::py::Truss1D2N ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Truss1D2N.py)

    set(pythontests
        examples::py::Brick8NCoupling
        examples::py::DamageBar
        examples::py::SparseDirectSolverMUMPS
        examples::py::Truss1D2N
    )

    if(MKLSolver_FOUND)
        add_test(examples::py::SparseDirectSolverMKLDSS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SparseDirectSolverMKLDSS.py)
        set(pythontests "${pythontests};examples::py::SparseDirectSolverMKLDSS")
        
        add_test(examples::py::SparseDirectSolverMKLPardiso ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SparseDirectSolverMKLPardiso.py)
        set(pythontests "${pythontests};examples::py::SparseDirectSolverMKLPardiso")
    endif()

    # set PYTHONPATH environment for all examples
    set_property(TEST ${pythontests} PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
endif()

if(ENABLE_METAMODEL)
    add_test(examples::py::NeuralNetwork ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/NeuralNetwork.py)
    set_property(TEST examples::py::NeuralNetwork PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
endif()

