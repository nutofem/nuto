# TODO: these tests have been commented out at some point in time, we should
# try to resurrect some of them:
# - MisesPlasticity
# - Plane2D_patchtest
# - Temperature1D

if(ENABLE_MECHANICS)
    add_test(examples::py::Brick8N ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Brick8N.py)
    add_test(examples::py::Brick8NCoupling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Brick8NCoupling.py)
    add_test(examples::py::SparseDirectSolverMUMPS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SparseDirectSolverMUMPS.py)
    add_test(examples::py::Temperature1D ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Truss1D2N.py)
    add_test(examples::py::Truss1D2N ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Truss1D2N.py)
    add_test(examples::py::Truss1D2NCoupling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Truss1D2NCoupling.py)
    add_test(examples::py::VisualizeTruss1D2N ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/VisualizeTruss1D2N.py)

    set(pythontests
        examples::py::Brick8N
        examples::py::Brick8NCoupling
        examples::py::SparseDirectSolverMUMPS
        examples::py::Temperature1D
        examples::py::Truss1D2N
        examples::py::Truss1D2NCoupling
        examples::py::VisualizeTruss1D2N
    )

    if(PARDISO_FOUND)
        add_test(examples::py::SparseDirectSolverMKLPardiso ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SparseDirectSolverMKLPardiso.py)
        set(pythontests examples::py::SparseDirectSolverMKLPardiso)
    endif()

    if(MKLSolver_FOUND)
        add_test(EXAMPLES::PYTHON::SparseDirectSolverMKLDSS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SparseDirectSolverMKLDSS.py)
        set(MYTESTS EXAMPLES::PYTHON::SparseDirectSolverMKLDSS)
    endif()

    # set PYTHONPATH environment for all examples
    set_property(TEST ${pythontests} PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")
endif()

if(ENABLE_METAMODEL)
    add_test(examples::py::NeuralNetwork ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/NeuralNetwork.py)
    set_property(TEST examples::py::NeuralNetwork PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")
endif()

