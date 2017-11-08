add_unit_test(Average)
add_unit_test(LinearInterpolation math/Interpolation.cpp)
add_unit_test(CubicSplineInterpolation math/Interpolation.cpp)
add_unit_test(EigenCompanion)
add_unit_test(SparseMatrixCSRGeneral math/SparseMatrixCSR.cpp)
add_unit_test(SparseMatrixCSRVector2General
        math/SparseMatrixCSRGeneral.cpp
        math/SparseMatrixCSR.cpp
    )
add_unit_test(Legendre)
add_unit_test(NaturalCoordinateMemoizer)
add_unit_test(NewtonRaphson
    math/SparseMatrixCSR.cpp
    math/SparseDirectSolverMUMPS.cpp
    base/Timer.cpp
    )
target_link_libraries(NewtonRaphson Mumps::Mumps)

add_unit_test(Gmres)
add_unit_test(SpatialContainer)
target_link_libraries(SpatialContainer Ann::Ann)

if(ARPACK_FOUND)
    add_unit_test(EigenSolverArpack)
    target_link_libraries(EigenSolverArpack Math Base)
endif()

if(ENABLE_MPI)
    add_unit_test(EpetraHelloWorld)
    target_link_libraries(EpetraHelloWorld trilinos_epetra)
    target_include_directories(EpetraHelloWorld SYSTEM PRIVATE ${Epetra_INCLUDE_DIRS})
endif()
