# -*- coding: utf-8 -*-
import nuto
import numpy as np

def symmetric_solve():
    A = nuto.DoubleSparseMatrixCSRSymmetric(5,9)
    A.AddValue(0,0,9)
    A.AddValue(0,1,1.5)
    A.AddValue(0,2,6)
    A.AddValue(0,3,0.75)
    A.AddValue(0,4,3)
    A.AddValue(1,1,0.5)
    A.AddValue(2,2,12)
    A.AddValue(3,3,0.625)
    A.AddValue(4,4,16)
    A.SetOneBasedIndexing()
    print("symmetric matrix, sparse CSR storage")
    A.Info()
    print("\nsymmetric matrix, full storage")
    print(A.ConvertToFullMatrix())

    # right hand side vector
    rhs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    print("\nright hand side vector")
    print(rhs)

    # solver
    solver = nuto.SparseDirectSolverMUMPS()

    print("\nsolving the symmetric problem")
    sol = np.zeros(5)
    solver.Solve(A,rhs,sol)
    print("\nsolution of the symmetric problem")
    print(sol)

def nonsymmetric_solve():
    A = nuto.DoubleSparseMatrixCSRGeneral(5,5,13)
    A.AddValue(0,0,9)
    A.AddValue(0,1,1.5)
    A.AddValue(0,2,6)
    A.AddValue(0,3,0.75)
    A.AddValue(0,4,3)
    A.AddValue(1,0,1.5)
    A.AddValue(1,1,0.5)
    A.AddValue(2,0,6)
    A.AddValue(2,2,12)
    A.AddValue(3,0,0.75)
    A.AddValue(3,3,0.625)
    A.AddValue(4,0,3)
    A.AddValue(4,4,16)
    A.SetOneBasedIndexing()
    print("\nnonsymmetric matrix, sparse CSR storage")
    A.Info()
    print("\nnonsymmetric matrix, full storage")
    print(A.ConvertToFullMatrix())

    # right hand side vector
    rhs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    print("\nright hand side vector")
    print(rhs)

    # solver
    solver = nuto.SparseDirectSolverMUMPS()

    # solve nonsymmetric problem
    print("\nsolving the nonsymmetric problem")
    sol = np.zeros(5)
    solver.Solve(A,rhs,sol)
    print("\nsolution of the nonsymmetric problem")
    print(sol)

symmetric_solve()
nonsymmetric_solve()
