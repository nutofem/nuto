# -*- coding: utf-8 -*-
import nuto

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
    print "symmetric matrix, sparse CSR storage"
    A.Info()
    A_full = nuto.DoubleFullMatrix(A)
    print "\nsymmetric matrix, full storage"
    A_full.Info(3)

    # right hand side vector
    rhs = nuto.DoubleFullVector([1,2,3,4,5])
    print "\nright hand side vector"
    rhs.Info()

    # solver
    solver = nuto.SparseDirectSolverMUMPS()
    solver.SetVerboseLevel(3)

    print "\nsolving the symmetric problem"
    sol = nuto.DoubleFullVector(5)
    solver.Solve(A,rhs,sol)
    print "\nsolution of the symmetric problem"
    sol.Info()

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
    print "\nnonsymmetric matrix, sparse CSR storage"
    A.Info()
    print "\nnonsymmetric matrix, full storage"
    A_full = nuto.DoubleFullMatrix(A)
    A_full.Info(3)

    # right hand side vector
    rhs = nuto.DoubleFullVector([1,2,3,4,5])
    print "\nright hand side vector"
    rhs.Info()

    # solver
    solver = nuto.SparseDirectSolverMUMPS()
    solver.SetVerboseLevel(3)

    # solve nonsymmetric problem
    print "\nsolving the nonsymmetric problem"
    sol = nuto.DoubleFullVector(5)
    solver.Solve(A,rhs,sol)
    print "\nsolution of the nonsymmetric problem"
    sol.Info()

symmetric_solve()
nonsymmetric_solve()
