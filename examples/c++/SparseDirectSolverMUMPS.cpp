#include <iostream>

#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
int main()
{
    try
    {
        // create symmetric sparse matrix
        NuTo::SparseMatrixCSRSymmetric<double> A_sy(5,9);
        A_sy.AddEntry(0,0,9);
        A_sy.AddEntry(0,1,1.5);
        A_sy.AddEntry(0,2,6);
        A_sy.AddEntry(0,3,0.75);
        A_sy.AddEntry(0,4,3);
        A_sy.AddEntry(1,1,0.5);
        A_sy.AddEntry(2,2,12);
        A_sy.AddEntry(3,3,0.625);
        A_sy.AddEntry(4,4,16);
        A_sy.SetOneBasedIndexing();
        std::cout << "symmetric matrix, sparse CSR storage" << std::endl;
        A_sy.Info();
        std::cout << std::endl << "symmetric matrix, full storage" << std::endl;
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> A_sy_full(A_sy);
        A_sy_full.Info(12,3);

        // nonsymmetric coefficient matrix
        NuTo::SparseMatrixCSRGeneral<double> A_nosy(5,5,13);
        A_nosy.AddEntry(0,0,9);
        A_nosy.AddEntry(0,1,1.5);
        A_nosy.AddEntry(0,2,6);
        A_nosy.AddEntry(0,3,0.75);
        A_nosy.AddEntry(0,4,3);
        A_nosy.AddEntry(1,0,1.5);
        A_nosy.AddEntry(1,1,0.5);
        A_nosy.AddEntry(2,0,6);
        A_nosy.AddEntry(2,2,12);
        A_nosy.AddEntry(3,0,0.75);
        A_nosy.AddEntry(3,3,0.625);
        A_nosy.AddEntry(4,0,3);
        A_nosy.AddEntry(4,1,1);
        A_nosy.AddEntry(4,4,16);
        A_nosy.SetOneBasedIndexing();
        std::cout << std::endl << "nonsymmetric matrix, sparse CSR storage" << std::endl;
        A_nosy.Info();
        std::cout << std::endl << "nonsymmetric matrix, full storage" << std::endl;
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> A_nosy_full(A_nosy);
        A_nosy_full.Info(12,3);

        // create right hand side vector
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rhs(5,1);
        rhs.SetValue(0,0,1);
        rhs.SetValue(1,0,2);
        rhs.SetValue(2,0,3);
        rhs.SetValue(3,0,4);
        rhs.SetValue(4,0,5);
        std::cout << std::endl << "right hand side vector" << std::endl;
        rhs.Info(12,3);

        // create solver
        NuTo::SparseDirectSolverMUMPS mumps;
        mumps.SetVerboseLevel(3);

        // solve symmetric problem
        std::cout << std::endl << "solving the symmetric problem" << std::endl;
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> sol_sy(5,1);
        mumps.Solve(A_sy,rhs,sol_sy);
        std::cout << std::endl << "solution of the symmetric problem" << std::endl;
        sol_sy.Info(12,3);

        // solve nonsymmetric problem
        std::cout << std::endl << "solving the nonsymmetric problem" << std::endl;
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> sol_nosy(5,1);
        mumps.Solve(A_nosy,rhs,sol_nosy);
        std::cout << std::endl << "solution of the nonsymmetric problem" << std::endl;
        sol_nosy.Info(12,3);

        // solve for Schur complement
        std::cout << std::endl << "solving the Schur complement of A with respect to indices 0 and 4" << std::endl;
        NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> schur_Indices(2,1);
        //attention - zero based indexing for the indices
        schur_Indices(0,0) = 0;
        schur_Indices(1,0) = 4;
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> schur_complement(2,2);
        mumps.SchurComplement(A_nosy,schur_Indices,schur_complement);
        schur_complement.Info(12,3); //correct solution is [0.6 3 ]
                                     //                    [0   16]

    }
    catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    return 0;
}
