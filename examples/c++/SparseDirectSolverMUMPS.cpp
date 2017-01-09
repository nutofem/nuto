#include <iostream>

#include "math/MathException.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseDirectSolverMUMPS.h"
int main()
{
    try
    {
        // create symmetric sparse matrix
        NuTo::SparseMatrixCSRSymmetric<double> A_sy(5,9);
        A_sy.AddValue(0,0,9);
        A_sy.AddValue(0,1,1.5);
        A_sy.AddValue(0,2,6);
        A_sy.AddValue(0,3,0.75);
        A_sy.AddValue(0,4,3);
        A_sy.AddValue(1,1,0.5);
        A_sy.AddValue(2,2,12);
        A_sy.AddValue(3,3,0.625);
        A_sy.AddValue(4,4,16);
        A_sy.SetOneBasedIndexing();
        std::cout << "symmetric matrix, sparse CSR storage" << std::endl;
        A_sy.Info();
        std::cout << std::endl << "symmetric matrix, full storage" << std::endl;
        std::cout << A_sy.ConvertToFullMatrix() << std::endl;

        // nonsymmetric coefficient matrix
        NuTo::SparseMatrixCSRGeneral<double> A_nosy(5,5,13);
        A_nosy.AddValue(0,0,9);
        A_nosy.AddValue(0,1,1.5);
        A_nosy.AddValue(0,2,6);
        A_nosy.AddValue(0,3,0.75);
        A_nosy.AddValue(0,4,3);
        A_nosy.AddValue(1,0,1.5);
        A_nosy.AddValue(1,1,0.5);
        A_nosy.AddValue(2,0,6);
        A_nosy.AddValue(2,2,12);
        A_nosy.AddValue(3,0,0.75);
        A_nosy.AddValue(3,3,0.625);
        A_nosy.AddValue(4,0,3);
        A_nosy.AddValue(4,1,1);
        A_nosy.AddValue(4,4,16);
        A_nosy.SetOneBasedIndexing();
        std::cout << std::endl << "nonsymmetric matrix, sparse CSR storage" << std::endl;
        A_nosy.Info();
        std::cout << std::endl << "nonsymmetric matrix, full storage" << std::endl;
        std::cout << A_nosy.ConvertToFullMatrix() << std::endl;

        // create right hand side vector
        Eigen::VectorXd rhs(5);
        rhs[0] = 1;
        rhs[1] = 2;
        rhs[2] = 3;
        rhs[3] = 4;
        rhs[4] = 5;
        std::cout << std::endl << "right hand side vector \n" << rhs << std::endl;

        // create solver
        NuTo::SparseDirectSolverMUMPS mumps;
        mumps.SetVerboseLevel(3);

        // solve symmetric problem
        std::cout << std::endl << "solving the symmetric problem" << std::endl;
        Eigen::VectorXd sol_sy(5);
        mumps.Solve(A_sy,rhs,sol_sy);
        std::cout << std::endl << "solution of the symmetric problem \n" << sol_sy << std::endl;

        // solve nonsymmetric problem
        std::cout << std::endl << "solving the nonsymmetric problem" << std::endl;
        Eigen::VectorXd sol_nosy(5);
        mumps.Solve(A_nosy,rhs,sol_nosy);
        std::cout << std::endl << "solution of the nonsymmetric problem \n" << sol_nosy << std::endl;

        // solve for Schur complement
        std::cout << std::endl << "solving the Schur complement of A with respect to indices 0 and 4" << std::endl;
        Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> schur_Indices(2,1);
        //attention - zero based indexing for the indices
        schur_Indices(0,0) = 0;
        schur_Indices(1,0) = 4;
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> schur_complement(2,2);
        mumps.SchurComplement(A_nosy,schur_Indices,schur_complement);
        std::cout << schur_complement << std::endl;
        //correct solution is [0.6 3 ]
        //                    [0   16]

    }
    catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    return 0;
}
