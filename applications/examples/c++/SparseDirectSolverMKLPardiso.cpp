#include <iostream>

#include "math/MathException.h"
#include "math/FullMatrix.h"
#include "math/FullVector.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseDirectSolverMKLPardiso.h"

int main()
{
    try
    {
        // create symmetric sparse matrix
        NuTo::SparseMatrixCSRSymmetric<double> A_sy(5, 9);
        A_sy.AddValue(0, 0, 9);
        A_sy.AddValue(0, 1, 1.5);
        A_sy.AddValue(0, 2, 6);
        A_sy.AddValue(0, 3, 0.75);
        A_sy.AddValue(0, 4, 3);
        A_sy.AddValue(1, 1, 0.5);
        A_sy.AddValue(2, 2, 12);
        A_sy.AddValue(3, 3, 0.625);
        A_sy.AddValue(4, 4, 16);
        A_sy.SetOneBasedIndexing();
        std::cout << "symmetric matrix, sparse CSR storage" << std::endl;
        A_sy.Info();
        std::cout << std::endl << "symmetric matrix, full storage" << std::endl;
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> A_sy_full(A_sy);
        A_sy_full.Info();

        // nonsymmetric coefficient matrix
        NuTo::SparseMatrixCSRGeneral<double> A_nosy(5, 5, 13);
        A_nosy.AddValue(0, 0, 9);
        A_nosy.AddValue(0, 1, 1.5);
        A_nosy.AddValue(0, 2, 6);
        A_nosy.AddValue(0, 3, 0.75);
        A_nosy.AddValue(0, 4, 3);
        A_nosy.AddValue(1, 0, 1.5);
        A_nosy.AddValue(1, 1, 0.5);
        A_nosy.AddValue(2, 0, 6);
        A_nosy.AddValue(2, 2, 12);
        A_nosy.AddValue(3, 0, 0.75);
        A_nosy.AddValue(3, 3, 0.625);
        A_nosy.AddValue(4, 0, 3);
        A_nosy.AddValue(4, 4, 16);
        A_nosy.SetOneBasedIndexing();
        std::cout << std::endl << "nonsymmetric matrix, sparse CSR storage" << std::endl;
        A_nosy.Info();
        std::cout << std::endl << "nonsymmetric matrix, full storage" << std::endl;
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> A_nosy_full(A_nosy);
        A_nosy_full.Info(3);

        // create right hand side vector
        NuTo::FullVector<double, Eigen::Dynamic> rhs(5, 1);
        rhs.SetValue(0, 0, 1);
        rhs.SetValue(1, 0, 2);
        rhs.SetValue(2, 0, 3);
        rhs.SetValue(3, 0, 4);
        rhs.SetValue(4, 0, 5);
        std::cout << std::endl << "right hand side vector" << std::endl;
        rhs.Info();

        // create solver
        NuTo::SparseDirectSolverMKLPardiso pardiso;
        pardiso.SetVerboseLevel(3);

        // solve symmetric problem
        std::cout << std::endl << "solving the symmetric problem" << std::endl;
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> sol_sy(5, 1);
        pardiso.Solve(A_sy, rhs, sol_sy);
        std::cout << std::endl << "solution of the symmetric problem" << std::endl;
        sol_sy.Info();

        // solve nonsymmetric problem
        std::cout << std::endl << "solving the nonsymmetric problem" << std::endl;
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> sol_nosy(5, 1);
        pardiso.Solve(A_nosy, rhs, sol_nosy);
        std::cout << std::endl << "solution of the nonsymmetric problem" << std::endl;
        sol_nosy.Info();
    }
    catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    return 0;
}
