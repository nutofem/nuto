
#include <iostream>

#include "nuto/base/Exception.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMKLPardiso.h"

int main()
{
    try
    {
        // create symmetric sparse matrix
        std::cout << "==================================" << std::endl;
        std::cout << "symmetric matrix, sparse CSR storage" << std::endl;
        NuTo::SparseMatrixCSRSymmetric<double> A_sy(5, 14);
        A_sy.AddValue(0, 1, 1);
        A_sy.AddValue(0, 2, 2);
        A_sy.AddValue(0, 3, 3);
        A_sy.AddValue(0, 4, 4);
        A_sy.AddValue(1, 1, 11);
        A_sy.AddValue(1, 2, 12);
        A_sy.AddValue(1, 3, 13);
        A_sy.AddValue(1, 4, 14);
        A_sy.AddValue(2, 2, 22);
        A_sy.AddValue(2, 3, 23);
        A_sy.AddValue(2, 4, 24);
        A_sy.AddValue(3, 3, 33);
        A_sy.AddValue(3, 4, 34);
        A_sy.AddValue(4, 4, 44);
        A_sy.SetOneBasedIndexing();
        A_sy.Info();
        std::cout << std::endl << "symmetric matrix, full storage" << std::endl;
        std::cout << A_sy.ConvertToFullMatrix() << std::endl;

        std::cout << "----------------------------------" << std::endl;
        std::cout << "* operator" << std::endl;
        NuTo::SparseMatrixCSRSymmetric<double> B_sy;
        B_sy = A_sy * 10.0;
        B_sy.Info();

        std::cout << "----------------------------------" << std::endl;
        std::cout << "*= operator" << std::endl;
        A_sy *= 10.0;
        std::cout << A_sy.ConvertToFullMatrix() << std::endl;


        std::cout << "==================================" << std::endl;
        std::cout << "general matrix, sparse CSR storage" << std::endl;
        NuTo::SparseMatrixCSRGeneral<double> A_ge(5, 5, 14);
        A_ge.AddValue(0, 1, 1);
        A_ge.AddValue(0, 2, 2);
        A_ge.AddValue(0, 3, 3);
        A_ge.AddValue(0, 4, 4);
        A_ge.AddValue(1, 1, 11);
        A_ge.AddValue(1, 2, 12);
        A_ge.AddValue(1, 3, 13);
        A_ge.AddValue(1, 4, 14);
        A_ge.AddValue(2, 2, 22);
        A_ge.AddValue(2, 3, 23);
        A_ge.AddValue(2, 4, 24);
        A_ge.AddValue(3, 3, 33);
        A_ge.AddValue(3, 4, 34);
        A_ge.AddValue(4, 4, 44);
        A_ge.SetOneBasedIndexing();
        A_ge.Info();
        std::cout << std::endl << "general matrix, full storage" << std::endl;
        std::cout << A_ge.ConvertToFullMatrix() << std::endl;

        std::cout << "----------------------------------" << std::endl;
        std::cout << "* operator" << std::endl;
        NuTo::SparseMatrixCSRGeneral<double> B_ge;
        B_ge = A_ge * 10.0;
        B_ge.Info();

        std::cout << "----------------------------------" << std::endl;
        std::cout << "*= operator" << std::endl;
        A_ge *= 10.0;
        A_ge.Info();

        std::cout << "==================================" << std::endl;
        std::cout << "multiply operators with full matrix" << std::endl;
        std::cout << "----------------------------------" << std::endl;
        std::cout << "* operator" << std::endl;
        std::cout << (A_sy * 0.1).ConvertToFullMatrix() << std::endl;
    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    return 0;
}
