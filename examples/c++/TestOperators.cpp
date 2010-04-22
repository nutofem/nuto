// $Id$

#include <iostream>

#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
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
        NuTo::SparseMatrixCSRSymmetric<double> A_sy(5,14);
        A_sy.AddEntry(0,1,1);
        A_sy.AddEntry(0,2,2);
        A_sy.AddEntry(0,3,3);
        A_sy.AddEntry(0,4,4);
        A_sy.AddEntry(1,1,11);
        A_sy.AddEntry(1,2,12);
        A_sy.AddEntry(1,3,13);
        A_sy.AddEntry(1,4,14);
        A_sy.AddEntry(2,2,22);
        A_sy.AddEntry(2,3,23);
        A_sy.AddEntry(2,4,24);
        A_sy.AddEntry(3,3,33);
        A_sy.AddEntry(3,4,34);
        A_sy.AddEntry(4,4,44);
        A_sy.SetOneBasedIndexing();
        A_sy.Info();
        std::cout << std::endl << "symmetric matrix, full storage" << std::endl;
        NuTo::FullMatrix<double> A_sy_full(A_sy);
        A_sy_full.Info();

        std::cout << "----------------------------------" << std::endl;
        std::cout << "* operator" << std::endl;
        NuTo::SparseMatrixCSRSymmetric<double> B_sy;
        B_sy = A_sy * 10.0;
        B_sy.Info();
        NuTo::FullMatrix<double> B_sy_full(B_sy);
        //~ B_sy_full.Info();
        
        std::cout << "----------------------------------" << std::endl;
        std::cout << "*= operator" << std::endl;
		A_sy *= 10.0;
        A_sy.Info();
		A_sy_full=A_sy;
        //~ A_sy_full.Info();


        std::cout << "==================================" << std::endl;
        std::cout << "general matrix, sparse CSR storage" << std::endl;
        NuTo::SparseMatrixCSRGeneral<double> A_ge(5,5,14);
        A_ge.AddEntry(0,1,1);
        A_ge.AddEntry(0,2,2);
        A_ge.AddEntry(0,3,3);
        A_ge.AddEntry(0,4,4);
        A_ge.AddEntry(1,1,11);
        A_ge.AddEntry(1,2,12);
        A_ge.AddEntry(1,3,13);
        A_ge.AddEntry(1,4,14);
        A_ge.AddEntry(2,2,22);
        A_ge.AddEntry(2,3,23);
        A_ge.AddEntry(2,4,24);
        A_ge.AddEntry(3,3,33);
        A_ge.AddEntry(3,4,34);
        A_ge.AddEntry(4,4,44);
        A_ge.SetOneBasedIndexing();
        A_ge.Info();
        std::cout << std::endl << "general matrix, full storage" << std::endl;
        NuTo::FullMatrix<double> A_ge_full(A_ge);
        A_ge_full.Info();

        std::cout << "----------------------------------" << std::endl;
        std::cout << "* operator" << std::endl;
        NuTo::SparseMatrixCSRGeneral<double> B_ge;
        B_ge = A_ge * 10.0;
        B_ge.Info();
        NuTo::FullMatrix<double> B_ge_full(B_ge);
        //~ B_ge_full.Info();
        
        std::cout << "----------------------------------" << std::endl;
        std::cout << "*= operator" << std::endl;
		A_ge *= 10.0;
        A_ge.Info();
		A_ge_full=A_ge;
        //~ A_ge_full.Info();

        std::cout << "==================================" << std::endl;
        std::cout << "multiply operators with full matrix" << std::endl;
        std::cout << "----------------------------------" << std::endl;
        std::cout << "* operator" << std::endl;
        B_sy_full = A_sy * 0.1;
        B_sy_full.Info();
        
        std::cout << "----------------------------------" << std::endl;
        std::cout << "*= operator" << std::endl;
		B_sy_full *= 10.0;
        B_sy_full.Info();

    }
    catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    return 0;
}
