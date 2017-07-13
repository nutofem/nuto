#include "BoostUnitTest.h"
#include "math/Gmres.h"
#include <iostream>


using MatrixType = Eigen::MatrixXd;
using VectorType = Eigen::VectorXd;
using SparseMatrixType = Eigen::SparseMatrix<double>;


BOOST_AUTO_TEST_CASE(SolveSystemIdentity)
{

    const int dim = 10;

    MatrixType A(dim, dim);
    A.setIdentity();
    A *= 1337;

    std::cout << "A \n" << A << "\n\n";

    VectorType x(dim);
    x.setZero();

    VectorType b(dim);
    b.setOnes();

    int maxNumRestarts = 100;
    int krylovDim = 1;
    double tolerance = 1.e-6;

    int numRestarts = NuTo::Gmres<MatrixType>(A, b, x, maxNumRestarts, tolerance, krylovDim);

    BOOST_CHECK((A * x - b).isMuchSmallerThan(1.e-4, 1.e-1));

    std::cout << "#restarts \t" << numRestarts << std::endl;
}


BOOST_AUTO_TEST_CASE(SolveSystemSpd)
{

    const int dim = 4;

    MatrixType A(dim, dim);
    A << 1, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2;

    std::cout << "A \n" << A << "\n\n";

    VectorType x(dim);
    x.setZero();

    VectorType b(dim);
    b.setOnes();


    int maxNumRestarts = 100;
    int krylovDim = 2;
    double tolerance = 1.e-6;

    int numRestarts = NuTo::Gmres<MatrixType>(A, b, x, maxNumRestarts, tolerance, krylovDim);

    BOOST_CHECK((A * x - b).isMuchSmallerThan(1.e-4, 1.e-1));

    std::cout << "#restarts \t" << numRestarts << std::endl;
}


BOOST_AUTO_TEST_CASE(SolveSystemSparse)
{

    const int dim = 4;

    SparseMatrixType A(dim, dim);
    // row 0
    A.insert(0, 0) = 1;
    A.insert(0, 1) = -1;

    // row 1
    A.insert(1, 0) = -1;
    A.insert(1, 1) = 2;
    A.insert(1, 2) = -1;

    // row 2
    A.insert(2, 1) = -1;
    A.insert(2, 2) = 2;
    A.insert(2, 3) = -1;

    // row 3
    A.insert(3, 2) = -1;
    A.insert(3, 3) = 2;

    std::cout << "A \n" << A << "\n\n";

    VectorType x(dim);
    x.setZero();

    VectorType b(dim);
    b.setOnes();

    int maxNumRestarts = 100;
    int krylovDim = 2;
    double tolerance = 1.e-6;

    int numRestarts = NuTo::Gmres<MatrixType>(A, b, x, maxNumRestarts, tolerance, krylovDim);
    BOOST_CHECK((A * x - b).isMuchSmallerThan(1.e-4, 1.e-1));
    std::cout << "#restarts \t" << numRestarts << "\t <- diagonal preconditioner" << std::endl;

    x.setZero();
    numRestarts = NuTo::Gmres<MatrixType, Eigen::IdentityPreconditioner>(A, b, x, maxNumRestarts, tolerance, krylovDim);
    BOOST_CHECK((A * x - b).isMuchSmallerThan(1.e-4, 1.e-1));
    std::cout << "#restarts \t" << numRestarts << "\t <- identity preconditioner" << std::endl;
}
