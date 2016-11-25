

#define BOOST_TEST_MODULE ArpackTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "nuto/base/Timer.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/math/EigenSolverArpack.h"
#include "nuto/math/EigenSolverArpackEnum.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"


// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}



struct ArpackTestFixture
{
    ArpackTestFixture()
    {
        BOOST_TEST_MESSAGE( "setup fixture" );
        auto A_Full = GetA();
        auto M_Full = GetM();

        A_Full.EigenVectorsSymmetric(mEigenValuesStandard,mEigenVectorsStandard);

        A_Full.GeneralizedEigenVectorsSymmetric(M_Full,mEigenValuesGeneral,mEigenVectorsGeneral);

        mEigenValuesStandardRef = mEigenValuesStandard.GetBlock(8 - mNumEigenValuesCompute, 0, mNumEigenValuesCompute, 1);
        mEigenValuesGeneralRef  = mEigenValuesGeneral.GetBlock (8 - mNumEigenValuesCompute, 0, mNumEigenValuesCompute, 1);

        mEigenSolver.SetShowTime(false);

    }

    static NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> GetA()
    {
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> a(8,8);
        a(0,0) = 8.;
        a(1,1) = 2.;
        a(2,2) = 4.;
        a(3,3) = 7.;
        a(4,4) = 10.;
        a(5,5) = 6.;
        a(6,6) = 5.;
        a(7,7) = 1.;
        a(0,1) = 1.;
        a(1,0) = 1.;
        a(0,3) = -1.;
        a(3,0) = -1.;
        a(2,4) = 2.;
        a(4,2) = 2.;
        a(6,1) = -3.;
        a(1,6) = -3.;
        return a;
    }

    static NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> GetM()
    {
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> m(8,8);
        m(0,0) = 1.;
        m(1,1) = 1.;
        m(2,2) = 1.;
        m(3,3) = 1.;
        m(4,4) = 1.;
        m(5,5) = 1.;
        m(6,6) = 1.;
        m(7,7) = 1.;

        m(0,1) = 0.5;
        m(1,0) = 0.5;
        return m;
    }

    //values for the standard eigenvalue problem
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mEigenVectorsStandard;
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mEigenValuesStandard;

    //values for the generalized eigenvalue problem
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mEigenVectorsGeneral;
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mEigenValuesGeneral;

    int mNumEigenValuesCompute = 3;

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mEigenValuesStandardRef;

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mEigenValuesGeneralRef;

    NuTo::EigenSolverArpack mEigenSolver;

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mEigenVectors;
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mEigenValues;

};

BOOST_FIXTURE_TEST_CASE(Info, ArpackTestFixture)
{
    std::cout << "FullMatrix\n" << GetA() << std::endl;
    std::cout << "eigenvalues of standard eigenvalue problem.\n" << mEigenValuesStandard.Trans() << std::endl;
    std::cout << "eigenvalues of generalized eigenvalue problem.\n" << mEigenValuesGeneral.Trans() << std::endl;
}


BOOST_FIXTURE_TEST_CASE(DSDRV1, ArpackTestFixture)
{

    NuTo::SparseMatrixCSRVector2Symmetric<double> A_symmetric(GetA());

    mEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV1);
    //largest magnitude
    mEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::eWhich::LA);
    mEigenSolver.Solve(A_symmetric, nullptr, mNumEigenValuesCompute, mEigenValues, mEigenVectors);

    std::cout << "largest eigenvalues of DSDRV1\n" << mEigenValues.Trans() << std::endl;
    BOOST_CHECK_SMALL((mEigenValues-mEigenValuesStandardRef).norm() , 1e-5);
}


BOOST_FIXTURE_TEST_CASE(DSDRV2, ArpackTestFixture)
{

    NuTo::SparseMatrixCSRVector2Symmetric<double> A_symmetric(GetA());

    mEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV2);
    mEigenSolver.SetShift(-10., 0);
    // for DSDRV2, the which flag refers to the modified problem (A-shift*I)*x = lambda * x
    //the eigenvalues of the original problem are x_orig = shift + 1/x_mod
    //-->for a shift of zero, the largest eigenvalue is computed with SA (smallest amplitude)
    //if you want to calculate all eigenvalues closest to an given value, set Shift to this value and which to 'LM' (fabs - both directions)
    mEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::eWhich::SA);
    mEigenSolver.Solve(A_symmetric, nullptr, mNumEigenValuesCompute, mEigenValues, mEigenVectors);

    std::cout << "largest eigenvalues of DSDRV2\n" << mEigenValues.Trans() << std::endl;
    BOOST_CHECK_SMALL((mEigenValues-mEigenValuesStandardRef).norm() , 1e-5);
}

BOOST_FIXTURE_TEST_CASE(DSDRV3, ArpackTestFixture)
{
    NuTo::SparseMatrixCSRVector2Symmetric<double> A_symmetric(GetA());
    NuTo::SparseMatrixCSRVector2Symmetric<double> M_symmetric(GetM());

    mEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV3);
    mEigenSolver.SetShift(0.0); //not used
//    //for DSDRV3, the which flag refers to the modified problem M^(-1)A*x = lambda * x
    mEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::eWhich::LA);
    mEigenSolver.Solve(A_symmetric, &M_symmetric, mNumEigenValuesCompute, mEigenValues, mEigenVectors);

    std::cout << "largest eigenvalues of DSDRV3\n" << mEigenValues.Trans() << std::endl;
    BOOST_CHECK_SMALL((mEigenValues-mEigenValuesGeneralRef).norm() , 1e-5);
}

BOOST_FIXTURE_TEST_CASE(DSDRV4, ArpackTestFixture)
{
    NuTo::SparseMatrixCSRVector2Symmetric<double> A_symmetric(GetA());
    NuTo::SparseMatrixCSRVector2Symmetric<double> M_symmetric(GetM());

    mEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV4);
    mEigenSolver.SetShift(-10.0);
    //for DSDRV4, the which flag refers to the modified problem (A-shift*M)^(-1) *M *x = lambda * x
    //the eigenvalues of the original problem are x_orig = shift + 1/x_mod
    //-->for a shift of zero, the largest eigenvalue is computed with SA (smallest amplitude)
    //if you want to calculate all eigenvalues closest to an given value, set Shift to this value and which to 'LM' (fabs - both directions)
    mEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::eWhich::SA);
    mEigenSolver.Solve(A_symmetric, &M_symmetric, mNumEigenValuesCompute, mEigenValues, mEigenVectors);

    std::cout << "largest eigenvalues of DSDRV4\n" << mEigenValues.Trans() << std::endl;
    BOOST_CHECK_SMALL((mEigenValues-mEigenValuesGeneralRef).norm() , 1e-5);
}

BOOST_FIXTURE_TEST_CASE(DSDRV5, ArpackTestFixture)
{
    NuTo::SparseMatrixCSRVector2Symmetric<double> A_symmetric(GetA());
    NuTo::SparseMatrixCSRVector2Symmetric<double> M_symmetric(GetM());

    mEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV5);
	mEigenSolver.SetShift(-10.0);
	mEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::eWhich::LA);
	mEigenSolver.Solve(A_symmetric, &M_symmetric, mNumEigenValuesCompute, mEigenValues, mEigenVectors);

    std::cout << "largest eigenvalues of DSDRV5\n" << mEigenValues.Trans() << std::endl;
    BOOST_CHECK_SMALL((mEigenValues-mEigenValuesGeneralRef).norm() , 1e-5);
}

BOOST_FIXTURE_TEST_CASE(DNDRV1, ArpackTestFixture)
{

    NuTo::SparseMatrixCSRVector2General<double> A_general(GetA());
    //set driver for unsymmetric solution
    mEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::eDriver::DNDRV1);
    mEigenSolver.SetShift(0.0); //not used
    //largest real
    mEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::eWhich::LR);
    mEigenSolver.Solve(A_general, nullptr, mNumEigenValuesCompute, mEigenValues, mEigenVectors);
    //sort the eigenvalues (I don't know why this is not done in ARPACK)
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValuesSorted = mEigenValues.SortRow(0).col(0);

    std::cout << "largest eigenvalues of DNDRV1\n" << eigenValuesSorted.Trans() << std::endl;

    BOOST_CHECK_SMALL((eigenValuesSorted-mEigenValuesStandardRef).norm() , 1e-5);
}

BOOST_FIXTURE_TEST_CASE(LargestEVSymmetric, ArpackTestFixture)
{
    NuTo::SparseMatrixCSRVector2Symmetric<double> A_symmetric(GetA());
    auto result = mEigenSolver.GetLargest(A_symmetric);

    BOOST_CHECK_CLOSE(result.first, mEigenValuesStandard.maxCoeff(), 1e-5);
}
BOOST_FIXTURE_TEST_CASE(SmallestEVSymmetric, ArpackTestFixture)
{
    NuTo::SparseMatrixCSRVector2Symmetric<double> A_symmetric(GetA());
    auto result = mEigenSolver.GetSmallest(A_symmetric);

    BOOST_CHECK_CLOSE(result.first, mEigenValuesStandard.minCoeff(), 1e-5);
}



BOOST_FIXTURE_TEST_CASE(LargestEVGeneral, ArpackTestFixture)
{
    NuTo::SparseMatrixCSRVector2General<double> A_general(GetA());
    auto result = mEigenSolver.GetLargest(A_general);

    BOOST_CHECK_CLOSE(result.first, mEigenValuesStandard.maxCoeff(), 1e-5);
}

BOOST_FIXTURE_TEST_CASE(SmallestEVGeneral, ArpackTestFixture)
{
    NuTo::SparseMatrixCSRVector2General<double> A_general(GetA());
    auto result = mEigenSolver.GetSmallest(A_general);

    BOOST_CHECK_CLOSE(result.first, mEigenValuesStandard.minCoeff(), 1e-5);
}
