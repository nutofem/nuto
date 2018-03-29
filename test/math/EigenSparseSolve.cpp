//#define EIGEN_USE_MKL_ALL
#include "BoostUnitTest.h"
#include "nuto/math/EigenSparseSolve.h"
#include "nuto/base/Exception.h"
#include <boost/test/data/test_case.hpp>

namespace bdata = boost::unit_test::data;

using namespace NuTo;

struct LinearSystem
{
    LinearSystem()
        : b(3, -3, -2)
        , expected_x(1, -2, -2)
    {
        A.resize(3, 3);

        // not very sparse, but hey...
        A.insert(0, 0) = 1.0;
        A.insert(0, 1) = -1.0;
        A.insert(1, 0) = -1.0;
        A.insert(1, 1) = 2.0;
        A.insert(1, 2) = -1.0;
        A.insert(2, 1) = -1.0;
        A.insert(2, 2) = 2.0;
        A.makeCompressed();
    }

    Eigen::SparseMatrix<double> A;
    Eigen::Vector3d b;
    Eigen::Vector3d expected_x;
};

auto builtInSolverNames = {"EigenSparseLU",       "EigenSparseQR",          "EigenSimplicialLLT",
                           "EigenSimplicialLDLT", "EigenConjugateGradient", "EigenLeastSquaresConjugateGradient",
                           "EigenBiCGSTAB"};

BOOST_DATA_TEST_CASE(builtInSolvers, bdata::make(builtInSolverNames), solver)
{
    LinearSystem sys;
    auto x = EigenSparseSolve(sys.A, sys.b, solver);
    BoostUnitTest::CheckVector(x, sys.expected_x, 3);
}

auto suiteSparseSolverNames = {"SuiteSparseLU", "SuiteSparseSupernodalLLT"};

BOOST_DATA_TEST_CASE(suiteSparse, bdata::make(suiteSparseSolverNames), solver)
{
    LinearSystem sys;
#ifdef HAVE_SUITESPARSE
    auto x = EigenSparseSolve(sys.A, sys.b, solver);
    BoostUnitTest::CheckVector(x, sys.expected_x, 3);
#else
    BOOST_CHECK_THROW(EigenSparseSolve(sys.A, sys.b, solver), Exception);
#endif
}

auto mumpsSolverNames = {"MumpsLU", "MumpsLDLT"};

BOOST_DATA_TEST_CASE(mumps, bdata::make(mumpsSolverNames), solver)
{
    LinearSystem sys;
#ifdef HAVE_MUMPS
    auto x = EigenSparseSolve(sys.A, sys.b, solver);
    BoostUnitTest::CheckVector(x, sys.expected_x, 3);
#else
    BOOST_CHECK_THROW(EigenSparseSolve(sys.A, sys.b, solver), Exception);
#endif
}

auto pardisoSolverNames = {"PardisoLU", "PardisoLLT", "PardisoLDLT"};

BOOST_DATA_TEST_CASE(pardisotest, bdata::make(pardisoSolverNames), solver)
{
    LinearSystem sys;
#ifdef HAVE_MKL
    auto x = EigenSparseSolve(sys.A, sys.b, solver);
    BoostUnitTest::CheckVector(x, sys.expected_x, 3);
#else
    BOOST_CHECK_THROW(EigenSparseSolve(sys.A, sys.b, solver), Exception);
#endif
}
