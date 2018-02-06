#include "BoostUnitTest.h"
#include "mechanics/solver/MumpsSolver.h"
#include "base/Exception.h"
#include <boost/test/data/test_case.hpp>
#include <iostream>

namespace bdata = boost::unit_test::data;

using namespace NuTo;

struct LinearSystem
{
    LinearSystem()
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

        b << 3.0, -3.0, -2.0;
        expected_x << 1.0, -2.0, -2.0;
    }

    Eigen::SparseMatrix<double> A;
    Eigen::Vector3d b;
    Eigen::Vector3d expected_x;
};

auto mumpsSolverNames = {"MumpsLU", "MumpsLDLT" };

BOOST_DATA_TEST_CASE(mumpsSolvers, bdata::make(mumpsSolverNames), solver)
{
    LinearSystem sys;
    auto x = MumpsSolver(sys.A, sys.b, solver);
    BoostUnitTest::CheckVector(x, sys.expected_x, 3);
}
