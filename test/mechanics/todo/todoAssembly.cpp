#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/CellInterface.h"

BOOST_AUTO_TEST_CASE(Assembly)
{
    fakeit::Mock<NuTo::CellInterface> cell;
    NuTo::DofVector<double> dofGradient;
    NuTo::DofVector<int> dofNumbering;
    NuTo::DofType dof("mock", 1, 0);


    dofGradient[dof] = Eigen::Vector3d({42, 157, 6174});
    dofNumbering[dof] = Eigen::Vector3i({1, 0, 2});
    Method(cell, Gradient) = dofGradient;
    Method(cell, DofNumbering) = dofNumbering;

    BOOST_TEST_MESSAGE("\nTODO: implement the actual assembly class...\n");
    BOOST_CHECK(false);

    // Eigen::Vector3d expected = Eigen::Vector3d({157, 42, 6174});
    // BoostUnitTest::CheckVector(assembly.stomething..., expected, 3);
}
